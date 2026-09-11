#!/usr/bin/env python3
"""Point a dataset's cached H&E back at the file it was loaded from.

The viewer records two things about an H&E that are not in the image itself:
``he_path``, the file it came from, and ``he_pixel_size_um``, the physical pixel
size its OME metadata declared. Both live in the store's ``viewer_session``
attrs, and both used to be lost on the first session restore --
``tab_he_registration._restore_session`` did not pass them through, and because
the H&E attrs are written as "computed value wins, ``None`` included" (so that
clearing an H&E actually clears it), one restore erased ``he_path`` from the
store for good.

What that costs is not cosmetic. With no path on record, ``he:load`` falls back
to its ``from_store`` block -- the copy in the viewer's cache -- so the recorded
notebook cell replays only against that cache, and with no pixel size Coarse
Align loses its best scale prior and falls back to a tissue-area estimate.

This puts both back, for datasets registered before the fix. It matches each
store's recorded ``he_filename`` against a directory of H&E images and
**requires the candidate's (height, width) to equal the recorded
``he_shape_yx``** before it writes anything: a filename match alone is not
evidence that it is the same image, and the shape is. Nothing is read but TIFF
metadata -- no pixels, no SpatialData load.

Usage::

    palms-relink-he /data/output-XETG... --he-dir /slides         # one dataset
    palms-relink-he /data/output-XETG*/ --he-dir /slides          # shell glob
    palms-relink-he /data --recursive --he-dir /slides            # walk a root
    palms-relink-he ... --dry-run                                 # write nothing

Close the viewer first: it writes its own session back at exit.

The ARMS overlay's ``arms_he_path`` is deliberately out of scope --
``tab_arms._restore_session`` always passed its path through, so it never had
this bug.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from palms.utils.zarr_safe import safe_group_update, store_lock

#: Image suffixes the viewer's own "Load H&E Image..." dialog offers.
HE_SUFFIXES = (".ome.tif", ".ome.tiff", ".tif", ".tiff", ".svs")


@dataclass
class Result:
    """What this tool decided about one dataset, and why."""

    dataset: Path
    #: relinked | linked | no-session | no-he | not-found | shape-mismatch | error
    status: str
    detail: str = ""
    he_path: str | None = None
    px_um: float | None = None

    @property
    def changed(self) -> bool:
        return self.status == "relinked"


def find_datasets(roots, recursive: bool) -> list[Path]:
    """The dataset directories to consider, in a stable order.

    A dataset is a directory holding ``sdata_cached.zarr``; that, rather than
    ``experiment.xenium``, is the test, because the attrs being repaired live in
    the cache and a crop export has no raw Xenium output at all.
    """
    found: list[Path] = []
    for root in roots:
        root = Path(root).resolve()
        if (root / "sdata_cached.zarr").is_dir():
            found.append(root)
        elif recursive and root.is_dir():
            # Two explicit levels rather than `rglob`: a run directory holding
            # `output-XETG.../`, or a directory of those. An unbounded walk here
            # would descend into every zarr store it finds on the way, which is
            # thousands of directories per dataset on a network share.
            for pattern in ("*/sdata_cached.zarr", "*/*/sdata_cached.zarr"):
                found.extend(sorted(p.parent for p in root.glob(pattern) if p.is_dir()))
    seen, unique = set(), []
    for path in found:
        if path not in seen:
            seen.add(path)
            unique.append(path)
    return unique


def index_he_files(he_dir: Path) -> dict[str, Path]:
    """Candidate H&E files by basename, non-recursively.

    Keyed by name alone: the store records ``he_filename``, which is exactly
    ``Path(path).name``, so that is the only thing there is to match on. A
    directory with two files of the same name cannot arise.
    """
    index: dict[str, Path] = {}
    for path in sorted(he_dir.iterdir()):
        if path.is_file() and path.name.lower().endswith(HE_SUFFIXES):
            index.setdefault(path.name, path)
    return index


def he_file_facts(path: Path):
    """``(height, width, pixel_size_um_or_None)`` for an H&E file, from metadata.

    Reads the series shape and the OME pixel size and nothing else, via the same
    ``registration.he_pixel_size_um`` the viewer records, so a relinked dataset
    gets the value a fresh load would have given it.
    """
    import tifffile

    from palms.utils.registration import he_pixel_size_um

    with tifffile.TiffFile(str(path)) as tif:
        shape = tuple(tif.series[0].shape)
        axes = (tif.series[0].axes or "").upper()
        px_um = he_pixel_size_um(tif)
    # (Y, X, S) for interleaved RGB, (C, Y, X) for planar; the two largest dims
    # are spatial either way, and `axes` settles it when it is present.
    if len(shape) == 2:
        height, width = shape
    elif "Y" in axes and "X" in axes and len(axes) == len(shape):
        height, width = shape[axes.index("Y")], shape[axes.index("X")]
    else:
        height, width = sorted(shape)[-2:]
    return int(height), int(width), px_um


def plan_dataset(data_path: Path, he_index: dict[str, Path], relink_all: bool) -> Result:
    """Decide what, if anything, this dataset needs — without writing."""
    from palms.utils.session import _read_prev_attrs

    cache_path = data_path / "sdata_cached.zarr"
    attrs = _read_prev_attrs(cache_path)
    if not attrs:
        return Result(data_path, "no-session", "no viewer_session in the cache")

    he_filename = attrs.get("he_filename")
    if not he_filename:
        return Result(data_path, "no-he", "no H&E has been loaded for this dataset")

    recorded = attrs.get("he_path")
    has_px = attrs.get("he_pixel_size_um") is not None
    if not relink_all and recorded and Path(recorded).exists() and has_px:
        return Result(data_path, "linked", f"{he_filename} -> {recorded}")

    candidate = he_index.get(he_filename)
    if candidate is None:
        return Result(data_path, "not-found",
                      f"no {he_filename} in the H&E directory")

    shape_yx = attrs.get("he_shape_yx")
    height, width, px_um = he_file_facts(candidate)
    if shape_yx is not None and tuple(int(v) for v in shape_yx) != (height, width):
        # The name matched and the image did not. Refusing is the whole point:
        # writing this path would silently point the dataset at a different
        # section, and every transform it carries is expressed in the old one's
        # pixels.
        return Result(data_path, "shape-mismatch",
                      f"{he_filename} is {height}x{width}, the store recorded "
                      f"{shape_yx[0]}x{shape_yx[1]}")

    return Result(data_path, "relinked", f"{he_filename} -> {candidate}",
                  he_path=str(candidate), px_um=px_um)


def apply_result(result: Result) -> None:
    """Write the repaired attrs, through the store's own safe update path."""
    cache_path = result.dataset / "sdata_cached.zarr"
    with store_lock(cache_path), \
            safe_group_update(cache_path, "viewer_session") as (session, _stage):
        session.attrs["he_path"] = result.he_path
        if result.px_um is not None:
            session.attrs["he_pixel_size_um"] = float(result.px_um)


# ── CLI ──────────────────────────────────────────────────────────────────────

_STATUS_ORDER = ("relinked", "shape-mismatch", "not-found", "error",
                 "linked", "no-he", "no-session")


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="palms-relink-he",
        description=(
            "Restore the recorded H&E file path and pixel size for datasets "
            "whose session lost them, by matching the recorded filename against "
            "a directory of H&E images."
        ),
        epilog=(
            "A candidate is used only when its height and width equal the shape "
            "the store recorded. Close the viewer before running this."
        ),
    )
    parser.add_argument("datasets", type=Path, nargs="+",
                        help="dataset directories (those holding sdata_cached.zarr)")
    parser.add_argument("--he-dir", type=Path, required=True,
                        help="directory holding the H&E images")
    parser.add_argument("--recursive", action="store_true",
                        help="treat each argument as a root to search for datasets")
    parser.add_argument("--relink-all", action="store_true",
                        help="also rewrite datasets whose recorded path still resolves")
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would change and write nothing")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = _parse_args(argv)

    he_dir = args.he_dir.resolve()
    if not he_dir.is_dir():
        print(f"error: no such directory: {he_dir}", file=sys.stderr)
        return 2
    he_index = index_he_files(he_dir)
    if not he_index:
        print(f"error: no H&E images in {he_dir}", file=sys.stderr)
        return 2

    datasets = find_datasets(args.datasets, args.recursive)
    if not datasets:
        print("error: no datasets found (looked for sdata_cached.zarr)", file=sys.stderr)
        return 2

    print(f"{len(he_index)} H&E file(s) in {he_dir}")
    print(f"{len(datasets)} dataset(s) to check\n")

    results = []
    for data_path in datasets:
        try:
            result = plan_dataset(data_path, he_index, args.relink_all)
        except Exception as exc:                       # noqa: BLE001 - reported, not raised
            # Named apart from "not-found": a store this tool could not read is
            # a different thing from one whose H&E file is absent, and collapsing
            # them would read as a missing slide when it is a broken cache.
            result = Result(data_path, "error", f"could not be read: {exc}")
        results.append(result)
        print(f"  [{result.status:>14}] {data_path.name}: {result.detail}")

    todo = [r for r in results if r.changed]
    print()
    for status in _STATUS_ORDER:
        n = sum(1 for r in results if r.status == status)
        if n:
            print(f"  {status:>14}: {n}")

    if not todo:
        print("\nNothing to do.")
        return 0
    if args.dry_run:
        print(f"\n--dry-run: {len(todo)} dataset(s) would be relinked, nothing written.")
        return 0

    failures = 0
    for result in todo:
        try:
            apply_result(result)
            px = f", pixel size {result.px_um} um" if result.px_um is not None else ""
            print(f"  relinked {result.dataset.name}{px}")
        except Exception as exc:                       # noqa: BLE001 - per dataset
            failures += 1
            print(f"  FAILED {result.dataset.name}: {exc}", file=sys.stderr)

    print(f"\nRelinked {len(todo) - failures} of {len(todo)} dataset(s).")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
