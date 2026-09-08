#!/usr/bin/env python
"""Build the committed end-to-end test fixture from a real crop export.

    python scripts/prepare_e2e_fixture.py /path/to/demo_data/crop_7

The fixture is ~12 MB of real 10x-derived data checked into the repository, so
``tests/e2e/`` runs everywhere without a download. This script is what produced
it: a binary fixture nobody can regenerate is a binary nobody can review.

Three things have to happen to the copy, and the last two are the ones that get
forgotten — the same list ``scripts/prepare_demo_dataset.sh`` carries, for the
same reason:

  1. **Strip.** Drop everything derived that the viewer rebuilds anyway
     (``analysis.py``, the notebook, ``plots/``, the log) plus the per-gene
     transcript cache and the store's own trash/staging debris. 26 MB -> ~12 MB.

  2. **Neutralise the recorded paths.** The provenance graph stores absolute
     paths — the ``preamble`` node's ``data_path = Path(r"...")`` and each
     ``clustering:<key>`` node's ``read_csv`` — and ``experiment.xenium``'s
     ``palms_crop.derived_from`` names the parent crop. Committing those would
     publish the author's home directory, and every launch would then re-emit
     the preamble for a *different* path and flag all 30 descendants stale.
     They are rewritten to ``PLACEHOLDER_PATH``, which ``tests/e2e/conftest.py``
     rewrites again to wherever it stages the copy. Both directions go through
     ``rename_dataset.repair`` — the shipped tool, not a private regex.

  3. **Scrub the ARMS scan filename.** It reaches the napari layer list, and a
     real scan filename is a slide identifier. This is the leak that is
     invisible until it is in a public repository.

The safety shape matches ``prepare_demo_dataset.sh``: the destination may not be
or contain the source, the source may not be inside the destination, and the
only thing this will delete is a previous copy of the fixture.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

REPO = Path(__file__).resolve().parent.parent
DEFAULT_DEST = REPO / "tests" / "data" / "crop_7"

#: What the committed fixture records as its own location. Not a real directory
#: anywhere: it must not resolve on a developer's machine, or a launch against
#: the *source* tree would look like it worked.
PLACEHOLDER_PATH = "/palms-e2e-fixture/crop_7"

#: A neutral stand-in for the ARMS scan filename. See point 3 above.
ARMS_NAME = "demo_he_section.svs"

#: Derived outputs the viewer regenerates on launch, plus caches a test does not
#: need. `transcript_cache/` is 15 MB of per-gene feathers; without it the
#: Transcripts tab takes its parquet fallback, which is a real path worth
#: covering and costs nothing on a 683 KB file.
STRIP = (
    "analysis.py",
    "analysis_notebook.ipynb",
    "palms.log",
    "xenium_viewer.log",
    "plots",
    "transcript_cache",
)

#: Crash-safety scaffolding inside the store. `.xv_trash` alone is 2.2 MB of
#: superseded element copies, and a committed `.xv.lock` is a lock nobody holds.
STRIP_IN_STORE = (".xv_trash", ".xv_staging", ".xv_journal", ".xv.lock")


def _refuse_overlap(src: Path, dest: Path) -> None:
    """Neither directory may contain the other."""
    if dest == src or src in dest.parents:
        sys.exit(f"error: destination {dest} is, or is inside, the source")
    if dest in src.parents:
        sys.exit(f"error: destination {dest} contains the source")


def _clear_previous(dest: Path) -> None:
    """Remove a previous copy of the fixture — and nothing else."""
    if not dest.exists():
        return
    if not (dest / "experiment.xenium").exists():
        sys.exit(f"error: {dest} exists and is not a Xenium dataset — refusing to remove it")
    print(f"Removing the previous fixture at {dest}")
    shutil.rmtree(dest)


def _drop(target: Path, label: str) -> None:
    if not target.exists():
        return
    shutil.rmtree(target) if target.is_dir() else target.unlink()
    print(f"  dropped {label}")


def _strip(dest: Path) -> None:
    for name in STRIP:
        _drop(dest / name, name)


def _strip_store_debris(dest: Path) -> None:
    """Drop the crash-safety scaffolding — *after* the last write.

    ``safe_group_update`` stages, journals and trashes on every call, so the
    scrub below recreates all four. Stripping before it left a `.xv_trash`
    holding the un-scrubbed copy of the very attrs the scrub exists to remove.
    """
    for name in STRIP_IN_STORE:
        _drop(dest / "sdata_cached.zarr" / name, f"sdata_cached.zarr/{name}")


def _neutralise_paths(dest: Path) -> None:
    """Rewrite every recorded absolute path to sit under PLACEHOLDER_PATH.

    Substitution is on the **parent** directory, not the dataset's own. A crop
    export records ``output_dir=Path('<parent>')`` — the directory it wrote the
    crop *into* — which no substitution scoped to the dataset can reach, so a
    plain ``palms-rename-dataset --repair`` leaves it pointing at where the
    parent used to be. Rewriting the parent covers the dataset too, since
    substitution is path-prefix-only: ``<parent>/crop_7`` matches ``<parent>``.
    """
    from palms.scripts.rename_dataset import (
        Report, choose_graph_items, infer_old_path, read_sidecar_items,
        read_session_attrs, repair,
    )

    cache = dest / "sdata_cached.zarr"
    items = choose_graph_items(
        read_sidecar_items(dest),
        read_session_attrs(cache).get("prov_graph") or [] if cache.exists() else [],
    )
    recorded = infer_old_path(items)
    if recorded is None:
        sys.exit("error: no data_path recorded in the preamble — nothing to neutralise")
    old = str(Path(recorded).parent)
    new = str(Path(PLACEHOLDER_PATH).parent)
    report = Report(old_path=Path(old), new_path=Path(new))
    repair(dest, old, new, report, dry_run=False)
    print(report.render())
    # repair() regenerates analysis.py and the notebook when they exist. They do
    # not, because _strip ran first — deliberately, so the fixture carries no
    # derived output and the tests can assert the launch recreates them.


def _neutralise_experiment_specs(dest: Path) -> None:
    """Rewrite ``palms_crop.derived_from``, which names the parent crop's path."""
    from palms.utils.zarr_safe import atomic_json

    path = dest / "experiment.xenium"
    specs = json.loads(path.read_text())
    crop = specs.get("palms_crop")
    if not isinstance(crop, dict) or "derived_from" not in crop:
        return
    chain = crop["derived_from"]
    chain = chain if isinstance(chain, list) else [chain]
    parent = Path(PLACEHOLDER_PATH).parent
    crop["derived_from"] = [str(parent / Path(p).name) for p in chain]
    atomic_json(path, specs)
    print(f"  experiment.xenium: derived_from -> {crop['derived_from']}")


def _scrub_arms_filename(dest: Path) -> str | None:
    """Replace the ARMS scan filename in the session attrs. Returns the old one."""
    from palms.utils.zarr_safe import safe_group_update

    store = dest / "sdata_cached.zarr"
    if not (store / "viewer_session").exists():
        print("  no viewer_session in the cache — nothing to scrub")
        return None
    # safe_group_update, not a hand-written zarr.json: the store is never written
    # behind the viewer's back, even one destined for tests/.
    with safe_group_update(store, "viewer_session") as (group, _staging):
        previous = group.attrs.get("arms_he_filename")
        group.attrs["arms_he_filename"] = ARMS_NAME
    # Returned so _assert_clean can look for it, but deliberately not printed:
    # it is the slide identifier this step exists to remove, and a pasted log is
    # a published log.
    print(f"  arms_he_filename {'replaced' if previous else 'set'}: {ARMS_NAME!r}")
    return previous


def _assert_clean(dest: Path, src: Path, arms_was: str | None) -> None:
    """Refuse to hand over a fixture that still names the machine it came from.

    Run here rather than printed as advice: "check it before committing" is the
    step that gets skipped, and what it is checking for is unrecoverable once
    pushed.
    """
    needles = {str(Path.home()), str(src), str(src.parent)}
    if arms_was:
        needles.add(str(arms_was).lstrip("_").split(".")[0])
    needles = {n.encode() for n in needles if n and n != "/"}

    guilty: dict[str, list[str]] = {}
    for path in sorted(dest.rglob("*")):
        if not path.is_file():
            continue
        blob = path.read_bytes()
        for needle in needles:
            if needle in blob:
                guilty.setdefault(needle.decode(), []).append(
                    str(path.relative_to(dest)))
    if guilty:
        lines = [f"  {needle!r} in: {', '.join(files[:4])}"
                 for needle, files in guilty.items()]
        sys.exit("error: the fixture still carries host-identifying strings:\n"
                 + "\n".join(lines))
    print("  no host paths, no scan identifier")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source", type=Path, help="the crop export to build the fixture from")
    parser.add_argument("dest", type=Path, nargs="?", default=DEFAULT_DEST,
                        help=f"where to write it (default: {DEFAULT_DEST})")
    args = parser.parse_args(argv)

    src = args.source.expanduser().resolve()
    dest = args.dest.expanduser().resolve()
    if not (src / "experiment.xenium").exists():
        sys.exit(f"error: {src} is not a Xenium dataset (no experiment.xenium)")
    _refuse_overlap(src, dest)
    _clear_previous(dest)

    print(f"Copying {src} -> {dest}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest, symlinks=True)

    print("Stripping derived outputs and caches")
    _strip(dest)
    print("Neutralising recorded paths")
    _neutralise_paths(dest)
    _neutralise_experiment_specs(dest)
    print("Scrubbing the ARMS scan filename")
    arms_was = _scrub_arms_filename(dest)
    print("Dropping the store's crash-safety scaffolding")
    _strip_store_debris(dest)
    print("Checking the result")
    _assert_clean(dest, src, arms_was)

    size = sum(f.stat().st_size for f in dest.rglob("*") if f.is_file())
    count = sum(1 for f in dest.rglob("*") if f.is_file())
    print(f"\nReady: {dest} ({size / 1e6:.1f} MB in {count} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
