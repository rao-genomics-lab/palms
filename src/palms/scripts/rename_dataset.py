#!/usr/bin/env python3
"""Rename or move a Xenium dataset directory without breaking its provenance.

Nothing in this codebase stores a *dataset name*: not the AnnData table, not the
SpatialData store, not the cache manifest. The folder name reaches only the
napari window title. Cache freshness is a content hash of ``experiment.xenium``
(``loader._source_fingerprint``), so moving a dataset never triggers a rebuild,
and every cache/sidecar location is re-derived from ``data_path`` at launch.

What a move *does* break is the handful of **absolute paths recorded inside the
provenance graph and the session attrs** — chiefly the ``preamble`` node's
``data_path = Path(r"...")`` and the ``clustering:<key>`` nodes' ``read_csv`` of
``analysis/clustering/<key>/clusters.csv``. Those are the notebook's source of
truth, so after a hand-rolled ``mv`` the exported notebook replays against a
path that no longer exists.

Worse, it does not stay quietly broken. ``app.py`` re-emits the preamble for the
current ``data_path`` on every launch; after a move that is an ``upsert`` with
changed code, and ``ProvGraph.upsert`` flags every transitive descendant stale.
The first launch after a manual rename therefore marks the *entire* notebook ⚠
even though nothing was recomputed. Repairing the graph before that launch is
the reason this tool exists.

Usage::

    palms-rename-dataset /data/old_name new_name      # rename in place
    palms-rename-dataset /data/old_name /other/place  # move
    palms-rename-dataset /data/new_name --repair      # fix an already-moved one
    palms-rename-dataset ... --dry-run                # report, write nothing

Close the viewer first. There is no way to detect an open viewer: the store's
flock is held only for the duration of a write, not for the session.
"""

from __future__ import annotations

import argparse
import errno
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from palms.utils.adata_persistence import sidecar_dir  # noqa: E402
from palms.utils.zarr_safe import atomic_json, safe_group_update, store_lock  # noqa: E402

#: Session attrs holding a filesystem path. All four are frequently ``None`` in
#: practice — ``_on_he_restored_from_sdata`` does not repopulate
#: ``he_state["he_path"]``, so it decays to None on the first re-save — but a
#: dataset saved in the same session that loaded the image does carry them.
PATH_ATTRS = ("he_path", "arms_he_path", "arms_geojson_path", "arms_csv_path")

#: Sidecar JSONs under ``viewer_cache/`` whose string values may hold paths.
CNV_SIDECAR_GLOBS = ("cnv_*_result.json", "cnv_*_params.json")

CODE_FILE = "analysis.py"
NOTEBOOK_FILE = "analysis_notebook.ipynb"


# ── path substitution ────────────────────────────────────────────────────────

def _prefix_pattern(old: str) -> "re.Pattern[str]":
    """Match *old* only where it is a whole path or a path prefix.

    The lookahead is what keeps ``/data/foo`` from rewriting ``/data/foobar``
    and ``/data/foo.bak``: a match must be followed by a separator, or by
    something that cannot continue a path component (a quote, whitespace, a
    comma, end of string).
    """
    return re.compile(re.escape(old) + r"(?=" + re.escape(os.sep) + r"|$|[^A-Za-z0-9_.\-])")


def sub_text(text: str, old: str, new: str) -> tuple[str, int]:
    """Rewrite every path-prefix occurrence of *old* inside *text*.

    Used for recorded code, where the path is embedded in a source line rather
    than being the whole value.
    """
    return _prefix_pattern(old).subn(new.replace("\\", "\\\\"), text)


def sub_value(value, old: str, new: str) -> tuple[object, int]:
    """Rewrite *value* if it is a string under *old*; recurse into dicts/lists.

    Anything that is not a string is returned unchanged, so a path that lives
    outside the dataset directory (an H&E on another volume, say) survives
    untouched — that file did not move.
    """
    if isinstance(value, str):
        return sub_text(value, old, new)
    if isinstance(value, dict):
        total = 0
        out = {}
        for key, item in value.items():
            out[key], n = sub_value(item, old, new)
            total += n
        return out, total
    if isinstance(value, list):
        total = 0
        out = []
        for item in value:
            fixed, n = sub_value(item, old, new)
            out.append(fixed)
            total += n
        return out, total
    return value, 0


# ── report ───────────────────────────────────────────────────────────────────

@dataclass
class Change:
    target: str
    count: int
    note: str = ""


@dataclass
class Report:
    old_path: Path
    new_path: Path
    moved: bool = False
    dry_run: bool = False
    changes: list[Change] = field(default_factory=list)

    def add(self, target: str, count: int, note: str = "") -> None:
        self.changes.append(Change(target, count, note))

    @property
    def total(self) -> int:
        """Path references rewritten. Regenerated files (count -1) are not one."""
        return sum(c.count for c in self.changes if c.count >= 0)

    def render(self) -> str:
        verb = "Would rewrite" if self.dry_run else "Rewrote"
        if self.moved:
            head = "Would move" if self.dry_run else "Moved"
            lines = [f"{head}: {self.old_path}", f"    ->  {self.new_path}"]
        else:
            head = "Would repair" if self.dry_run else "Repaired"
            lines = [f"{head}: {self.new_path}",
                     f"    recorded path was: {self.old_path}"]
        if not self.changes:
            lines.append("")
            lines.append("No recorded paths needed rewriting.")
            return "\n".join(lines)
        lines.append("")
        for c in self.changes:
            if c.count < 0:                      # a regenerated derived output
                lines.append(f"  {'Would regenerate' if self.dry_run else 'Regenerated'}"
                             f" {c.target} ({c.note})")
            else:
                suffix = f"  ({c.note})" if c.note else ""
                lines.append(f"  {verb} {c.count:>3} path(s) in {c.target}{suffix}")
        lines.append("")
        lines.append(f"{verb} {self.total} path reference(s) in total.")
        return "\n".join(lines)


# ── pre-flight ───────────────────────────────────────────────────────────────

class PreflightError(RuntimeError):
    """A condition that makes moving or repairing unsafe."""


def is_dataset_dir(path: Path) -> bool:
    """A raw Xenium output, or a cache-only export from the Crop Dataset tool.

    Both are worth renaming, and the second is exactly why this does not simply
    test for raw files: a Crop Dataset export's zarr *is* the data.

    ``loader.has_raw_xenium_source`` is the single definition of "raw output is
    present" (issue #17) — deliberately conservative, ``True`` unless *none* of
    its markers exists, because partial raw output is broken raw output. Reusing
    it here rather than re-testing ``cells.zarr.zip`` is the point of that fix:
    a predicate applied in one place is not a guarantee.
    """
    from palms.loader import has_raw_xenium_source

    if (path / "experiment.xenium").exists():
        return True
    return (path / "sdata_cached.zarr").exists() and not has_raw_xenium_source(path)


def copykat_is_running(data_path: Path) -> bool:
    """A detached CopyKAT worker holds absolute paths and writes to them.

    Checked by looking at the sentinel directly rather than through
    ``tab_cnv._copykat_run_state``, which would drag Qt into a CLI.
    """
    return (data_path / "plots" / "copykat_RUNNING.txt").exists()


def preflight(data_path: Path, dest: Path | None, *, dry_run: bool) -> None:
    """Refuse before touching anything, rather than half-applying."""
    from palms.utils.cache_repair import verify
    from palms.utils.zarr_safe import recover_pending

    if not data_path.is_dir():
        raise PreflightError(f"not a directory: {data_path}")
    if not is_dataset_dir(data_path):
        raise PreflightError(
            f"{data_path} does not look like a Xenium dataset "
            f"(no experiment.xenium, and no cache-only export)"
        )
    if dest is not None and dest.exists():
        raise PreflightError(f"destination already exists: {dest}")
    if copykat_is_running(data_path):
        raise PreflightError(
            "a CopyKAT run is in progress (plots/copykat_RUNNING.txt exists).\n"
            "That worker is a detached process holding absolute paths; let it "
            "finish, or remove the sentinel if you know it died."
        )

    cache_path = data_path / "sdata_cached.zarr"
    if not cache_path.exists():
        return
    report = verify(cache_path)
    if not report.exists:
        return
    if not report.readable_metadata:
        raise PreflightError(
            "the zarr cache's root metadata is unreadable; repair it first "
            "(viewer: Tools -> Cache -> Verify / Repair).\n" + report.summary()
        )
    if report.pending_ops and not dry_run:
        # Finish or unwind an interrupted safe write *before* the move, so no
        # journal is left naming the old location.
        recovered = recover_pending(cache_path)
        if recovered:
            print(f"  Recovered {len(recovered)} interrupted write(s) before moving.")
    elif report.pending_ops:
        print(f"  Note: {len(report.pending_ops)} interrupted write(s) would be "
              f"recovered before moving.")


# ── the move ─────────────────────────────────────────────────────────────────

def resolve_destination(data_path: Path, new: Path) -> Path:
    """A bare name renames in place; anything with a separator is a full path."""
    if new.parent == Path("."):
        return data_path.parent / new.name
    return new.resolve() if not new.is_absolute() else new


def move(data_path: Path, dest: Path) -> None:
    """Rename only.

    No ``shutil.move`` fallback on purpose: a cross-device fallback would copy
    a multi-gigabyte zarr non-atomically, and a half-copied store is a much
    worse outcome than a clear message telling the user to move it themselves.
    """
    try:
        os.rename(data_path, dest)
    except OSError as exc:
        if exc.errno == errno.EXDEV:
            raise PreflightError(
                f"{data_path} and {dest} are on different filesystems.\n"
                f"Copy it yourself (mv / rsync -a), then re-run:\n"
                f"    palms-rename-dataset {dest} --repair --from {data_path}"
            ) from exc
        raise


# ── repair ───────────────────────────────────────────────────────────────────

_PREAMBLE_PATH = re.compile(r'data_path\s*=\s*Path\(\s*r?"([^"]*)"\s*\)')


def infer_old_path(items: list[dict]) -> str | None:
    """The ``data_path`` the graph was last recorded with, from its preamble."""
    for node in items:
        if node.get("id") == "preamble":
            match = _PREAMBLE_PATH.search(node.get("code") or "")
            if match:
                return match.group(1)
    return None


def read_session_attrs(cache_path: Path) -> dict:
    """The ``viewer_session`` attrs, or {} when there is no session yet.

    Reuses ``session._read_prev_attrs`` rather than opening zarr here, so this
    module never touches the store outside ``safe_group_update``. It swallows
    read errors and returns {}, which is safe here only because ``preflight``
    has already refused a store with unreadable metadata.
    """
    from palms.utils.session import _read_prev_attrs

    return _read_prev_attrs(cache_path)


def read_sidecar_items(data_path: Path) -> list[dict]:
    """The provenance graph from ``viewer_cache/prov_graph.json``, or []."""
    from palms.tabs._helpers import PROV_GRAPH_SIDECAR

    path = sidecar_dir(data_path) / PROV_GRAPH_SIDECAR
    try:
        items = json.loads(path.read_text())
    except (OSError, ValueError):
        return []
    return items if isinstance(items, list) else []


def choose_graph_items(sidecar_items: list[dict], attr_items: list[dict]) -> list[dict]:
    """Which serialization to regenerate the derived files from.

    Same precedence as ``app._load_prov_graph_items``, reimplemented rather than
    imported because that module pulls in napari. The sidecar is written on
    every recorded step and the attr only at exit, so the sidecar wins unless it
    is *smaller*, which means something wrote a partial graph.
    """
    if sidecar_items and len(sidecar_items) >= len(attr_items):
        return sidecar_items
    return attr_items or sidecar_items


def repair(data_path: Path, old: str, new: str, report: Report, *, dry_run: bool) -> None:
    """Rewrite every recorded reference to *old* under *data_path*."""
    cache_path = data_path / "sdata_cached.zarr"

    sidecar_items = read_sidecar_items(data_path)
    attrs = read_session_attrs(cache_path) if cache_path.exists() else {}
    attr_items = attrs.get("prov_graph") or []

    fixed_sidecar, n_sidecar = sub_value(sidecar_items, old, new)
    fixed_attr_items, n_attr_graph = sub_value(attr_items, old, new)

    changed_attrs: dict = {}
    n_attr_paths = 0
    for key in PATH_ATTRS:
        if attrs.get(key):
            fixed, n = sub_value(attrs[key], old, new)
            if n:
                changed_attrs[key] = fixed
                n_attr_paths += n

    if n_sidecar:
        report.add("viewer_cache/prov_graph.json", n_sidecar,
                   f"{len(sidecar_items)} node(s)")
    if n_attr_graph:
        report.add("sdata_cached.zarr/viewer_session (prov_graph attr)", n_attr_graph,
                   f"{len(attr_items)} node(s)")
    if n_attr_paths:
        report.add("sdata_cached.zarr/viewer_session (path attrs)", n_attr_paths)

    # CNV sidecars. Currently write-only — tab_cnv re-derives every one of these
    # from ctx.data_path — so this is housekeeping, not a correctness fix. Do
    # not promote it into a load-bearing read path.
    cnv_fixed: list[tuple[Path, object, int]] = []
    sidecars = sidecar_dir(data_path)
    if sidecars.is_dir():
        for pattern in CNV_SIDECAR_GLOBS:
            for path in sorted(sidecars.glob(pattern)):
                try:
                    payload = json.loads(path.read_text())
                except (OSError, ValueError):
                    continue
                fixed, n = sub_value(payload, old, new)
                if n:
                    cnv_fixed.append((path, fixed, n))
                    report.add(f"viewer_cache/{path.name}", n)

    # Reported before the dry-run bail-out: a dry run that stayed silent about
    # the two files it would overwrite would be understating what it does.
    for name in (CODE_FILE, NOTEBOOK_FILE):
        if (data_path / name).exists():
            report.add(name, -1, "from the repaired graph")

    if dry_run:
        return

    # ── apply ────────────────────────────────────────────────────────────────
    with store_lock(cache_path if cache_path.exists() else None):
        if n_sidecar:
            from palms.tabs._helpers import PROV_GRAPH_SIDECAR
            atomic_json(sidecar_dir(data_path, create=True) / PROV_GRAPH_SIDECAR,
                        fixed_sidecar)
        if (n_attr_graph or n_attr_paths) and cache_path.exists():
            with safe_group_update(cache_path, "viewer_session") as (session, _stage):
                for key, value in changed_attrs.items():
                    session.attrs[key] = value
                if n_attr_graph:
                    session.attrs["prov_graph"] = fixed_attr_items
        for path, payload, _n in cnv_fixed:
            atomic_json(path, payload)

    # ── derived outputs, regenerated from the repaired graph ─────────────────
    items = choose_graph_items(
        fixed_sidecar if isinstance(fixed_sidecar, list) else [],
        fixed_attr_items if isinstance(fixed_attr_items, list) else [],
    )
    if items:
        _regenerate_derived(data_path, items)


def _atomic_replace(path: Path, write: "Callable[[Path], None]") -> None:
    """Produce *path* via a temp file and ``os.replace``, like ``atomic_json``.

    Not just crash-safety. Writing in place mutates the *inode*, so a hardlinked
    snapshot of the dataset (``cp -al``, which is how these get backed up
    cheaply) would silently receive the edit too — found by exactly that probe.
    ``os.replace`` swaps in a new inode instead, leaving any other link alone.
    """
    tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    write(tmp)
    os.replace(tmp, path)


def _regenerate_derived(data_path: Path, items: list[dict]) -> None:
    """Rewrite analysis.py and the notebook from the repaired graph.

    Both are derived outputs — the viewer rebuilds ``analysis.py`` from the graph
    on every session restore (``app.py``), with this exact recipe — so
    regenerating is more honest than patching their text, which could leave them
    disagreeing with the graph they claim to render. Only files that already
    exist are rewritten; this tool does not create outputs the user never asked
    for.
    """
    from palms.utils.prov_graph import ProvGraph, graph_to_cells

    graph = ProvGraph.from_list(items)

    code_path = data_path / CODE_FILE
    if code_path.exists():
        source = "\n".join(c.source for c in graph_to_cells(graph)
                           if c.cell_type == "code") + "\n"
        _atomic_replace(code_path, lambda p: p.write_text(source))

    nb_path = data_path / NOTEBOOK_FILE
    if nb_path.exists():
        from palms.utils.notebook_export import write_graph_notebook
        _atomic_replace(nb_path, lambda p: write_graph_notebook(graph, p))


# ── CLI ──────────────────────────────────────────────────────────────────────

def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="palms-rename-dataset",
        description=(
            "Rename or move a Xenium dataset directory, rewriting the absolute "
            "paths recorded in its provenance graph and session state."
        ),
        epilog=(
            "Close the viewer before running this. Renaming a dataset the "
            "viewer has open will not corrupt the store, but the viewer will "
            "write its own session back to the old location at exit."
        ),
    )
    parser.add_argument("dataset", type=Path,
                        help="the dataset directory (the one holding experiment.xenium)")
    parser.add_argument("new_name", type=Path, nargs="?", default=None,
                        help="new name, or a full destination path. Omit with --repair.")
    parser.add_argument("--repair", action="store_true",
                        help="do not move; only fix the paths recorded inside a "
                             "dataset that has already been moved")
    parser.add_argument("--from", dest="from_path", type=Path, default=None,
                        help="with --repair: the dataset's previous path. Inferred "
                             "from the recorded preamble when omitted.")
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would change and write nothing")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = _parse_args(argv)
    data_path = args.dataset.resolve()

    if args.repair and args.new_name is not None:
        print("error: give a new name or --repair, not both", file=sys.stderr)
        return 2
    if not args.repair and args.new_name is None:
        print("error: give a new name, or --repair to fix an already-moved "
              "dataset", file=sys.stderr)
        return 2

    try:
        if args.repair:
            return _run_repair(data_path, args)
        return _run_move(data_path, args)
    except PreflightError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _run_repair(data_path: Path, args) -> int:
    preflight(data_path, None, dry_run=args.dry_run)

    if args.from_path is not None:
        old = str(args.from_path)
    else:
        cache_path = data_path / "sdata_cached.zarr"
        items = read_sidecar_items(data_path) or (
            read_session_attrs(cache_path).get("prov_graph") or []
            if cache_path.exists() else []
        )
        old = infer_old_path(items)
        if old is None:
            raise PreflightError(
                "could not infer the previous path: no preamble node was found "
                "in this dataset's provenance graph. Pass --from explicitly."
            )
        print(f"Inferred previous path from the recorded preamble: {old}")

    new = str(data_path)
    if old == new:
        print(f"Nothing to do: the recorded path already matches {new}")
        return 0

    report = Report(old_path=Path(old), new_path=data_path, dry_run=args.dry_run)
    repair(data_path, old, new, report, dry_run=args.dry_run)
    print(report.render())
    return 0


def _run_move(data_path: Path, args) -> int:
    dest = resolve_destination(data_path, args.new_name)
    preflight(data_path, dest, dry_run=args.dry_run)

    report = Report(old_path=data_path, new_path=dest, moved=True,
                    dry_run=args.dry_run)
    if args.dry_run:
        repair(data_path, str(data_path), str(dest), report, dry_run=True)
        print(report.render())
        return 0

    move(data_path, dest)
    repair(dest, str(data_path), str(dest), report, dry_run=False)
    print(report.render())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
