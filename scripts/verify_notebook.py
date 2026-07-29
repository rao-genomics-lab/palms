#!/usr/bin/env python
"""Replay a recorded session's notebook against the raw data and measure it.

The viewer's reproducibility claim has two halves. The first — that the code it
records is the code it ran — is structural, enforced by ``utils/steps.py`` and
asserted in the test suite. The second cannot be asserted, only run: does
executing the exported notebook, from the raw Xenium output, in a clean kernel,
reproduce what the viewer showed? This script runs it.

    python scripts/verify_notebook.py /path/to/xenium/output/ --out report.json

What it does:

1. Reads the provenance graph out of ``<cache>/viewer_session`` attrs. No GUI,
   no napari, no SpatialData load — just zarr attrs.
2. Derives the notebook with ``notebook_export.write_graph_notebook`` and
   appends one *injected* cell that dumps the replayed results to disk. That
   cell is clearly marked and is the only thing added.
3. Executes the notebook in a fresh kernel with ``allow_errors=False``, timing
   every cell.
4. Compares the replayed ``adata.obs`` against the clusterings the viewer
   persisted in the zarr table (``clustering_*``), and the replayed ranked genes
   against ``uns['rank_genes_groups']``.
5. Writes a JSON report: per-clustering ARI and cluster counts, top-N gene
   agreement, per-cell wall-clock, package versions — and the ids of every
   **comment-only node**, which replay as silent no-ops. That list is the
   remaining recording work, enumerated by measurement rather than by reading
   code.

Exit status is 0 when every comparison passed, 1 when any failed, and 2 when the
notebook could not be executed at all.
"""
from __future__ import annotations

import argparse
import ast
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from xenium_viewer.utils import notebook_export  # noqa: E402
from xenium_viewer.utils.prov_graph import ProvGraph  # noqa: E402

CACHE_DIRNAME = "sdata_cached.zarr"
SESSION_GROUP = "viewer_session"
SIDECAR_DIRNAME = "viewer_cache"    # utils/adata_persistence.py
PROV_GRAPH_SIDECAR = "prov_graph.json"
TABLE_PATH = "tables/table"
CLUSTERING_PREFIX = "clustering_"   # utils/adata_persistence.py
DEFAULT_TOP_N = 10

# Appended to the exported notebook — the only cell this script adds, and the
# only way a separate kernel can hand its results back. Written defensively:
# which names a real session's graph binds depends on which analyses were run.
_DUMP_CELL = '''\
# ── injected by scripts/verify_notebook.py (not part of the recorded analysis) ──
import numpy as _np
from pathlib import Path as _Path

_verify_out = _Path({out!r})
# Everything is compared as strings, and an all-string frame is the one thing
# parquet will always accept out of a real obs table.
adata.obs.astype(str).to_parquet(_verify_out / "replay_obs.parquet")
if "rank_df" in globals():
    rank_df.to_parquet(_verify_out / "replay_rank.parquet")
if "nhood_zscore" in globals():
    _np.save(_verify_out / "replay_nhood.npy", _np.asarray(nhood_zscore))
'''


# ── reading the recorded session ─────────────────────────────────────────────

def cache_path(data_path: Path, explicit: Path | None = None) -> Path:
    path = explicit if explicit is not None else data_path / CACHE_DIRNAME
    if not path.exists():
        raise SystemExit(f"no viewer cache at {path}")
    return path


def read_graph(data_path: Path, cache: Path) -> tuple[ProvGraph, str]:
    """Load the provenance graph, sidecar first, session attr second.

    The viewer rewrites ``viewer_cache/prov_graph.json`` on every recorded step;
    the zarr attr is only written on a dataset switch or at exit. Reading the
    attr alone measured a graph 16 minutes behind the results sitting in the
    same store — the replay then "failed" to produce clusterings the recorded
    session never contained. Returns the graph and which source it came from.
    """
    sidecar = data_path / SIDECAR_DIRNAME / PROV_GRAPH_SIDECAR
    if sidecar.exists():
        try:
            items = json.loads(sidecar.read_text())
        except (OSError, ValueError) as exc:
            print(f"warning: {sidecar} unreadable ({exc}); falling back to the "
                  f"session attrs")
        else:
            if items:
                return ProvGraph.from_list(list(items)), str(sidecar)

    import zarr
    store = zarr.open(str(cache), mode="r")
    if SESSION_GROUP not in store:
        raise SystemExit(f"{cache} has no {SESSION_GROUP}/ — nothing was recorded")
    items = dict(store[SESSION_GROUP].attrs).get("prov_graph")
    if not items:
        raise SystemExit(
            f"{cache}/{SESSION_GROUP} records no provenance graph. Run an "
            f"analysis in the viewer and save the session first."
        )
    return ProvGraph.from_list(list(items)), f"{cache}/{SESSION_GROUP} attrs"


def read_viewer_obs(cache: Path):
    """The viewer's persisted ``obs`` table, without loading X."""
    import zarr
    from anndata.io import read_elem
    store = zarr.open(str(cache), mode="r")
    return read_elem(store[f"{TABLE_PATH}/obs"])


def read_viewer_rank_genes(cache: Path) -> tuple[dict, str | None]:
    """Top gene names per group from the viewer's persisted ``uns``."""
    import zarr
    from anndata.io import read_elem
    store = zarr.open(str(cache), mode="r")
    uns_group = store[f"{TABLE_PATH}/uns"]
    if "rank_genes_groups" not in uns_group:
        return {}, None
    rgg = read_elem(uns_group["rank_genes_groups"])
    groupby = read_elem(uns_group["rank_genes_groupby"]) \
        if "rank_genes_groupby" in uns_group else None
    names = rgg.get("names")
    if names is None:
        return {}, groupby
    # anndata stores `names` as a structured array (one field per group).
    if hasattr(names, "dtype") and names.dtype.names:
        return {g: list(names[g]) for g in names.dtype.names}, groupby
    return {str(g): list(v) for g, v in dict(names).items()}, groupby


def comment_only_nodes(graph: ProvGraph) -> list[str]:
    """Node ids that claim to be code but parse to no executable statement.

    These are Phase 0.3's punch list. A notebook containing one runs it
    successfully and does nothing — the step it claims to document is simply
    absent from the replay, and no error is raised anywhere.

    ``NOTE`` nodes are excluded: they are viewer state that has no code
    equivalent (the canvas background, an overlay), declared as such and
    rendered as markdown. Counting them here buried the real gaps among them.
    """
    from xenium_viewer.utils.prov_graph import NOTE

    out = []
    for node_id in graph.topo_sort():
        node = graph.get(node_id)
        if node.kind == NOTE:
            continue
        try:
            tree = ast.parse(node.code)
        except SyntaxError:
            out.append(node_id)
            continue
        if not tree.body:
            out.append(node_id)
    return out


def note_nodes(graph: ProvGraph) -> list[str]:
    """Node ids declared as viewer state — no code, and none expected."""
    from xenium_viewer.utils.prov_graph import NOTE

    return [nid for nid in graph.topo_sort() if graph.get(nid).kind == NOTE]


# ── the replay ───────────────────────────────────────────────────────────────

def build_notebook(graph: ProvGraph, work_dir: Path) -> tuple[Path, list]:
    """Write the notebook; also return each cell's originating node id.

    The node ids are what makes a failure actionable: nbclient reports a cell
    index, and "cell 4 failed" says nothing about which recorded step is broken.
    """
    from xenium_viewer.utils.prov_graph import graph_to_cells
    derived = graph_to_cells(graph)
    node_ids = [cell.node_id for cell in derived]
    node_ids.append(None)  # the injected dump cell belongs to no node
    cells = [(cell.cell_type, cell.source) for cell in derived]
    cells.append(("code", _DUMP_CELL.format(out=str(work_dir))))
    nb_path = work_dir / "verify_notebook.ipynb"
    notebook_export.write_notebook(cells, nb_path)
    return nb_path, node_ids


def execute(nb_path: Path, data_path: Path, timeout: int,
            timings: list[dict], node_ids: list, cursor: dict) -> float:
    """Execute the notebook, appending per-cell timings; return wall-clock.

    Runs with the *dataset* as the working directory, matching the environment
    a user would replay it in. *timings* is passed in rather than returned so
    that a failing run still reports how far it got and how long each completed
    cell took — a replay that dies on cell 4 of 9 is exactly when that detail
    matters most. *cursor* is filled in with the failing cell, so the caller can
    name the node that broke rather than only quoting nbclient's traceback.

    ``on_cell_executed`` fires for the failing cell too, immediately before
    nbclient raises, so the failure is taken from ``on_cell_error`` and the
    timing entry it already appended is flagged rather than trusted.
    """
    started: dict[int, float] = {}

    def node_of(cell_index):
        return node_ids[cell_index] if cell_index < len(node_ids) else None

    def on_start(cell, cell_index, **_):
        started[cell_index] = time.perf_counter()

    def on_executed(cell, cell_index, **_):
        begin = started.pop(cell_index, None)
        if begin is None:
            return
        source = (cell.source or "").strip().splitlines()
        timings.append({
            "cell": cell_index,
            "node": node_of(cell_index),
            "seconds": round(time.perf_counter() - begin, 3),
            "first_line": source[0][:100] if source else "",
        })

    def on_error(cell, cell_index, **_):
        cursor["cell"] = cell_index
        cursor["node"] = node_of(cell_index)
        for entry in timings:
            if entry["cell"] == cell_index:
                entry["failed"] = True

    t0 = time.perf_counter()
    notebook_export.execute_notebook(
        nb_path, cwd=data_path, timeout=timeout,
        on_cell_start=on_start, on_cell_executed=on_executed,
        on_cell_error=on_error,
    )
    return round(time.perf_counter() - t0, 3)


# ── comparison ───────────────────────────────────────────────────────────────

def _align(viewer_obs, replay_obs):
    """Row-align the two obs tables, preferring an explicit ``cell_id``.

    The viewer's table index and the notebook's may differ (the viewer may have
    been launched on a cropped or re-segmented export), so alignment is on cell
    identity, not on row order.
    """
    import pandas as pd

    def keyed(frame):
        if "cell_id" in frame.columns:
            return frame.set_index(pd.Index(frame["cell_id"].astype(str)))
        return frame.set_index(frame.index.astype(str))

    left, right = keyed(viewer_obs), keyed(replay_obs)
    shared = left.index.intersection(right.index)
    return left.loc[shared], right.loc[shared], len(shared)


def compare_clusterings(viewer_obs, replay_obs) -> list[dict]:
    from sklearn.metrics import adjusted_rand_score

    left, right, n_shared = _align(viewer_obs, replay_obs)
    results = []
    for column in viewer_obs.columns:
        if not column.startswith(CLUSTERING_PREFIX):
            continue
        key = column[len(CLUSTERING_PREFIX):]
        if key not in right.columns:
            results.append({
                "clustering": key, "status": "not_in_replay",
                "note": "the viewer persisted this clustering but the notebook "
                        "did not produce it",
            })
            continue
        viewer_labels = left[column].astype(str)
        replay_labels = right[key].astype(str)
        mask = (viewer_labels != "nan") & (replay_labels != "nan")
        viewer_labels, replay_labels = viewer_labels[mask], replay_labels[mask]
        if len(viewer_labels) == 0:
            results.append({"clustering": key, "status": "empty"})
            continue
        ari = float(adjusted_rand_score(viewer_labels, replay_labels))
        results.append({
            "clustering": key,
            "status": "ok" if ari == 1.0 else "diverged",
            "ari": ari,
            "identical_labels": bool((viewer_labels.values == replay_labels.values).all()),
            "n_cells_compared": int(len(viewer_labels)),
            "n_cells_shared": int(n_shared),
            "n_clusters_viewer": int(viewer_labels.nunique()),
            "n_clusters_replay": int(replay_labels.nunique()),
        })
    return results


def compare_rank_genes(viewer_names: dict, replay_rank, top_n: int) -> dict:
    """Top-N gene name agreement per group."""
    if not viewer_names:
        return {"status": "no_viewer_result"}
    if replay_rank is None:
        return {"status": "not_in_replay"}

    replayed = {
        str(group): list(sub.head(top_n)["names"])
        for group, sub in replay_rank.groupby("group", observed=True)
    }
    per_group = {}
    for group, names in viewer_names.items():
        wanted = [str(n) for n in list(names)[:top_n]]
        got = [str(n) for n in replayed.get(str(group), [])][:top_n]
        per_group[str(group)] = {
            "identical": wanted == got,
            "n_shared": len(set(wanted) & set(got)),
            "n_compared": len(wanted),
            "viewer_top": wanted,
            "replay_top": got,
        }
    all_identical = all(g["identical"] for g in per_group.values())
    return {
        "status": "ok" if all_identical else "diverged",
        "top_n": top_n,
        "groups": per_group,
    }


def package_versions() -> dict:
    """The same pins the ``environment`` node records, plus the replay's own."""
    from xenium_viewer.utils.environment import RECORDED_PACKAGES, package_versions
    return package_versions(RECORDED_PACKAGES + ("nbclient", "nbformat"))


# ── entry point ──────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Replay a recorded session's notebook and compare it to the "
                    "viewer's own results.",
    )
    parser.add_argument("data_path", type=Path, help="Xenium output directory")
    parser.add_argument("--out", type=Path, default=Path("verify_report.json"),
                        help="where to write the JSON report")
    parser.add_argument("--cache", type=Path, default=None,
                        help="zarr cache (default: <data_path>/sdata_cached.zarr)")
    parser.add_argument("--work-dir", type=Path, default=None,
                        help="keep the notebook and dumped results here "
                             "(default: a temp dir, deleted afterwards)")
    parser.add_argument("--timeout", type=int, default=7200,
                        help="per-cell execution timeout in seconds")
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N,
                        help="how many ranked genes per group to compare")
    parser.add_argument("--dry-run", action="store_true",
                        help="build the notebook and list the comment-only nodes, "
                             "but do not execute it (seconds instead of an hour)")
    args = parser.parse_args(argv)

    data_path = args.data_path.resolve()
    cache = cache_path(data_path, args.cache)
    graph, graph_source = read_graph(data_path, args.cache or cache)
    skipped = comment_only_nodes(graph)
    notes = note_nodes(graph)

    print(f"Provenance graph: {len(graph)} nodes from {graph_source} "
          f"({len(skipped)} comment-only, which replay as no-ops; "
          f"{len(notes)} viewer-state notes, which are not code by design)")

    if args.work_dir is not None:
        work_dir = args.work_dir.resolve()
        work_dir.mkdir(parents=True, exist_ok=True)
        cleanup = None
    else:
        cleanup = tempfile.mkdtemp(prefix="xv-verify-")
        work_dir = Path(cleanup)

    report = {
        "data_path": str(data_path),
        "cache": str(cache),
        "graph_source": graph_source,
        "n_nodes": len(graph),
        "node_ids": graph.topo_sort(),
        "comment_only_nodes": skipped,
        "note_nodes": notes,
        "versions": package_versions(),
    }

    try:
        nb_path, node_ids = build_notebook(graph, work_dir)
        print(f"Notebook: {nb_path}")

        if args.dry_run:
            report["execution"] = {"status": "skipped (--dry-run)"}
            report["passed"] = None
            args.out.write_text(json.dumps(report, indent=2, default=str))
            print(f"\nComment-only nodes that would be skipped: {len(skipped)}")
            for node_id in skipped:
                print(f"    {node_id}")
            if cleanup is not None:
                print("\n(the notebook above is in a temp dir and is about to be "
                      "deleted — pass --work-dir to keep it)")
            print(f"\nReport: {args.out}")
            return 0

        print("Executing in a clean kernel — this replays the whole analysis "
              "from the raw output…")
        timings: list[dict] = []
        cursor: dict = {}
        try:
            total = execute(nb_path, data_path, args.timeout, timings, node_ids,
                            cursor)
        except Exception as exc:           # noqa: BLE001 — reported, not raised
            report["execution"] = {
                "status": "failed",
                "failed_node": cursor.get("node"),
                "failed_cell": cursor.get("cell"),
                "cells_completed": len([c for c in timings if not c.get("failed")]),
                "cells": timings,
                "error": f"{type(exc).__name__}: {exc}",
            }
            args.out.write_text(json.dumps(report, indent=2, default=str))
            print(f"\nFAILED at node {cursor.get('node')!r} "
                  f"(cell {cursor.get('cell')}), after "
                  f"{report['execution']['cells_completed']} cell(s) ran:")
            print(f"  {type(exc).__name__}: {str(exc).splitlines()[-1][:300]}")
            if cleanup is not None:
                print("Re-run with --work-dir to keep the failing notebook.")
            print(f"Report: {args.out}")
            return 2

        report["execution"] = {
            "status": "ok",
            "total_seconds": total,
            "cells": timings,
            "slowest": sorted(timings, key=lambda c: -c["seconds"])[:5],
        }

        import pandas as pd
        replay_obs = pd.read_parquet(work_dir / "replay_obs.parquet")
        rank_path = work_dir / "replay_rank.parquet"
        replay_rank = pd.read_parquet(rank_path) if rank_path.exists() else None

        viewer_obs = read_viewer_obs(cache)
        viewer_names, groupby = read_viewer_rank_genes(cache)

        report["clusterings"] = compare_clusterings(viewer_obs, replay_obs)
        report["rank_genes"] = compare_rank_genes(viewer_names, replay_rank, args.top_n)
        report["rank_genes"]["groupby"] = groupby
    finally:
        if cleanup is not None:
            shutil.rmtree(cleanup, ignore_errors=True)

    failures = [c for c in report["clusterings"] if c.get("status") != "ok"]
    if report["rank_genes"].get("status") == "diverged":
        failures.append({"rank_genes": "diverged"})
    report["passed"] = not failures

    args.out.write_text(json.dumps(report, indent=2, default=str))
    _print_summary(report, args.out)
    return 0 if report["passed"] else 1


def _print_summary(report: dict, out_path: Path) -> None:
    print()
    print(f"Replayed in {report['execution']['total_seconds']}s")
    for entry in report["clusterings"]:
        if entry.get("status") == "ok":
            print(f"  ✓ {entry['clustering']}: ARI = {entry['ari']:.6f} "
                  f"({entry['n_clusters_viewer']} clusters, "
                  f"{entry['n_cells_compared']} cells)")
        else:
            print(f"  ✗ {entry['clustering']}: {entry.get('status')} "
                  f"{entry.get('ari', '')}")
    rank = report["rank_genes"]
    if rank.get("status") == "ok":
        print(f"  ✓ rank genes: top-{rank['top_n']} identical in all "
              f"{len(rank['groups'])} groups")
    elif rank.get("status") == "diverged":
        bad = [g for g, v in rank["groups"].items() if not v["identical"]]
        print(f"  ✗ rank genes: {len(bad)} group(s) differ: {bad}")
    else:
        print(f"  – rank genes: {rank.get('status')}")

    if any(e.get("status") == "not_in_replay" for e in report["clusterings"]):
        print("\n  Results with no node behind them usually mean the graph on disk "
              "predates them.\n  It is written on every recorded step to "
              "viewer_cache/prov_graph.json; a session\n  recorded by an older "
              "build only reached disk at viewer exit, so close the viewer "
              "and re-run.")

    skipped = report["comment_only_nodes"]
    print(f"\nComment-only nodes skipped on replay: {len(skipped)}")
    for node_id in skipped:
        print(f"    {node_id}")
    print(f"\nReport: {out_path}")


if __name__ == "__main__":
    sys.exit(main())
