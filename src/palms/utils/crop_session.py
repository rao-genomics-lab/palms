"""The `viewer_session` a cropped export needs in order to open like a dataset.

An export used to ship elements and nothing else, so the first launch found no
`viewer_session`, restored nothing, and showed none of the work that made the
dataset worth cropping. Restoring now falls back to the store's own elements, so
this group is no longer what makes an export *usable* — it is what makes it
*complete*: the filenames a layer is titled with, the UI rows for patch and
external overlays, and the provenance graph explaining where the data came from.

What is deliberately absent: the affines
----------------------------------------
No `affine_3x3`, no `arms_affine_3x3`, no `coarse_affine`. The exported element
already carries its placement, `utils/registration_seed` reads it back, and a
second copy is a second thing to disagree with — re-registering inside the export
writes the element immediately but the session only at exit, so the two would
drift within one session. One writer for geometry: the element.

The flips follow from that. The element transform is written as ``fine @ flip``,
so the flip is already inside it; the session says `False` and means "no
*further* flip", which is exactly what re-deriving one would otherwise apply
twice.

What is deliberately nulled: the source paths
---------------------------------------------
`he_path`, `arms_he_path`, `arms_geojson_path`, `arms_csv_path`. Those files
describe the *source* slide, and the export holds a slice of it. Keeping
`arms_geojson_path` is not merely untidy — `tab_arms._on_arms_restored` re-reads
it whenever the sdata tiles are empty, which would pull all of the source's tiles,
in source coordinates, into a crop that covers a fraction of it.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

log = logging.getLogger(__name__)


def _strip_affine(entries, carried_elements) -> list:
    """UI rows for overlays that actually travelled, without their affine copy."""
    out = []
    for entry in (entries or []):
        if not isinstance(entry, dict):
            continue
        name = entry.get("element_name")
        if name and name not in carried_elements:
            continue
        row = {k: v for k, v in entry.items() if k != "affine_matrix"}
        out.append(row)
    return out


def rewrite_graph_paths(items, old_path: str, new_path: str) -> list:
    """Point a carried provenance graph at the export instead of the source.

    Without this the export's first launch re-emits `preamble` for its own
    `data_path`, `upsert` sees a changed node, and every descendant is flagged
    stale — the whole notebook comes up warning about nothing.

    Substitution is prefix-only, reusing `rename_dataset.sub_value` rather than a
    second implementation: its rule that ``/d/foo`` must not match ``/d/foobar``
    is the property that keeps a path merely *near* the dataset from being
    rewritten, and it is already under test.
    """
    from palms.scripts.rename_dataset import sub_value

    rewritten, _ = sub_value(items, old_path, new_path)
    return rewritten


# Every path a recorded cell names, as a raw or plain string literal. Recorded
# code writes them as ``r"…"``; nothing recorded spans a line.
_QUOTED = re.compile(r'r?["\']([^"\'\n]+)["\']')


def recorded_paths_inside(items, export_path: str) -> dict:
    """``{node id: [path, …]}`` for every recorded path under *export_path*.

    Only paths inside the export are the export's business. One pointing
    somewhere else was never rewritten and still names the file it always did —
    an H&E slide beside the source dataset, a downloaded reference — so it is
    not this function's problem, and reporting it would be noise.
    """
    export_path = str(export_path).rstrip("/")
    found = {}
    for item in items:
        code = item.get("code") or ""
        hits = [m.group(1) for m in _QUOTED.finditer(code)
                if m.group(1) == export_path or m.group(1).startswith(export_path + "/")]
        if hits:
            found[item.get("id")] = hits
    return found


def repair_missing_reads(items, export_path: str, present_root: Path,
                         also_expected: tuple = ()) -> tuple[list, list]:
    """Fix carried nodes that now read files the export does not contain.

    ``rewrite_graph_paths`` repoints *every* recorded path at the export, which
    is right for the preamble's ``data_path`` and wrong for anything under a
    directory the export does not copy. The concrete case is 10x's own
    clustering CSVs: a parent's ``clustering:<key>`` backstop node legitimately
    reads ``<parent>/analysis/clustering/<key>/clusters.csv``, the rewrite turns
    that into ``<export>/analysis/clustering/<key>/clusters.csv``, and an export
    carries ``experiment.xenium``, ``transcripts.parquet`` and the zarr — never
    ``analysis/``. The node then points at a file that has never existed in that
    dataset, and the notebook dies there.

    Copying ``analysis/`` instead would be wrong: the crop is a cell subset, so
    the parent's CSV reindexes to NaN for every cell the crop does not contain.
    The labels the export *does* have are in its own store, so a
    ``clustering:<key>`` node is rewritten to read them back — the same cell
    ``_record_clustering`` emits when no CSV backs a key, saying in the notebook
    that it reloads rather than recomputes.

    *present_root* is the staging directory: existence is checked against what
    the export actually holds, not against a list of what it is believed to
    hold. *also_expected* names top-level entries the caller has not written
    *yet* — the session is written between the zarr and ``transcripts.parquet``
    on purpose, so a check made there would otherwise call the parquet missing
    every time. Returns ``(items, unrepaired)``, where each entry of
    *unrepaired* is ``(node id, path)`` for a dangling path with no such
    substitution — reported rather than silently left, since the notebook will
    fail on it.
    """
    from palms.utils.clustering_code import reload_clustering_code

    export_path = str(export_path).rstrip("/")
    cache_path = f"{export_path}/sdata_cached.zarr"
    out, unrepaired = [], []

    def _present(path: str) -> bool:
        relative = Path(path).relative_to(export_path)
        if relative.parts and relative.parts[0] in also_expected:
            return True
        return (present_root / relative).exists()

    for item in items:
        item = dict(item)
        dangling = [
            path for path in recorded_paths_inside([item], export_path).get(item.get("id"), [])
            if not _present(path)
        ]
        node_id = item.get("id") or ""
        if dangling and node_id.startswith("clustering:"):
            key = node_id.split(":", 1)[1]
            item["code"] = reload_clustering_code(
                key, cache_path,
                reason=(f"Clustering '{key}' came with the source dataset, whose "
                        "analysis/ folder this crop does not carry."),
            )
            item["stale"] = False
        elif dangling:
            unrepaired.extend((node_id, path) for path in dangling)
        out.append(item)

    return out, unrepaired


def build_session_attrs(ctx, *, carried_elements, cluster_labels,
                        graph_items, he_shape_yx=None, arms_shape_yx=None,
                        roi_count: int = 0) -> dict:
    """The attrs an exported `viewer_session` should hold. Pure — no I/O."""
    he_state = getattr(ctx, "he_state", None) or {}
    arms_state = getattr(ctx, "arms_state", None) or {}
    state = getattr(ctx, "state", None) or {}

    attrs = {
        # Titles only. The layer name is what patch overlays link against by
        # name, so dropping it would break the link even though the pixels moved.
        "he_filename": he_state.get("he_filename"),
        "arms_he_filename": arms_state.get("he_filename"),
        # See the module docstring: the export is a slice, not the file.
        "he_path": None,
        "arms_he_path": None,
        "arms_geojson_path": None,
        "arms_csv_path": None,
        # The element carries the placement; these say "no further flip".
        "flip_v": False,
        "flip_h": False,
        "arms_flip_v": False,
        "arms_flip_h": False,
        "he_shape_yx": list(he_shape_yx) if he_shape_yx else None,
        "arms_he_shape_yx": list(arms_shape_yx) if arms_shape_yx else None,
        "cluster_labels": cluster_labels or None,
        "roi_count": int(roi_count),
        "has_rank_genes": state.get("rank_genes_df") is not None,
        "rank_genes_groupby": state.get("rank_genes_groupby"),
        "has_roi_deg": False,
        "has_arms_tile_deg": False,
        "marker_genes_json": state.get("marker_genes_json"),
        "segmentation_source": "xenium",
        "external_images_ui": _strip_affine(
            getattr(ctx, "external_images_state", None), carried_elements),
        "patch_overlays_ui": _strip_affine(
            getattr(ctx, "patch_overlays_state", None), carried_elements),
        # The export is written in the current format, so every one-time
        # migration is already satisfied. Leaving these unset makes the first
        # launch re-run migrations against a store that never needed them.
        "migrated_to_adata": True,
        "migrated_landmarks_to_sdata": True,
        "migrated_rank_genes_to_adata": True,
        "migrated_deg_to_sdata": True,
    }
    if graph_items:
        attrs["prov_graph"] = graph_items
    return attrs


def write_export_session(staging_dir: Path, attrs: dict, graph_items) -> None:
    """Write the session group and the provenance sidecar into a staging export.

    Both copies, because CLAUDE.md's rule is to persist the graph wherever an
    artifact is persisted, and `app._load_prov_graph_items` prefers the sidecar.
    They are written from the same list, so the "sidecar is never legitimately
    smaller" check cannot trip on an export.
    """
    from palms.utils.zarr_safe import atomic_json, safe_group_update

    staging_dir = Path(staging_dir)
    cache = staging_dir / "sdata_cached.zarr"

    json_safe, dropped = _json_safe(attrs)
    if dropped:
        log.warning("export session keys could not be serialized and were dropped: %s",
                    ", ".join(sorted(dropped)))

    with safe_group_update(cache, "viewer_session") as (session, _stage):
        for key, value in json_safe.items():
            session.attrs[key] = value

    if graph_items:
        from palms.tabs._helpers import PROV_GRAPH_SIDECAR

        sidecar_dir = staging_dir / "viewer_cache"
        sidecar_dir.mkdir(parents=True, exist_ok=True)
        atomic_json(sidecar_dir / PROV_GRAPH_SIDECAR, graph_items)


def _json_safe(obj):
    """Reuse session.py's serialisability filter so both agree on what survives."""
    from palms.utils.session import _json_safe as _sess_json_safe

    return _sess_json_safe(obj)
