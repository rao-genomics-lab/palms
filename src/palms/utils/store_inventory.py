"""What a viewer dataset has on disk, and which of it may be deleted.

A dataset accumulates viewer-created data in four places the user never sees:
the SpatialData cache, ``viewer_cache/`` beside it, ``transcript_cache/`` from
``palms-preprocess``, and sibling backup stores. This module inventories all
four *plus* the raw 10x output, so the Dataset tab can show where the space went
and offer to reclaim the parts the viewer made.

Two rules keep that safe, and both live here rather than in the Qt code:

* :func:`assert_deletable` is the single choke point. A path may be deleted only
  if it resolves inside one of :func:`deletable_roots` — directories the viewer
  itself created. The raw Xenium output is *listed* but can never be selected,
  and that is a property the tests check over every node the inventory produces,
  not a promise the UI is trusted to keep.
* Anything unrecognised defaults to **not** deletable: an unknown entry in the
  dataset directory is raw output, an unknown element or obs column is left
  alone. The ``loader._USER_*`` allow-lists decide the yes cases only.

Like :func:`cache_repair.verify` this is filesystem-only — no SpatialData, no
AnnData — so it still reports on a store too broken to open. And it never
writes: no mkdir, no rename, no unlink. An inventory that can mutate the store
is a defect by definition, which is why a test greps this file for those calls.
"""

from __future__ import annotations

import json
import logging
import os
import stat
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Optional

from palms import loader
from palms.utils import adata_persistence, cache_repair, zarr_safe
from palms.utils.cache_repair import human_bytes

log = logging.getLogger(__name__)

# ── Node kinds ───────────────────────────────────────────────────────────────
GROUP = "GROUP"        # a heading with children (a type dir, the trash)
RAW = "RAW"            # original 10x output — never deletable
ELEMENT = "ELEMENT"    # a SpatialData element directory
OBS = "OBS"            # a column inside a table's obs
UNS = "UNS"
OBSM = "OBSM"
SESSION = "SESSION"    # a viewer_session subgroup, array or attr
SIDECAR = "SIDECAR"    # a file in viewer_cache/ (or the legacy in-store spot)
DERIVED = "DERIVED"    # a whole regenerable directory, e.g. transcript_cache/
BACKUP = "BACKUP"      # a sibling sdata_cached_*_*.zarr
TRASH = "TRASH"        # .xv_trash, or one copy inside it

# How a deletion can be undone. Stated per node because the answer differs.
RECOVER_TRASH = "trash"      # one copy kept in .xv_trash
RECOVER_REBUILD = "rebuild"  # restorable from a backup during a cache rebuild
RECOVER_NONE = "none"        # gone

#: Kinds allowed to name a whole deletable root (see :func:`assert_deletable`).
_ROOT_KINDS = frozenset({BACKUP, TRASH, DERIVED})

#: Kinds whose deletion is a table rewrite, so they carry no path at all.
_PATHLESS_KINDS = frozenset({OBS, UNS, OBSM})

TRANSCRIPT_CACHE_DIRNAME = "transcript_cache"

#: Viewer outputs written into the dataset directory itself. Outside every
#: deletable root, so they are shown read-only rather than special-cased into
#: the predicate — a file-level exemption would puncture the containment model
#: for the sake of a few small files. Listed here rather than left to the raw
#: section because calling the viewer's own log "original 10x output, never
#: modified by the viewer" is simply false.
DERIVED_IN_DATA_DIR = {
    "analysis.py": "exported analysis code",
    "analysis_notebook.ipynb": "exported analysis notebook",
    "plots": "figures the viewer exported, plus CopyKAT run sentinels",
    # Without this key a 220 MB DegaFile export falls through to the raw section
    # and is labelled "original 10x output, never modified by the viewer" —
    # false, and the reason the size report would stop adding up.
    "degafiles": "Celldega DegaFile export · re-runnable from Tools -> Publish",
    "palms.log": "the viewer's own log",
    # Datasets opened before the rename carry the old filename. Dropping the key
    # would not hide the file — it would fall through to the raw section and be
    # labelled "original 10x output, never modified by the viewer", which is
    # both false and un-deletable.
    "xenium_viewer.log": "the viewer's own log (written before the rename)",
}


class NotDeletable(RuntimeError):
    """A path is not inside a directory the viewer created."""


# ── The elements that must never be deletable ────────────────────────────────
# Deleting the table bricks the dataset; the rest break Crop Export and the
# Segmentation tab's revert-to-Xenium path. Blocking them also means the tab
# never has to tear down the load-bearing napari layers.
_CORE_REASONS = {
    "tables/table": "core — deleting this destroys the dataset",
    "labels/cell_labels": "core — cell rasters and Segmentation revert need it",
    "labels/nucleus_labels": "core — required by Crop Export",
    "images/morphology_focus": "core — required by Crop Export",
    "points/transcripts": "core — required by Crop Export",
}
#: Derived from the reason map so the two cannot drift apart.
CORE_ELEMENTS = tuple(_CORE_REASONS)

_UNRECOGNISED = "not recognised as viewer-created — left alone"

# ── Table contents ───────────────────────────────────────────────────────────
# uns keys a viewer action produced. `cnv_run_info` is the legacy flat form.
_DELETABLE_UNS = frozenset(loader._USER_UNS_KEYS) | {
    "rank_genes_groupby", "cnv_runs", "cnv_run_info", "cnv",
}
# Rankings are keyed per clustering (`rank_genes_<key>`), so membership in a
# fixed set is not enough — a name test is. Without it a keyed ranking showed up
# as "not recognised as viewer-created" and could not be deleted at all.
_DELETABLE_UNS_PREFIXES = ("rank_genes",)


def _is_deletable_uns(name: str) -> bool:
    return name in _DELETABLE_UNS or name.startswith(_DELETABLE_UNS_PREFIXES)
# X_umap is a copy of analysis/umap/ from the raw output; obsm['spatial'] is
# structural and every spatial step depends on it.
_DELETABLE_OBSM = frozenset({"X_umap", "X_cnv"})

_OBS_DETAIL = {
    "clustering_": "clustering",
    "cluster_labels_": "cluster labels",
    "cnv_score": "CNV score",
    "copykat_leiden_res": "CopyKAT clustering",
}

#: Columns the Xenium loader itself produces. Never paired as a clustering twin
#: (below), however a clustering happens to be named.
_STRUCTURAL_OBS = frozenset({
    "region", "instance_id", "cell_id", "cell_labels", "z_level",
    "transcript_counts", "control_probe_counts", "control_codeword_counts",
    "unassigned_codeword_counts", "deprecated_codeword_counts",
    "genomic_control_counts", "total_counts", "cell_area", "nucleus_area",
    "nucleus_count", "segmentation_method",
})


def _clustering_twin_of(column: str, columns: set[str]) -> str:
    """The ``clustering_<column>`` this bare column is the raw half of, or "".

    A Leiden run leaves *two* obs columns of the same data: the recorded step
    writes ``adata.obs[$key]`` so the notebook reproduces it, and
    ``save_clustering_to_adata`` writes ``clustering_<key>`` for the viewer.
    Offering only the prefixed one meant "delete this clustering" left an
    identical copy behind — so the bare twin is deletable, and cascades with it.
    """
    if column in _STRUCTURAL_OBS or column.startswith(loader._USER_OBS_PREFIXES):
        return ""
    twin = f"{adata_persistence.CLUSTERING_PREFIX}{column}"
    return twin if twin in columns else ""

# ── Sidecars ─────────────────────────────────────────────────────────────────
_SIDECAR_DETAIL = {
    "adata_norm_cache.h5ad": "normalised expression · recomputed on demand",
    "roi_deg_cache.parquet": "ROI differential expression · recomputed on demand",
    "arms_tile_deg_cache.parquet": "ARMS tile DEG · recomputed on demand",
    "cnv_copykat_input.h5ad": "CopyKAT worker input · not needed once a run finishes",
    "cnv_copykat_params.json": "CopyKAT worker parameters",
    "dega_staging": ("DegaFile export staging · symlinks to the raw output plus "
                     "what celldega unpacked · only makes the next export slower"),
}
#: Sidecars that look ordinary but must not be offered. The provenance graph is
#: the notebook's source of truth, and it *wins* over the copy in the session
#: attrs on load — so deleting it silently loses every step recorded since the
#: last session save. Spelled out rather than imported from ``tabs._helpers``,
#: which owns the constant but pulls in Qt.
#: Matched as a prefix, so the dated ``prov_graph.backup_*.json`` copies are
#: protected too — they exist precisely to survive a graph going wrong.
_BLOCKED_SIDECAR_PREFIX = "prov_graph"
_BLOCKED_SIDECAR_REASON = "the provenance graph — the notebook's source of truth"


def _blocked_sidecar_reason(name: str) -> str:
    if name.startswith(_BLOCKED_SIDECAR_PREFIX):
        return _BLOCKED_SIDECAR_REASON
    return ""

# ── The session-state trap ───────────────────────────────────────────────────
# ``save_session`` rebuilds viewer_session from ctx.state / ctx.he_state /
# ctx.arms_state at exit, and carries every unrecognised attr forward from the
# previous attrs. Deleting a session node on disk without clearing its in-memory
# mirror therefore puts it straight back on the next save. Each entry is
# (name of the dict on ViewerContext, key in that dict).
SESSION_MEMORY: dict[str, tuple[tuple[str, str], ...]] = {
    "session:group/he": (
        ("he_state", "affine_3x3"), ("he_state", "coarse_affine"),
        ("he_state", "he_filename"), ("he_state", "he_path"),
        ("he_state", "he_shape_yx"), ("he_state", "flip_v"), ("he_state", "flip_h"),
    ),
    "session:group/arms": (
        ("arms_state", "affine_3x3"), ("arms_state", "he_filename"),
        ("arms_state", "he_path"), ("arms_state", "he_shape_yx"),
        ("arms_state", "flip_v"), ("arms_state", "flip_h"),
        ("arms_state", "geojson_path"), ("arms_state", "csv_path"),
    ),
    "session:attr/cluster_labels": (("state", "cluster_labels"),),
    "session:attr/marker_genes_json": (("state", "marker_genes_json"),),
    "session:attr/external_images_ui": (("state", "external_images_ui"),),
    "session:attr/patch_overlays_ui": (("state", "patch_overlays_ui"),),
}

#: Session attrs offered as rows. Everything else in the group's attrs is a
#: migration marker or a derived flag, where deletion means nothing.
_SESSION_ATTRS = {
    "cluster_labels": "per-cluster names you typed",
    "marker_genes_json": "saved marker-gene sets",
    "external_images_ui": "contrast/opacity for registered images",
    "patch_overlays_ui": "display settings for patch overlays",
}
#: Attrs listed read-only. `prov_graph` is blocked for the same reason
#: `prov_graph.json` is.
_BLOCKED_SESSION_ATTRS = {
    "prov_graph": "the provenance graph — the notebook's source of truth",
    "segmentation_source": "which segmentation is active — not data",
}


def session_memory_keys(key: str) -> tuple[tuple[str, str], ...]:
    """In-memory ``(ctx dict, key)`` pairs to clear alongside session node *key*."""
    return SESSION_MEMORY.get(key, ())


# ── Model ────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class Node:
    """One row of the inventory.

    ``path`` is the filesystem object whose removal *is* this node's deletion,
    and it is deliberately ``None`` for table contents and session attrs:
    dropping an obs column is a rewrite of the whole table, so handing out a
    path would let an executor ``unlink`` it and corrupt ``tables/table``.
    """

    key: str
    kind: str
    name: str
    detail: str = ""
    size_bytes: Optional[int] = None
    path: Optional[Path] = None
    deletable: bool = False
    blocked_reason: str = ""      # "" iff deletable
    recoverable: str = RECOVER_NONE
    cascade: tuple[str, ...] = ()
    parent: str = ""              # key of the parent row, "" at top level


@dataclass(frozen=True)
class Section:
    title: str
    note: str = ""
    nodes: tuple[Node, ...] = ()

    @property
    def total_bytes(self) -> int:
        """Bytes in this section, counting only top-level rows to avoid
        double-counting children that live inside their parent."""
        return sum(n.size_bytes or 0 for n in self.nodes if not n.parent)


@dataclass(frozen=True)
class Plan:
    """A vetted batch of deletions, in the order the executor must apply them."""

    nodes: tuple[Node, ...] = ()
    added: tuple[Node, ...] = ()          # pulled in by cascade, not ticked
    total_bytes: int = 0
    unrecoverable: tuple[Node, ...] = ()
    warnings: tuple[str, ...] = ()
    dropped: tuple[str, ...] = ()         # selected keys that no longer exist

    @property
    def is_empty(self) -> bool:
        return not self.nodes

    def of_kinds(self, *kinds: str) -> tuple[Node, ...]:
        return tuple(n for n in self.nodes if n.kind in kinds)


TABLE_KINDS = (OBS, UNS, OBSM)

# Table edits run *first*. With custom segmentation active, deleting
# tables/custom_table and then re-persisting obs makes `_persist_table` route to
# `_persist_custom_table`, which recreates the element we just removed. Backups
# run last because they are the Cache tab's recovery source: if anything earlier
# fails, the safety net is still there.
_KIND_ORDER = {
    OBS: 0, UNS: 0, OBSM: 0,
    SESSION: 1,
    ELEMENT: 2,
    SIDECAR: 3, DERIVED: 3,
    TRASH: 4,
    BACKUP: 5,
}


# ── Sizes ────────────────────────────────────────────────────────────────────

def _entry_size(path: Path, *, skip: frozenset = frozenset()) -> Optional[int]:
    """Bytes *path* itself occupies, or ``None`` if it cannot be measured.

    Not shared with :func:`zarr_safe._dir_size`, which walks directories only
    (a 12 GB parquet would report 0 B), follows symlinks when sizing, and cannot
    prune the store's internal directories.

    A symlink counts as 0 and is never descended, so no size in this tree can be
    inflated by something outside the viewer's own directories — the same rule
    :func:`assert_deletable` enforces for deletion.
    """
    try:
        info = os.lstat(path)
    except OSError:
        return None
    if stat.S_ISLNK(info.st_mode):
        return 0
    if not stat.S_ISDIR(info.st_mode):
        return info.st_size
    total = 0
    stack = [str(path)]
    while stack:
        try:
            entries = list(os.scandir(stack.pop()))
        except OSError:
            continue
        for entry in entries:
            if entry.name in skip:
                continue
            try:
                if entry.is_symlink():
                    continue
                if entry.is_dir(follow_symlinks=False):
                    stack.append(entry.path)
                else:
                    total += entry.stat(follow_symlinks=False).st_size
            except OSError:
                continue
    return total


def _recoverable_from_trash(size: Optional[int]) -> str:
    """Whether ``safe_delete_element`` will really keep a copy of this.

    ``prune_trash`` drops a backup above its budget, so a big element is not
    recoverable however much the UI would like to say it is. Derived from the
    constant so the row stays honest if the budget changes.
    """
    if size is None or size > zarr_safe.DEFAULT_MAX_TRASH_BYTES:
        return RECOVER_NONE
    return RECOVER_TRASH


# ── Section 1: the raw Xenium output ─────────────────────────────────────────

def _viewer_owned_names(data_path: Path, cache_path: Optional[Path]) -> set[str]:
    """Top-level names in *data_path* that the viewer created."""
    names = {
        adata_persistence.SIDECAR_DIRNAME,
        TRANSCRIPT_CACHE_DIRNAME,
        *DERIVED_IN_DATA_DIR,
    }
    if cache_path is not None:
        names.add(Path(cache_path).name)
    for pattern in cache_repair.BACKUP_PATTERNS:
        names.update(p.name for p in data_path.glob(pattern))
    return names


def _raw_section(data_path: Path, cache_path: Optional[Path]) -> Section:
    """Everything in the dataset directory the viewer did not put there.

    By exclusion, never by a hardcoded Xenium file list: a vendor file nobody
    has heard of must show up read-only, not become selectable.
    """
    from palms.loader import is_cache_only

    owned = _viewer_owned_names(data_path, cache_path)
    # In a Crop Dataset export the "raw" files (experiment.xenium,
    # transcripts.parquet) were written by the viewer itself. They stay
    # read-only — they are still the only copy — but calling them 10x output
    # the viewer never touched is simply untrue, and it is the same wrong
    # mental model that let Force Rebuild loose on these datasets.
    cache_only = is_cache_only(data_path)
    reason = (
        "the only copy of this dataset's source data — a Crop Dataset export "
        "has no 10x output to fall back on"
        if cache_only else
        "original 10x output — the viewer never modifies it"
    )
    nodes: list[Node] = []
    try:
        entries = sorted(data_path.iterdir())
    except OSError as exc:
        return Section("Original Xenium output", f"Could not be listed: {exc}")
    for entry in entries:
        if entry.name in owned:
            continue
        nodes.append(Node(
            key=f"raw:{entry.name}",
            kind=RAW,
            name=entry.name,
            size_bytes=_entry_size(entry),
            path=entry,
            deletable=False,
            blocked_reason=reason,
        ))
    total = human_bytes(sum(n.size_bytes or 0 for n in nodes))
    title = "Dataset source files" if cache_only else "Original Xenium output"
    return Section(
        title,
        f"{total} · read-only, never modified by the viewer",
        tuple(nodes),
    )


# ── Section 2: the SpatialData cache ─────────────────────────────────────────

def _element_block_reason(element: str) -> str:
    """"" if this element may be deleted, else why not."""
    if element in _CORE_REASONS:
        return _CORE_REASONS[element]
    etype, _, name = element.partition("/")
    if name in (adata_persistence.CUSTOM_LABELS_KEY, adata_persistence.CUSTOM_TABLE_KEY):
        return ""
    keys: list = []
    if etype == "shapes":
        keys = loader._USER_SHAPE_KEYS
    elif etype == "images":
        keys = loader._USER_IMAGE_KEYS
    if loader._is_user_element(name, keys):
        return ""
    return _UNRECOGNISED


def _element_detail(element: str) -> str:
    _, _, name = element.partition("/")
    for labels in (loader._SHAPE_LABELS, loader._IMAGE_LABELS):
        if name in labels:
            return labels[name]
    if name.startswith("ext_"):
        return "registered external image"
    if name.startswith("patch_"):
        return "patch overlay"
    if name.endswith(loader._USER_SHAPE_SUFFIXES):
        return "registration landmarks"
    if name == adata_persistence.CUSTOM_LABELS_KEY:
        return "custom segmentation raster"
    if name == adata_persistence.CUSTOM_TABLE_KEY:
        return "custom segmentation cell table"
    return ""


def _cluster_counts(cache_path: Path, table: str) -> dict[str, int]:
    """``{clustering column: distinct cluster count}``, read straight from zarr.

    A count that cannot be read is a missing detail, never a failed inventory —
    ``read_obs_columns`` needs a readable ``_index`` and a zarr-v3 layout, and
    the whole point of this module is to work on a store that has neither.
    """
    try:
        columns = cache_repair.read_obs_columns(
            cache_path, ("clustering_",), table=table)
    except Exception as exc:
        log.debug("could not read cluster counts from %s: %s", cache_path, exc)
        return {}
    counts: dict[str, int] = {}
    for column, (_index, values) in columns.items():
        try:
            counts[column] = len({v for v in values.tolist() if v not in ("", None)})
        except Exception:
            continue
    return counts


def _obs_nodes(cache_path: Path, element: str, parent: str) -> list[Node]:
    table = element.partition("/")[2]
    obs_dir = cache_path / element / "obs"
    if not obs_dir.is_dir():
        return []
    counts = _cluster_counts(cache_path, element)
    entries = sorted(_public_entries(obs_dir))
    columns = {e.name for e in entries}
    nodes: list[Node] = []
    for entry in entries:
        column = entry.name
        deletable = column.startswith(loader._USER_OBS_PREFIXES)
        detail = ""
        for prefix, label in _OBS_DETAIL.items():
            if column.startswith(prefix):
                detail = label
                break
        if column in counts:
            detail = f"{counts[column]} clusters"
        twin = _clustering_twin_of(column, columns)
        if twin:
            deletable = True
            detail = f"the raw scanpy output behind {twin}"
        nodes.append(Node(
            key=f"obs:{table}/{column}",
            kind=OBS,
            name=column,
            detail=detail,
            size_bytes=_entry_size(entry),
            path=None,          # deletion is a table rewrite, not an unlink
            deletable=deletable,
            blocked_reason="" if deletable else _UNRECOGNISED,
            recoverable=RECOVER_REBUILD if deletable else RECOVER_NONE,
            parent=parent,
        ))
    return nodes


def _uns_obsm_nodes(cache_path: Path, element: str, parent: str) -> list[Node]:
    table = element.partition("/")[2]
    nodes: list[Node] = []
    for sub, kind, is_allowed in (("uns", UNS, _is_deletable_uns),
                                  ("obsm", OBSM, _DELETABLE_OBSM.__contains__)):
        directory = cache_path / element / sub
        if not directory.is_dir():
            continue
        for entry in sorted(_public_entries(directory)):
            deletable = is_allowed(entry.name)
            nodes.append(Node(
                key=f"{sub}:{table}/{entry.name}",
                kind=kind,
                name=entry.name,
                detail="UMAP coordinates" if entry.name == "X_umap" else "",
                size_bytes=_entry_size(entry),
                path=None,
                deletable=deletable,
                blocked_reason="" if deletable else _UNRECOGNISED,
                recoverable=RECOVER_REBUILD if deletable else RECOVER_NONE,
                parent=parent,
            ))
    return nodes


#: zarr's own metadata, which is not content and must never be offered — the
#: group's ``zarr.json`` *is* the group. (v2's markers all start with a dot and
#: are dropped by the prefix rule below.)
_ZARR_METADATA = frozenset({"zarr.json"})


def _public_entries(directory: Path) -> list[Path]:
    """Entries of *directory* that name real content.

    Both files and directories: a 0-d ``uns`` scalar's on-disk shape varies
    between anndata/zarr versions, and accepting only directories silently
    dropped keys.
    """
    try:
        return [p for p in directory.iterdir()
                if not p.name.startswith(("_", "."))
                and p.name not in _ZARR_METADATA]
    except OSError:
        return []


def _cache_section(cache_path: Optional[Path]) -> Section:
    if cache_path is None or not cache_path.is_dir():
        return Section(
            "Viewer cache",
            "No zarr cache for this dataset — nothing here to reclaim.",
        )
    elements = cache_repair._disk_elements(cache_path)
    nodes: list[Node] = []
    by_type: dict[str, list[Node]] = {}
    for element in elements:
        etype, _, name = element.partition("/")
        reason = _element_block_reason(element)
        size = _entry_size(cache_path / element)
        group_key = f"group:{etype}"
        node = Node(
            key=f"element:{element}",
            kind=ELEMENT,
            name=name,
            detail=_element_detail(element) if not reason else reason,
            size_bytes=size,
            path=cache_path / element,
            deletable=not reason,
            blocked_reason=reason,
            recoverable=_recoverable_from_trash(size) if not reason else RECOVER_NONE,
            parent=group_key,
        )
        by_type.setdefault(etype, []).append(node)

    for etype in zarr_safe.ELEMENT_TYPES:
        children = by_type.get(etype)
        if not children:
            continue
        nodes.append(Node(
            key=f"group:{etype}",
            kind=GROUP,
            name=etype,
            size_bytes=sum(c.size_bytes or 0 for c in children),
            deletable=False,
            blocked_reason="a group of elements — select the ones you want",
        ))
        nodes.extend(children)
        if etype == "tables":
            for child in children:
                element = f"tables/{child.name}"
                nodes.extend(_obs_nodes(cache_path, element, child.key))
                nodes.extend(_uns_obsm_nodes(cache_path, element, child.key))

    described = {}
    try:
        described = cache_repair.describe_store(cache_path)
    except OSError:
        pass
    # The section total is the sum of the visible rows, not describe_store's own
    # walk: a header that disagrees with its rows reads as a bug in a
    # disk-usage view (and it saves a second full walk of a 30 GB store).
    total = sum(n.size_bytes or 0 for n in nodes if n.kind == GROUP)
    note = f"{cache_path.name} · {human_bytes(total)}"
    if described:
        note += f" · {human_bytes(described['free_bytes'])} free on this disk"
    return Section("Viewer cache", note, tuple(nodes))


# ── Section 3: session state ─────────────────────────────────────────────────

def _session_attrs(cache_path: Path) -> dict:
    """viewer_session attrs, parsed with json so a broken store still reports."""
    marker = cache_path / "viewer_session" / "zarr.json"
    try:
        return json.loads(marker.read_text()).get("attributes", {}) or {}
    except (OSError, ValueError):
        return {}


def _session_section(cache_path: Optional[Path]) -> Section:
    if cache_path is None or not (cache_path / "viewer_session").is_dir():
        return Section("Session state", "No saved session for this dataset.")
    session_dir = cache_path / "viewer_session"
    nodes: list[Node] = []

    _GROUP_DETAIL = {
        "he": "H&E registration affine",
        "arms": "ARMS overlay affine",
        "rois": "legacy ROI arrays (ROIs now live in shapes/rois)",
        "clusterings": "legacy clustering parquets (now in the table)",
    }
    for entry in sorted(_public_entries(session_dir)):
        if entry.is_dir():
            key = f"session:group/{entry.name}"
            detail = _GROUP_DETAIL.get(entry.name, "")
        else:
            key = f"session:file/{entry.name}"
            detail = "legacy session file"
        nodes.append(Node(
            key=key,
            kind=SESSION,
            name=entry.name,
            detail=detail,
            size_bytes=_entry_size(entry),
            path=entry,
            deletable=True,
            recoverable=RECOVER_NONE,
        ))

    attrs = _session_attrs(cache_path)
    for name, detail in _SESSION_ATTRS.items():
        if attrs.get(name) in (None, {}, []):
            continue
        nodes.append(Node(
            key=f"session:attr/{name}",
            kind=SESSION,
            name=name,
            detail=detail,
            path=None,          # an attrs rewrite, not an unlink
            deletable=True,
            recoverable=RECOVER_NONE,
        ))
    for name, reason in _BLOCKED_SESSION_ATTRS.items():
        if attrs.get(name) in (None, {}, []):
            continue
        nodes.append(Node(
            key=f"session:attr/{name}",
            kind=SESSION,
            name=name,
            detail=reason,
            path=None,
            deletable=False,
            blocked_reason=reason,
        ))
    return Section(
        "Session state",
        "viewer_session/ · ROIs, registration and the names you typed",
        tuple(nodes),
    )


# ── Section 4: derived caches ────────────────────────────────────────────────

def _derived_section(data_path: Path, cache_path: Optional[Path]) -> Section:
    nodes: list[Node] = []

    # viewer_cache/ is a directory the viewer created outright, so everything in
    # it is viewer output — no allow-list needed, only the blocked list.
    # sidecar_dir defaults to not creating the directory, and it must stay that
    # way here: an inventory of a dataset must not bring it into existence.
    sidecar_home = adata_persistence.sidecar_dir(data_path)
    for entry in sorted(_public_entries(sidecar_home)):
        reason = _blocked_sidecar_reason(entry.name)
        nodes.append(Node(
            key=f"sidecar:{entry.name}",
            kind=SIDECAR,
            name=(f"{adata_persistence.SIDECAR_DIRNAME}/{entry.name}"
                  + ("/" if entry.is_dir() else "")),
            detail=reason or _sidecar_detail(entry.name),
            size_bytes=_entry_size(entry),
            path=entry,
            deletable=not reason,
            blocked_reason=reason,
            recoverable=RECOVER_NONE,
        ))

    # Legacy sidecars still inside the store, so pre-migration datasets are not
    # missing rows. They are inside the cache root, hence deletable.
    if cache_path is not None and cache_path.is_dir():
        for name in cache_repair._find_sidecars(cache_path):
            reason = _blocked_sidecar_reason(name)
            nodes.append(Node(
                key=f"sidecar:store/{name}",
                kind=SIDECAR,
                name=f"{cache_path.name}/{name}",
                detail=reason or _sidecar_detail(name) or "legacy in-store sidecar",
                size_bytes=_entry_size(cache_path / name),
                path=cache_path / name,
                deletable=not reason,
                blocked_reason=reason,
                recoverable=RECOVER_NONE,
            ))

    transcripts = data_path / TRANSCRIPT_CACHE_DIRNAME
    if transcripts.is_dir():
        count = len(list(transcripts.glob("*.feather")))
        nodes.append(Node(
            key=f"derived:{TRANSCRIPT_CACHE_DIRNAME}",
            kind=DERIVED,
            name=f"{TRANSCRIPT_CACHE_DIRNAME}/",
            detail=(f"{count} per-gene feather files · without them transcript "
                    "loading falls back to a ~5s/gene parquet scan until "
                    "palms-preprocess is re-run"),
            size_bytes=_entry_size(transcripts),
            path=transcripts,
            deletable=True,
            recoverable=RECOVER_NONE,
        ))

    for name, detail in DERIVED_IN_DATA_DIR.items():
        entry = data_path / name
        if not entry.exists():
            continue
        nodes.append(Node(
            key=f"derived:{name}",
            kind=DERIVED,
            name=f"{name}/" if entry.is_dir() else name,
            detail=detail,
            size_bytes=_entry_size(entry),
            path=entry,
            deletable=False,
            blocked_reason=("written into the dataset folder, outside the "
                            "viewer's deletable directories — remove it by hand"),
        ))

    total = human_bytes(sum(n.size_bytes or 0 for n in nodes))
    return Section(
        "Derived caches",
        f"{total} · regenerable analysis output",
        tuple(nodes),
    )


def _sidecar_detail(name: str) -> str:
    if name in _SIDECAR_DETAIL:
        return _SIDECAR_DETAIL[name]
    if name.startswith("adata_cnv_cache_"):
        backend = name[len("adata_cnv_cache_"):].removesuffix(".h5ad")
        return f"CNV profile matrix ({backend}) · needed for the heatmap"
    if name.startswith("cnv_") and name.endswith("_result.json"):
        backend = name[len("cnv_"):-len("_result.json")]
        return f"CNV run metadata ({backend})"
    return "viewer output"


# ── Section 5: backups and trash ─────────────────────────────────────────────

def _backup_section(data_path: Path, cache_path: Optional[Path]) -> Section:
    nodes: list[Node] = []
    if cache_path is not None:
        backups = cache_repair.find_backups(cache_path)
    else:
        backups = sorted(
            (p for pattern in cache_repair.BACKUP_PATTERNS
             for p in data_path.glob(pattern) if p.is_dir()),
            reverse=True,
        )
    for backup in backups:
        nodes.append(Node(
            key=f"backup:{backup.name}",
            kind=BACKUP,
            name=backup.name,
            detail="a whole previous cache · the Cache tab recovers from these",
            size_bytes=_entry_size(backup),
            path=backup,
            deletable=True,
            recoverable=RECOVER_NONE,
        ))

    if cache_path is not None and (cache_path / zarr_safe.TRASH_DIR).is_dir():
        trash_dir = cache_path / zarr_safe.TRASH_DIR
        copies = zarr_safe.list_trash(cache_path)
        count = sum(len(v) for v in copies.values())
        nodes.append(Node(
            key="trash:all",
            kind=TRASH,
            name=zarr_safe.TRASH_DIR,
            detail=(f"{count} previous version(s) of overwritten elements · "
                    "the Cache tab recovers from these"),
            size_bytes=_entry_size(trash_dir),
            path=trash_dir,
            deletable=True,
            recoverable=RECOVER_NONE,
        ))
        for element, paths in sorted(copies.items()):
            for path in paths:
                nodes.append(Node(
                    key=f"trash:{element}/{path.name}",
                    kind=TRASH,
                    name=f"{element} — {path.name}",
                    detail="one previous version",
                    size_bytes=_entry_size(path),
                    path=path,
                    deletable=True,
                    recoverable=RECOVER_NONE,
                    parent="trash:all",
                ))

    total = human_bytes(sum(n.size_bytes or 0 for n in nodes if not n.parent))
    return Section(
        "Backups & trash",
        f"{total} · ⚠ this is what the Cache tab recovers from",
        tuple(nodes),
    )


# ── Cascades ─────────────────────────────────────────────────────────────────

def _cascade_candidates(node: Node) -> tuple[str, ...]:
    """Keys that must go with *node*, before checking they exist."""
    if node.kind == ELEMENT:
        if node.name.startswith("ext_"):
            # Matches tab_external_images.on_remove, which deletes both.
            return (f"element:shapes/{node.name}_xenium_lm",
                    f"element:shapes/{node.name}_image_lm")
        # The Segmentation tab needs both halves; one alone leaves a labels
        # layer with no table and a segmentation_source pointing at nothing.
        if node.name == adata_persistence.CUSTOM_LABELS_KEY:
            return (f"element:tables/{adata_persistence.CUSTOM_TABLE_KEY}",)
        if node.name == adata_persistence.CUSTOM_TABLE_KEY:
            return (f"element:labels/{adata_persistence.CUSTOM_LABELS_KEY}",)
        return ()
    if node.kind == OBS and node.name.startswith(adata_persistence.CLUSTERING_PREFIX):
        table = node.key.partition(":")[2].partition("/")[0]
        suffix = node.name[len(adata_persistence.CLUSTERING_PREFIX):]
        # The typed-in names, and the bare twin a Leiden run leaves behind.
        # _link_cascades drops whichever of these does not exist.
        return (f"obs:{table}/{adata_persistence.CLUSTER_LABELS_PREFIX}{suffix}",
                f"obs:{table}/{suffix}")
    return ()


def _link_cascades(sections: list[Section]) -> list[Section]:
    """Resolve cascade candidates against the real key set.

    Done as a second pass so a cascade can never name a node that is not there.
    """
    keys = {n.key for section in sections for n in section.nodes}
    return [
        replace(section, nodes=tuple(
            replace(node, cascade=tuple(
                k for k in _cascade_candidates(node) if k in keys and k != node.key
            ))
            for node in section.nodes
        ))
        for section in sections
    ]


def build_inventory(data_path, cache_path=None) -> list[Section]:
    """Every viewer-visible thing on disk for this dataset, in five sections.

    All five are always returned, empty ones carrying an explanatory note, so
    the tree code and the tests never need a special case for ``--no-cache``.
    """
    data_path = Path(data_path)
    cache_path = Path(cache_path) if cache_path is not None else None
    return _link_cascades([
        _raw_section(data_path, cache_path),
        _cache_section(cache_path),
        _session_section(cache_path),
        _derived_section(data_path, cache_path),
        _backup_section(data_path, cache_path),
    ])


# ── The safety predicate ─────────────────────────────────────────────────────

def deletable_roots(data_path, cache_path=None) -> tuple[Path, ...]:
    """Directories the viewer created for this dataset, fully resolved.

    Nothing outside these may ever be unlinked. This tuple *is* constraint 1:
    the raw 10x output is not inside any of them, so no amount of UI or
    executor bugs can reach it.
    """
    data_path = Path(data_path)
    try:
        anchor = data_path.resolve(strict=True)
    except OSError:
        return ()

    candidates: list[Path] = [
        data_path / adata_persistence.SIDECAR_DIRNAME,
        data_path / TRANSCRIPT_CACHE_DIRNAME,
    ]
    if cache_path is not None:
        cache_path = Path(cache_path)
        candidates.append(cache_path)
        candidates += cache_repair.find_backups(cache_path)
    else:
        candidates += sorted(
            p for pattern in cache_repair.BACKUP_PATTERNS
            for p in data_path.glob(pattern) if p.is_dir()
        )

    roots: list[Path] = []
    for candidate in candidates:
        try:
            # strict=True drops what does not exist, and resolve collapses ".."
            # and follows symlinks, both of which matter before containment.
            resolved = candidate.resolve(strict=True)
        except OSError:
            continue
        # A root that *is* the dataset directory, or contains it, would make the
        # whole 10x output deletable — exactly the failure this exists to
        # prevent. Reachable in practice: a cache path of "." or a
        # sdata_cached.zarr symlinked to its own parent.
        if resolved == anchor or anchor.is_relative_to(resolved):
            log.warning("refusing %s as a deletable root: it contains the dataset",
                        resolved)
            continue
        if resolved not in roots:
            roots.append(resolved)
    return tuple(roots)


def assert_deletable(path, roots: Iterable[Path], *, kind: Optional[str] = None) -> None:
    """Raise :class:`NotDeletable` unless *path* is inside a viewer-created root.

    The single choke point: nothing may touch the filesystem without passing
    through here first. Policy about whole roots lives *inside*, keyed on
    ``kind`` — a choke point that asks its callers for permission is not one.
    """
    roots = tuple(roots)
    if not roots:
        raise NotDeletable(
            f"{path}: this dataset has no viewer-created directories to delete from")
    try:
        resolved = Path(path).resolve()
    except OSError as exc:
        raise NotDeletable(f"{path}: could not be resolved ({exc})") from None

    if not any(resolved == root or resolved.is_relative_to(root) for root in roots):
        raise NotDeletable(
            f"{path} is not inside a directory the viewer created "
            f"({', '.join(str(r) for r in roots)}) — refusing to delete it")
    if resolved in roots and kind not in _ROOT_KINDS:
        raise NotDeletable(
            f"{path} is a whole viewer directory; only a "
            f"{'/'.join(sorted(_ROOT_KINDS))} node may name one (got kind={kind!r})")


def assert_node_deletable(node: Node, roots: Iterable[Path]) -> None:
    """The node-level choke point the executor calls before anything else."""
    if not node.deletable:
        raise NotDeletable(f"{node.name}: {node.blocked_reason or 'not deletable'}")
    if node.path is None:
        if node.kind not in _PATHLESS_KINDS and not node.key.startswith("session:attr/"):
            raise NotDeletable(
                f"{node.key}: kind {node.kind} needs a path and has none")
        return
    assert_deletable(node.path, roots, kind=node.kind)


# ── Planning a deletion ──────────────────────────────────────────────────────

def _subsumed(node: Node, others: Iterable[Node]) -> bool:
    """True when another chosen node's path already contains this one.

    Selecting ``.xv_trash`` *and* a copy inside it would otherwise rmtree the
    parent and then fail on a path that no longer exists.
    """
    if node.path is None:
        return False
    for other in others:
        if other.path is None or other.key == node.key or other.path == node.path:
            continue
        if node.path.is_relative_to(other.path):
            return True
    return False


def plan_deletion(sections: Iterable[Section], selected_keys: Iterable[str]) -> Plan:
    """Expand, vet and order a user selection into an applicable batch.

    Takes the sections rather than nodes because cascades are *keys*: resolving
    them needs the index.
    """
    index = {n.key: n for section in sections for n in section.nodes}
    warnings: list[str] = []
    dropped: list[str] = []

    queue = list(dict.fromkeys(selected_keys))
    seed = set(queue)
    chosen: dict[str, Node] = {}
    added: list[str] = []
    while queue:
        key = queue.pop(0)
        if key in chosen:
            continue
        node = index.get(key)
        if node is None:
            dropped.append(key)
            warnings.append(f"{key} is no longer present — skipped")
            continue
        if not node.deletable:
            # The UI cannot tick a blocked row, so reaching here is a
            # programming error: fail before any dialog is shown.
            raise NotDeletable(f"{node.name}: {node.blocked_reason or 'not deletable'}")
        chosen[key] = node
        for target in node.cascade:
            child = index.get(target)
            if child is not None and child.deletable and target not in chosen:
                queue.append(target)
                if target not in seed:
                    added.append(target)

    picked = list(chosen.values())
    kept = [n for n in picked if not _subsumed(n, picked)]
    kept.sort(key=lambda n: (_KIND_ORDER.get(n.kind, 9), n.key))

    warnings += _plan_warnings(kept)
    unrecoverable = tuple(sorted(
        (n for n in kept if n.recoverable == RECOVER_NONE),
        key=lambda n: -(n.size_bytes or 0),
    ))
    return Plan(
        nodes=tuple(kept),
        added=tuple(index[k] for k in added if index[k] in kept),
        total_bytes=sum(n.size_bytes or 0 for n in kept),
        unrecoverable=unrecoverable,
        warnings=tuple(dict.fromkeys(warnings)),
        dropped=tuple(dropped),
    )


def _plan_warnings(nodes: Iterable[Node]) -> list[str]:
    out: list[str] = []
    nodes = list(nodes)
    if any(n.kind == OBS and n.name.startswith(adata_persistence.CLUSTERING_PREFIX)
           for n in nodes):
        out.append("the provenance graph keeps its clustering step, so the "
                   "notebook still replays to recreate the column — that is correct")
    if any(n.kind in (BACKUP, TRASH) for n in nodes):
        out.append("the Cache tab recovers from backups and trash; removing them "
                   "leaves it nothing to recover from")
    if any(n.kind == ELEMENT and n.name == adata_persistence.CUSTOM_TABLE_KEY
           for n in nodes):
        out.append("custom segmentation goes with it — the viewer reverts to the "
                   "Xenium segmentation")
    cnv_uns = [n for n in nodes if n.kind == UNS and n.name.startswith("cnv")]
    if cnv_uns and not any(n.kind == SIDECAR and "cnv" in n.name for n in nodes):
        out.append("a CNV backend is re-materialised from its cnv_*_result.json "
                   "sidecar, so removing only the uns entry will not hide it")
    return out


def bundles(sections: Iterable[Section]) -> dict[str, tuple[str, ...]]:
    """Named groups of keys the dialog may offer as one checkbox.

    Not cascades: these are things a user may reasonably want whole or in part,
    so the choice stays theirs.
    """
    index = {n.key: n for section in sections for n in section.nodes}
    found: dict[str, tuple[str, ...]] = {}

    backends: set[str] = set()
    for key, node in index.items():
        if node.kind == SIDECAR and node.name.endswith("_result.json"):
            backends.add(Path(node.name).name[len("cnv_"):-len("_result.json")])
        if node.kind == OBS and node.name.startswith("cnv_score_"):
            backends.add(node.name[len("cnv_score_"):])
    for backend in sorted(backends):
        members = tuple(
            key for key, node in index.items()
            if node.deletable and (
                (node.kind == OBS and node.name == f"cnv_score_{backend}")
                or (node.kind == UNS and node.name in ("cnv_runs", "cnv_run_info"))
                or (node.kind == SIDECAR
                    and node.name.endswith((f"adata_cnv_cache_{backend}.h5ad",
                                            f"cnv_{backend}_result.json")))
            )
        )
        if members:
            found[f"CNV backend: {backend}"] = tuple(sorted(members))

    ranked = tuple(
        key for key, node in index.items()
        if node.deletable and (
            (node.kind == UNS and node.name.startswith("rank_genes"))
            or (node.kind == SIDECAR and node.name.endswith("adata_norm_cache.h5ad"))
        )
    )
    if ranked:
        found["Ranked genes + the normalised cache"] = tuple(sorted(ranked))
    return found


# ── The confirmation text ────────────────────────────────────────────────────

_KIND_TITLES = {
    OBS: "Table columns (obs):",
    UNS: "Table results (uns):",
    OBSM: "Table matrices (obsm):",
    SESSION: "Session state:",
    ELEMENT: "Store elements:",
    SIDECAR: "Sidecar files:",
    DERIVED: "Derived caches:",
    TRASH: "Trash:",
    BACKUP: "Backup stores:",
}


def describe_plan(plan: Plan) -> str:
    """The confirm-dialog body.

    Every path is spelled out: several of these are irreversible, and the
    dataset directory is the one thing a user can cross-check against ``du``.
    """
    if plan.is_empty:
        return "Nothing selected."
    out = [
        f"Delete {len(plan.nodes)} item(s), reclaiming about "
        f"{human_bytes(plan.total_bytes)}?",
        "",
    ]
    for kind, title in _KIND_TITLES.items():       # dict order == executor order
        group = plan.of_kinds(kind)
        if not group:
            continue
        out.append(title)
        for node in group:
            size = f"  ({human_bytes(node.size_bytes)})" if node.size_bytes else ""
            out.append(f"    {node.name}{size}")
            if node.path is not None:
                out.append(f"        {node.path}")
        out.append("")
    if plan.added:
        out.append("Also removed, because they only make sense together:")
        out += [f"    {n.name}" for n in plan.added]
        out.append("")
    if plan.unrecoverable:
        out.append("⚠ Not recoverable — no backup copy is kept of these:")
        out += [f"    {n.name}" for n in plan.unrecoverable]
        out.append("")
    out += [f"⚠ {w}" for w in plan.warnings]
    return "\n".join(out).rstrip()
