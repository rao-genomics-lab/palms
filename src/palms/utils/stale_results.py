"""Which stored artifacts a stale provenance step left behind.

A node is ``stale`` when an ancestor was re-recorded with different code after
that node last ran (:meth:`prov_graph.ProvGraph.upsert`), so the artifact on
disk was produced by an input that no longer exists. Until now that was
display-only: the Notebook tab drew a badge and nothing acted on it.

This module is the **id bridge** between the two halves, and it is the only
place the correspondence is written down. Provenance ids are colon-namespaced
around the artifact they produce (``clustering:<key>``, ``rank_genes:<key>``,
``cnv:<backend>``); :mod:`store_inventory` keys name a row on disk
(``obs:table/clustering_<key>``, ``uns:table/rank_genes_<key>``,
``sidecar:adata_cnv_cache_<backend>.h5ad``). Nothing stored on a ``ProvNode``
says which is which — ``Step.outputs`` names in-memory namespace bindings, not
files — so the mapping has to be declared.

Pure, like :mod:`store_inventory` itself: no Qt, no filesystem, no zarr. It is
handed a graph and the sections the Dataset tab already scanned, and returns
keys for that tab's existing planner to vet. **It never decides that something
may be deleted** — every key it returns was already marked ``deletable`` by the
inventory, and ``plan_deletion`` / ``assert_node_deletable`` still have the last
word.

Three rules carry the safety:

* **An unkeyed slot is spared while any sibling is fresh.** Squidpy results are
  written to a single ``uns`` name — ``nhood_enrichment``, ``co_occurrence``,
  ``ligrec`` (:func:`adata_persistence.save_nhood_to_adata` and friends) — while
  the *nodes* are keyed per clustering (``nhood:<key>``). A stale ``nhood:A`` and
  a fresh ``nhood:B`` therefore name the same bytes, and clearing on the stale
  one alone would destroy a result that is perfectly current. :data:`SHARED_UNS`
  lists every such slot in one place so the three cannot drift apart.
* **Unmapped means untouched.** A node whose id this module does not recognise
  contributes nothing and is reported in ``unmatched``. That is not a gap to be
  closed by a fallback rule — a guessed mapping deletes the wrong file — and it
  is the correct answer for the terminals, whose figures live in
  ``<data_path>/plots/`` (outside every :func:`store_inventory.deletable_roots`)
  and whose CSV exports went to a path the user chose in a file dialog and which
  was never recorded anywhere.
* **The recipe is not the result.** Nothing here removes a node from the graph.
  Clearing a stale result leaves the step in place, so the notebook still
  replays and recreates it — which is exactly what
  ``store_inventory._plan_warnings`` already tells a user who deletes a
  clustering column by hand.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

from palms.utils import store_inventory
from palms.utils.adata_persistence import (
    ARMS_DEG_CACHE,
    CLUSTER_LABELS_PREFIX,
    CLUSTERING_PREFIX,
    ROI_DEG_CACHE,
)

#: Spelled out rather than imported from :mod:`gene_analysis`, which owns it but
#: pulls in scanpy and matplotlib at import time — the same trade
#: ``store_inventory`` makes for the provenance-sidecar prefix. A test asserts
#: the two agree, so they cannot drift.
RANK_GENES_PREFIX = "rank_genes"

#: ``uns`` slots written under one fixed name by steps that are keyed per
#: clustering. Family prefix → the slot every member of that family overwrites.
#: Cleared only when *every* node in the family is stale; see the module
#: docstring.
SHARED_UNS = {
    "nhood": "nhood_enrichment",
    "cooccur": "co_occurrence",
    "ligrec": "ligrec",
}

#: Written once per ranking run and pointing at the most recent one, so it
#: belongs to the whole ``rank_genes`` family rather than to any single key.
#: ``rank_genes_groups`` is scanpy's unkeyed default, still present in caches
#: written before the keying landed (``gene_analysis.resolve_rank_key``).
_SHARED_RANK_UNS = ("rank_genes_groupby", "rank_genes_groups")

#: Same shape for CNV: the per-backend artifacts are keyed, the run registry and
#: the profile matrix are not. ``cnv_run_info`` is the legacy flat form.
_SHARED_CNV_UNS = ("cnv_runs", "cnv_run_info", "cnv")
_SHARED_CNV_OBSM = ("X_cnv",)

#: Recorded steps that persist nothing of their own, with the reason. Kept
#: explicit rather than falling through to "unrecognised", because the two are
#: different answers: this one is known to have no artifact, and — for
#: ``cnv:copykat_propagated`` — it must also stay out of the ``cnv`` family
#: test, or a stale propagation step alone would clear the run registry that
#: describes a CNV run which is still current.
_NO_ARTIFACT = {
    "preamble": "defines data_path and binds adata; writes nothing",
    "environment": "records package versions; writes nothing",
    "spatial_neighbors": "builds the graph on adata_norm, stored only inside "
                         "adata_norm_cache.h5ad",
    "rois": "ROI polygons are a zarr element, out of scope for this action",
    "cnv:copykat_propagated": "propagates an existing run; its clusterings carry "
                              "the columns",
}

#: Namespaces whose artifacts this action deliberately does not remove, with the
#: reason shown to the user. Figures go to ``<data_path>/plots/``, which is
#: outside every deletable root by design; exports went to a path chosen in a
#: file dialog and recorded nowhere; the rest are zarr elements or session
#: state, out of scope here.
_OUT_OF_SCOPE = {
    "plot": "figures under plots/ — outside the viewer's deletable directories",
    "export": "written to a path you chose in a save dialog, not recorded",
    "viewer": "viewer state, not an analysis result",
    "preview": "a template preview; nothing was stored",
    "he": "H&E registration lives in the session, not in a result cache",
    "extimg": "an image element, out of scope for this action",
    "roi_expression": "held in memory only",
}


@dataclass(frozen=True)
class StaleSelection:
    """What clearing the stale results would and would not touch.

    ``keys`` is the answer; the other three fields exist so the tab can say why
    a stale step is missing from it. A silent omission and a deliberate sparing
    look identical otherwise, and the sparing is the one a user has to be able
    to check.
    """

    keys: tuple[str, ...] = ()
    #: (node id, the inventory keys it contributed), for the ones that matched.
    matched: tuple[tuple[str, tuple[str, ...]], ...] = ()
    #: (node id, reason) for a stale node with nothing this action can remove.
    unmatched: tuple[tuple[str, str], ...] = ()
    #: (slot name, the fresh node still using it) for each spared shared slot.
    spared: tuple[tuple[str, str], ...] = ()

    @property
    def is_empty(self) -> bool:
        return not self.keys


def stale_ids(graph) -> tuple[str, ...]:
    """Ids of every stale node, in insertion order."""
    if graph is None:
        return ()
    return tuple(n.id for n in sorted(graph.nodes(), key=lambda n: n.seq) if n.stale)


def _family(node_id: str) -> str:
    """The namespace an id belongs to — the part before the first colon."""
    return node_id.partition(":")[0]


def _suffix(node_id: str) -> str:
    """Everything after the first colon: the clustering key, backend or name."""
    return node_id.partition(":")[2]


def _family_all_stale(graph, family: str) -> tuple[bool, str]:
    """Whether every node in *family* is stale, and the first fresh one if not.

    ``_NO_ARTIFACT`` members are excluded: a step that stores nothing cannot be
    the reason a shared slot is still current, in either direction.
    """
    fresh = [
        n.id for n in graph.nodes()
        if _family(n.id) == family and n.id not in _NO_ARTIFACT and not n.stale
    ]
    return (not fresh), (fresh[0] if fresh else "")


def artifact_names(node_id: str, graph) -> tuple[tuple[str, str], ...]:
    """``(inventory kind, short name)`` pairs the step *node_id* wrote.

    Names, not keys: an ``obs`` column's key embeds the table it lives in
    (``obs:table/x`` vs ``obs:custom_table/x``), and a sidecar's differs between
    ``viewer_cache/`` and the legacy in-store location. Resolving names against
    the real inventory in :func:`select_stale` handles all of those without this
    function having to know which store it is looking at.

    Returns ``()`` for anything unmapped — see the module docstring on why there
    is no fallback rule.
    """
    if node_id in _NO_ARTIFACT:
        return ()
    family, suffix = _family(node_id), _suffix(node_id)
    out: list[tuple[str, str]] = []

    if family == "clustering":
        # The bare scanpy twin and cluster_labels_<key> come along through
        # store_inventory._cascade_candidates, which already pairs them.
        out.append((store_inventory.OBS, f"{CLUSTERING_PREFIX}{suffix}"))
    elif family == "annotation":
        out.append((store_inventory.OBS, f"{CLUSTER_LABELS_PREFIX}{suffix}"))
    elif family == "rank_genes":
        out.append((store_inventory.UNS, f"{RANK_GENES_PREFIX}_{suffix}"))
        if _family_all_stale(graph, "rank_genes")[0]:
            out += [(store_inventory.UNS, n) for n in _SHARED_RANK_UNS]
    elif family in SHARED_UNS:
        if _family_all_stale(graph, family)[0]:
            out.append((store_inventory.UNS, SHARED_UNS[family]))
    elif family == "cnv":
        out.append((store_inventory.OBS, f"cnv_score_{suffix}"))
        out.append((store_inventory.SIDECAR, f"adata_cnv_cache_{suffix}.h5ad"))
        out.append((store_inventory.SIDECAR, f"cnv_{suffix}_result.json"))
        if _family_all_stale(graph, "cnv")[0]:
            out += [(store_inventory.UNS, n) for n in _SHARED_CNV_UNS]
            out += [(store_inventory.OBSM, n) for n in _SHARED_CNV_OBSM]
    elif node_id == "normalize":
        out.append((store_inventory.SIDECAR, "adata_norm_cache.h5ad"))
    elif node_id == "roi_deg":
        out.append((store_inventory.SIDECAR, ROI_DEG_CACHE))
    elif node_id == "arms:tile_deg":
        out.append((store_inventory.SIDECAR, ARMS_DEG_CACHE))

    return tuple(out)


def _short_name(node) -> str:
    """The name :func:`artifact_names` matches against, per inventory kind.

    Table contents carry the bare column/key in ``name``. A sidecar's ``name``
    is decorated with the directory it was found in, which differs between the
    current and legacy locations, so its key is the stable half.
    """
    if node.kind == store_inventory.SIDECAR:
        return node.key.partition(":")[2].rpartition("/")[2]
    return node.name


#: What an id maps to that is simply not present. Distinct from "unmapped": the
#: first is a step whose result has already gone (or was never persisted on this
#: machine), the second is one this module has no rule for. Reporting both as
#: "nothing to remove" reads as a gap in the mapping when it is not one.
_NOT_IN_STORE = "its stored result is not in this dataset"


def _reason_unmatched(node_id: str, had_targets: bool = False) -> str:
    if node_id in _NO_ARTIFACT:
        return _NO_ARTIFACT[node_id]
    if had_targets:
        return _NOT_IN_STORE
    family = _family(node_id)
    if family in _OUT_OF_SCOPE:
        return _OUT_OF_SCOPE[family]
    if family == "arms":
        return "ARMS registration lives in the session, not in a result cache"
    return "no stored result this action knows how to remove"


def select_stale(graph, sections: Iterable) -> StaleSelection:
    """Inventory keys for the results of every stale step in *graph*.

    Only keys that are actually present in *sections* and already marked
    ``deletable`` are returned, so a blocked or absent row can never reach
    ``store_inventory.plan_deletion`` — which treats one as a bug rather than a
    user error, and refuses the whole batch.
    """
    ids = stale_ids(graph)
    if not ids:
        return StaleSelection()

    sections = list(sections)
    # (kind, short name) → keys, over the deletable rows only.
    index: dict[tuple[str, str], list[str]] = {}
    for section in sections:
        for node in section.nodes:
            if not node.deletable:
                continue
            index.setdefault((node.kind, _short_name(node)), []).append(node.key)

    keys: list[str] = []
    matched: list[tuple[str, tuple[str, ...]]] = []
    unmatched: list[tuple[str, str]] = []
    spared: dict[str, str] = {}

    for node_id in ids:
        found: list[str] = []
        targets = artifact_names(node_id, graph)
        for target in targets:
            found += index.get(target, [])
        # Deduplicate while keeping the order the mapping declared them in.
        found = list(dict.fromkeys(found))
        if found:
            matched.append((node_id, tuple(found)))
            keys += found
        else:
            unmatched.append((node_id, _reason_unmatched(node_id, bool(targets))))

        family = _family(node_id)
        if family in SHARED_UNS or family in ("rank_genes", "cnv"):
            all_stale, fresh = _family_all_stale(graph, family)
            if not all_stale:
                slot = SHARED_UNS.get(family, f"{family} shared state")
                spared[slot] = fresh

    return StaleSelection(
        keys=tuple(dict.fromkeys(keys)),
        matched=tuple(matched),
        unmatched=tuple(unmatched),
        spared=tuple(sorted(spared.items())),
    )


def looks_like_a_path_rewrite(graph) -> bool:
    """Whether the staleness looks like a moved dataset rather than a real edit.

    ``app.py`` re-emits the ``preamble`` node for the current ``data_path`` on
    every launch, so the first launch after a dataset is moved or renamed by
    hand revises it and flags **every** descendant stale — for nothing; the
    results are fine and ``palms-rename-dataset --repair`` is the fix. Left
    unsaid, "clear all stale results" would then be a one-click delete of every
    analysis in the dataset.

    A heuristic, and only ever used to warn: the signature is that nothing with
    a dependency survived, which no ordinary re-run produces once a session has
    more than one branch.
    """
    if graph is None:
        return False
    with_deps = [n for n in graph.nodes() if n.deps]
    if len(with_deps) < 2:
        return False
    return all(n.stale for n in with_deps)


def describe_selection(sel: StaleSelection, plots_note: Optional[str] = None) -> str:
    """The tab's report text: what was ticked, what was spared, what was left."""
    lines: list[str] = []
    if sel.matched:
        lines.append(f"Stale steps with results on disk ({len(sel.matched)}):")
        for node_id, keys in sel.matched:
            lines.append(f"  • {node_id}")
            for key in keys:
                lines.append(f"      {key}")
    if sel.spared:
        lines += ["", "Left alone — a step that is still current shares them:"]
        for slot, fresh in sel.spared:
            lines.append(f"  • uns['{slot}'] — still in use by {fresh}")
    if sel.unmatched:
        lines += ["", f"Stale steps with nothing to remove here ({len(sel.unmatched)}):"]
        for node_id, reason in sel.unmatched:
            lines.append(f"  • {node_id} — {reason}")
        if plots_note:
            lines += ["", plots_note]
    if not lines:
        lines.append("No stale steps in the provenance graph.")
    return "\n".join(lines)
