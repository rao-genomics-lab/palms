"""The bridge from a stale provenance node to the rows it wrote on disk.

Pure — no Qt, no filesystem, no zarr — like the module it tests. The graph is
built with the real :class:`ProvGraph` and the sections with the real
:class:`store_inventory.Node`, so a change to either shape breaks these tests
rather than passing against a mock that no longer resembles the thing.

The property test at the end is the one that matters: it asserts over *every*
key the selector can produce that the key exists and is deletable, rather than
over a remembered list of cases — the idiom ``test_store_inventory.py`` uses for
the same reason.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from palms.utils import stale_results, store_inventory  # noqa: E402
from palms.utils.prov_graph import SETUP, TERMINAL, ProvGraph  # noqa: E402


# ── Fixtures ─────────────────────────────────────────────────────────────────

def _graph():
    """preamble → normalize → two clusterings, each with a ranking and a nhood."""
    g = ProvGraph()
    g.upsert("preamble", "data_path = ...", kind=SETUP)
    g.upsert("normalize", "sc.pp.normalize_total(adata)", deps=["preamble"], kind=SETUP)
    g.upsert("spatial_neighbors", "sq.gr.spatial_neighbors_knn(...)",
             deps=["normalize"])
    for key in ("leiden_r1.0", "leiden_r0.5"):
        g.upsert(f"clustering:{key}", f"adata.obs['{key}'] = ...", deps=["normalize"])
        g.upsert(f"rank_genes:{key}", "sc.tl.rank_genes_groups(...)",
                 deps=["normalize", f"clustering:{key}"])
        g.upsert(f"nhood:{key}", "sq.gr.nhood_enrichment(...)",
                 deps=["spatial_neighbors", f"clustering:{key}"])
    return g


def _obs(name, *, deletable=True):
    return store_inventory.Node(
        key=f"obs:table/{name}", kind=store_inventory.OBS, name=name,
        deletable=deletable,
        blocked_reason="" if deletable else "not recognised as viewer-created",
    )


def _uns(name, *, deletable=True):
    return store_inventory.Node(
        key=f"uns:table/{name}", kind=store_inventory.UNS, name=name,
        deletable=deletable,
        blocked_reason="" if deletable else "not recognised as viewer-created",
    )


def _obsm(name):
    return store_inventory.Node(
        key=f"obsm:table/{name}", kind=store_inventory.OBSM, name=name, deletable=True)


def _sidecar(filename, *, store=False, deletable=True):
    key = f"sidecar:store/{filename}" if store else f"sidecar:{filename}"
    return store_inventory.Node(
        key=key, kind=store_inventory.SIDECAR,
        name=f"viewer_cache/{filename}", deletable=deletable,
        blocked_reason="" if deletable else "blocked",
    )


def _sections(*nodes):
    return [store_inventory.Section("Viewer cache", nodes=tuple(nodes))]


def _full_sections():
    """Everything the graph above could possibly have written."""
    return _sections(
        _obs("clustering_leiden_r1.0"), _obs("leiden_r1.0"),
        _obs("cluster_labels_leiden_r1.0"),
        _obs("clustering_leiden_r0.5"), _obs("leiden_r0.5"),
        _obs("cnv_score_infercnv"),
        _obs("total_counts", deletable=False),
        _uns("rank_genes_leiden_r1.0"), _uns("rank_genes_leiden_r0.5"),
        _uns("rank_genes_groupby"), _uns("nhood_enrichment"),
        _uns("co_occurrence"), _uns("ligrec"), _uns("cnv_runs"),
        _obsm("X_cnv"),
        _sidecar("adata_norm_cache.h5ad"),
        _sidecar("roi_deg_cache.parquet"),
        _sidecar("arms_tile_deg_cache.parquet"),
        _sidecar("adata_cnv_cache_infercnv.h5ad"),
        _sidecar("cnv_infercnv_result.json"),
    )


# ── The constant that is duplicated rather than imported ─────────────────────

def test_the_rank_genes_prefix_matches_the_module_that_owns_it():
    """stale_results spells it out to avoid importing scanpy; assert they agree."""
    gene_analysis = pytest.importorskip("palms.utils.gene_analysis")
    assert stale_results.RANK_GENES_PREFIX == gene_analysis.RANK_GENES_PREFIX


# ── stale_ids ────────────────────────────────────────────────────────────────

def test_no_stale_nodes_selects_nothing():
    g = _graph()
    assert stale_results.stale_ids(g) == ()
    assert stale_results.select_stale(g, _full_sections()).is_empty


def test_a_missing_graph_is_not_an_error():
    assert stale_results.stale_ids(None) == ()
    assert stale_results.select_stale(None, _full_sections()).is_empty


def test_re_running_an_upstream_step_stales_its_descendants():
    g = _graph()
    g.upsert("normalize", "sc.pp.normalize_total(adata, target_sum=100)",
             deps=["preamble"], kind=SETUP)
    ids = set(stale_results.stale_ids(g))
    assert "clustering:leiden_r1.0" in ids
    assert "rank_genes:leiden_r1.0" in ids
    assert "normalize" not in ids       # the re-recorded node is fresh again


# ── The mapping, one case per row of the table ───────────────────────────────

def test_a_stale_clustering_selects_its_column():
    g = _graph()
    g.upsert("clustering:leiden_r1.0", "changed", deps=["normalize"])
    g.upsert("normalize", "changed too", deps=["preamble"], kind=SETUP)
    sel = stale_results.select_stale(g, _full_sections())
    assert "obs:table/clustering_leiden_r1.0" in sel.keys
    # The bare twin and the typed names are pulled in by the inventory's own
    # cascade, not by this module — it must not name them itself.
    assert "obs:table/leiden_r1.0" not in sel.keys
    assert "obs:table/cluster_labels_leiden_r1.0" not in sel.keys


def test_a_stale_ranking_selects_its_keyed_slot():
    g = _graph()
    g.upsert("clustering:leiden_r1.0", "changed", deps=["normalize"])
    sel = stale_results.select_stale(g, _full_sections())
    assert "uns:table/rank_genes_leiden_r1.0" in sel.keys
    # rank_genes:leiden_r0.5 is still fresh, so the shared pointer stays.
    assert "uns:table/rank_genes_groupby" not in sel.keys


def test_a_stale_cnv_run_selects_its_backend_artifacts():
    g = _graph()
    g.upsert("cnv:infercnv", "cnv_pipeline(...)", deps=["normalize"])
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sel = stale_results.select_stale(g, _full_sections())
    assert "obs:table/cnv_score_infercnv" in sel.keys
    assert "sidecar:adata_cnv_cache_infercnv.h5ad" in sel.keys
    assert "sidecar:cnv_infercnv_result.json" in sel.keys
    # It is the only cnv node and it is stale, so the shared registry goes too.
    assert "uns:table/cnv_runs" in sel.keys
    assert "obsm:table/X_cnv" in sel.keys


def test_the_deg_sidecars_are_mapped():
    g = _graph()
    g.upsert("roi_deg", "welch(...)", deps=["normalize"])
    g.upsert("arms:tile_deg", "welch(...)", deps=["normalize"])
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    keys = stale_results.select_stale(g, _full_sections()).keys
    assert "sidecar:roi_deg_cache.parquet" in keys
    assert "sidecar:arms_tile_deg_cache.parquet" in keys


def test_the_normalisation_cache_goes_only_when_normalize_itself_is_stale():
    """upsert clears `stale` on the node it re-records, so changing normalize
    leaves normalize fresh; only a change further up can stale it."""
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    assert "sidecar:adata_norm_cache.h5ad" not in \
        stale_results.select_stale(g, _full_sections()).keys

    g = _graph()
    g.upsert("preamble", 'data_path = Path(r"/elsewhere")', kind=SETUP)
    assert "sidecar:adata_norm_cache.h5ad" in \
        stale_results.select_stale(g, _full_sections()).keys


def test_a_legacy_in_store_sidecar_is_matched_by_its_basename():
    """Its `name` is decorated with the directory, so the key is the stable half."""
    g = _graph()
    g.upsert("roi_deg", "welch(...)", deps=["normalize"])
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sections = _sections(_sidecar("roi_deg_cache.parquet", store=True))
    assert "sidecar:store/roi_deg_cache.parquet" in \
        stale_results.select_stale(g, sections).keys


# ── The unkeyed-slot rule ────────────────────────────────────────────────────

def test_a_shared_uns_slot_is_spared_while_a_sibling_is_fresh():
    """The defect this rule exists for: one nhood key stale, the other current."""
    g = _graph()
    g.upsert("clustering:leiden_r1.0", "changed", deps=["normalize"])
    sel = stale_results.select_stale(g, _full_sections())
    assert "nhood:leiden_r1.0" in dict(sel.unmatched)
    assert "uns:table/nhood_enrichment" not in sel.keys
    assert ("nhood_enrichment", "nhood:leiden_r0.5") in sel.spared


def test_a_shared_uns_slot_goes_when_every_member_is_stale():
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sel = stale_results.select_stale(g, _full_sections())
    assert "uns:table/nhood_enrichment" in sel.keys
    assert "uns:table/rank_genes_groupby" in sel.keys
    assert sel.spared == ()


def test_a_stale_propagation_step_alone_does_not_clear_the_cnv_registry():
    """cnv:copykat_propagated stores nothing, so it must not vote on the family."""
    g = _graph()
    g.upsert("cnv:copykat", "copykat(...)", deps=["normalize"])
    g.upsert("cnv:copykat_propagated", "propagate(...)", deps=["cnv:copykat"])
    g.upsert("cnv:copykat", "copykat(..., changed)", deps=["normalize"])
    sel = stale_results.select_stale(g, _full_sections())
    assert "cnv:copykat_propagated" in dict(sel.unmatched)
    assert "uns:table/cnv_runs" not in sel.keys


# ── Unmapped means untouched ─────────────────────────────────────────────────

@pytest.mark.parametrize("node_id,kind", [
    ("plot:nhood:leiden_r1.0", TERMINAL),
    ("export:rank_genes:leiden_r1.0", TERMINAL),
    ("he:landmark_register", TERMINAL),
    ("extimg:load:ext_dapi", TERMINAL),
    ("roi_expression:EPCAM", TERMINAL),
])
def test_out_of_scope_namespaces_contribute_no_keys(node_id, kind):
    g = _graph()
    g.upsert(node_id, "code", deps=["normalize"], kind=kind)
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sel = stale_results.select_stale(g, _full_sections())
    reasons = dict(sel.unmatched)
    assert node_id in reasons and reasons[node_id]


def test_an_unrecognised_id_is_reported_rather_than_guessed():
    g = _graph()
    g.upsert("something:new", "code", deps=["normalize"])
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sel = stale_results.select_stale(g, _full_sections())
    assert dict(sel.unmatched)["something:new"] == \
        "no stored result this action knows how to remove"


def test_a_row_that_is_not_deletable_is_never_selected():
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sections = _sections(_obs("clustering_leiden_r1.0", deletable=False))
    assert stale_results.select_stale(g, sections).keys == ()


def test_a_mapped_step_whose_result_is_gone_says_so_rather_than_unmapped():
    """"Already gone" and "no rule for this id" are different answers.

    Reporting both as "nothing to remove" reads as a hole in the mapping table,
    which sends the reader looking for a bug that is not there.
    """
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sections = _sections(_obs("clustering_leiden_r1.0"))   # no uns rows at all
    reasons = dict(stale_results.select_stale(g, sections).unmatched)
    assert reasons["rank_genes:leiden_r1.0"] == stale_results._NOT_IN_STORE
    assert reasons["nhood:leiden_r1.0"] == stale_results._NOT_IN_STORE
    # ...while a genuinely unmapped id still says so.
    g.upsert("something:new", "code", deps=["normalize"])
    g.upsert("preamble", "moved", kind=SETUP)
    reasons = dict(stale_results.select_stale(g, sections).unmatched)
    assert "knows how to remove" in reasons["something:new"]


def test_an_artifact_that_is_not_on_disk_is_simply_absent():
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sel = stale_results.select_stale(g, _sections())
    assert sel.keys == ()
    assert len(sel.unmatched) == len(stale_results.stale_ids(g))


# ── The property that matters ────────────────────────────────────────────────

def test_every_selected_key_exists_and_is_deletable():
    """Over every node the selector can produce, not a remembered list.

    plan_deletion raises NotDeletable on a blocked key and refuses the whole
    batch, which tab_dataset reports as a bug rather than a user error — so a
    selector that can emit one is a defect, whatever the id happens to be.
    """
    g = _graph()
    for extra in ("cnv:infercnv", "roi_deg", "arms:tile_deg",
                  "cooccur:leiden_r1.0", "ligrec:leiden_r1.0",
                  "annotation:leiden_r1.0", "plot:volcano:leiden_r1.0"):
        g.upsert(extra, "code", deps=["normalize"])
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)

    sections = _full_sections()
    by_key = {n.key: n for s in sections for n in s.nodes}
    sel = stale_results.select_stale(g, sections)

    assert sel.keys, "the fixture should produce a non-empty selection"
    for key in sel.keys:
        assert key in by_key, f"{key} is not a row in the inventory"
        assert by_key[key].deletable, f"{key} is blocked and must not be selected"


def test_every_matched_key_is_attributed_to_a_stale_node():
    g = _graph()
    g.upsert("normalize", "changed", deps=["preamble"], kind=SETUP)
    sel = stale_results.select_stale(g, _full_sections())
    ids = set(stale_results.stale_ids(g))
    attributed = {k for _, keys in sel.matched for k in keys}
    assert set(sel.keys) == attributed
    assert {nid for nid, _ in sel.matched} <= ids
    assert {nid for nid, _ in sel.unmatched} <= ids


# ── The rename guard ─────────────────────────────────────────────────────────

def test_a_re_emitted_preamble_looks_like_a_path_rewrite():
    g = _graph()
    g.upsert("preamble", 'data_path = Path(r"/new/place")', kind=SETUP)
    assert stale_results.looks_like_a_path_rewrite(g)


def test_one_re_run_step_does_not_look_like_a_path_rewrite():
    g = _graph()
    g.upsert("clustering:leiden_r1.0", "changed", deps=["normalize"])
    assert not stale_results.looks_like_a_path_rewrite(g)


def test_a_fresh_or_missing_graph_is_not_a_path_rewrite():
    assert not stale_results.looks_like_a_path_rewrite(None)
    assert not stale_results.looks_like_a_path_rewrite(_graph())
    assert not stale_results.looks_like_a_path_rewrite(ProvGraph())


# ── The report ───────────────────────────────────────────────────────────────

def test_the_report_names_what_was_spared_and_why():
    g = _graph()
    g.upsert("clustering:leiden_r1.0", "changed", deps=["normalize"])
    sel = stale_results.select_stale(g, _full_sections())
    text = stale_results.describe_selection(sel, plots_note="figures stay in plots/")
    assert "nhood_enrichment" in text
    assert "nhood:leiden_r0.5" in text          # the fresh step holding it
    assert "figures stay in plots/" in text


def test_the_report_says_so_when_there_is_nothing_stale():
    text = stale_results.describe_selection(stale_results.StaleSelection())
    assert "No stale steps" in text


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
