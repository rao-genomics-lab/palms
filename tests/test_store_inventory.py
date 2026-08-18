"""The dataset inventory, and the predicate that keeps the raw output safe.

``assert_deletable`` is the whole of the promise that the Dataset tab can never
touch the original 10x output. The important test here is not any single refusal
but :func:`test_every_deletable_node_passes_assert_deletable`: it asserts the
*property* over every node the inventory produces, so the guarantee is
structural rather than a list of cases somebody remembered to write down.
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace  # noqa: F401  (kept for parity with sibling tests)

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")

from xenium_viewer.utils import store_inventory as si  # noqa: E402
from xenium_viewer.utils import zarr_safe  # noqa: E402


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


@pytest.fixture
def dataset(tiny_sdata, make_table):
    """A dataset directory with raw output, viewer data, sidecars and a backup.

    Built on ``tiny_sdata``, whose store is ``tmp_path/sdata_cached.zarr`` — so
    the dataset directory is ``tmp_path`` with no extra plumbing.
    """
    from spatialdata.models import Image2DModel

    cache = Path(tiny_sdata.path)
    data_path = cache.parent

    # ── raw 10x output (bytes only; no real formats are needed here) ──
    (data_path / "experiment.xenium").write_text('{"run_name": "test"}')
    (data_path / "transcripts.parquet").write_bytes(b"PAR1" * 64)
    (data_path / "cells.zarr.zip").write_bytes(b"PK\x03\x04" * 16)
    (data_path / "morphology_focus").mkdir()
    (data_path / "morphology_focus" / "morphology_focus_0000.ome.tif").write_bytes(b"II*\0")
    clustering = data_path / "analysis" / "clustering" / "gene_expression_graphclust"
    clustering.mkdir(parents=True)
    (clustering / "clusters.csv").write_text("Barcode,Cluster\na,1\n")

    # ── a deletable element plus the landmarks it cascades into ──
    image = Image2DModel.parse(np.zeros((3, 4, 4), dtype="uint8"))
    zarr_safe.safe_write_element(tiny_sdata, "ext_slide2", image)

    import geopandas as gpd
    from shapely.geometry import Point
    from spatialdata.models import ShapesModel
    landmarks = ShapesModel.parse(
        gpd.GeoDataFrame({"geometry": [Point(1, 1), Point(2, 2)], "radius": [1.0, 1.0]})
    )
    zarr_safe.safe_write_element(tiny_sdata, "ext_slide2_xenium_lm", landmarks)

    # ── table contents, persisted once ──
    adata = tiny_sdata["table"]
    n = adata.n_obs
    adata.obs["clustering_leiden_r1.0"] = pd.Categorical(["1", "2"] * (n // 2))
    adata.obs["cluster_labels_leiden_r1.0"] = pd.Categorical(["A", "B"] * (n // 2))
    adata.uns["rank_genes_groupby"] = "clustering_leiden_r1.0"
    adata.obsm["X_umap"] = np.zeros((n, 2), dtype="float32")
    # A second write of the table also seeds .xv_trash for free.
    zarr_safe.safe_write_element(tiny_sdata, "table", adata)

    # ── sidecars, the preprocess cache and a whole-cache backup ──
    sidecars = data_path / "viewer_cache"
    sidecars.mkdir()
    (sidecars / "adata_norm_cache.h5ad").write_bytes(b"\x89HDF" * 32)
    (sidecars / "prov_graph.json").write_text("[]")
    (sidecars / "cnv_infercnv_result.json").write_text("{}")
    transcripts = data_path / "transcript_cache"
    transcripts.mkdir()
    (transcripts / "ACTA2.feather").write_bytes(b"ARROW1" * 8)
    (transcripts / "MYH11.feather").write_bytes(b"ARROW1" * 8)
    backup = data_path / "sdata_cached_backup_20260101_000000.zarr"
    backup.mkdir()
    (backup / "zarr.json").write_text("{}")

    return SimpleNamespace(data_path=data_path, cache=cache, sdata=tiny_sdata)


def _all_nodes(sections):
    return [n for section in sections for n in section.nodes]


def _by_key(sections):
    return {n.key: n for n in _all_nodes(sections)}


# ── sections and the raw output ───────────────────────────────────────────────

def test_all_five_sections_are_present_even_without_a_cache(tmp_path):
    (tmp_path / "experiment.xenium").write_text("{}")
    (tmp_path / "cells.zarr.zip").write_bytes(b"PK")
    sections = si.build_inventory(tmp_path, None)
    assert [s.title for s in sections] == [
        "Original Xenium output", "Viewer cache", "Session state",
        "Derived caches", "Backups & trash",
    ]
    for section in sections[1:]:
        if not section.nodes:
            assert section.note, f"{section.title} must explain why it is empty"


def test_a_crop_export_is_not_described_as_untouched_10x_output(tmp_path):
    """Its experiment.xenium and transcripts.parquet were written by the viewer.

    Still not deletable — they are the only copy — but the old wording claimed
    the viewer never modifies them, which is the same wrong mental model that
    let Force Rebuild loose on these datasets.
    """
    (tmp_path / "experiment.xenium").write_text("{}")
    (tmp_path / "transcripts.parquet").write_bytes(b"PAR1")

    raw = si.build_inventory(tmp_path, None)[0]
    assert raw.title == "Dataset source files"
    assert raw.nodes
    for node in raw.nodes:
        assert node.deletable is False
        assert "only copy" in node.blocked_reason


def test_raw_files_are_listed_and_none_of_them_is_deletable(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    raw = [n for n in _all_nodes(sections) if n.kind == si.RAW]
    names = {n.name for n in raw}
    assert {"experiment.xenium", "transcripts.parquet", "cells.zarr.zip",
            "morphology_focus", "analysis"} <= names
    for node in raw:
        assert node.deletable is False
        assert node.blocked_reason


def test_an_unknown_top_level_entry_is_classified_raw(dataset):
    (dataset.data_path / "weird_vendor_file.bin").write_bytes(b"\x00" * 8)
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "raw:weird_vendor_file.bin"]
    assert node.kind == si.RAW and node.deletable is False


def test_the_viewer_directories_are_not_in_the_raw_section(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    raw_names = {n.name for n in _all_nodes(sections) if n.kind == si.RAW}
    assert not raw_names & {
        "sdata_cached.zarr", "viewer_cache", "transcript_cache",
        "sdata_cached_backup_20260101_000000.zarr",
    }


def test_build_inventory_creates_nothing_and_changes_nothing(dataset):
    import shutil
    # Remove viewer_cache first, so a stray sidecar_dir(create=True) is caught.
    shutil.rmtree(dataset.data_path / "viewer_cache")
    before = sorted(p.relative_to(dataset.data_path)
                    for p in dataset.data_path.rglob("*"))
    si.build_inventory(dataset.data_path, dataset.cache)
    after = sorted(p.relative_to(dataset.data_path)
                   for p in dataset.data_path.rglob("*"))
    assert before == after


# ── elements ─────────────────────────────────────────────────────────────────

def test_core_elements_is_derived_from_the_reason_map():
    assert set(si.CORE_ELEMENTS) == set(si._CORE_REASONS)


def test_every_core_element_present_is_listed_sized_and_blocked(dataset):
    nodes = _by_key(si.build_inventory(dataset.data_path, dataset.cache))
    table = nodes["element:tables/table"]
    assert table.deletable is False
    assert table.blocked_reason
    assert table.size_bytes and table.size_bytes > 0


def test_a_user_element_is_deletable_and_recoverable_from_trash(dataset):
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "element:images/ext_slide2"]
    assert node.deletable is True
    assert node.recoverable == si.RECOVER_TRASH
    assert node.blocked_reason == ""


def test_an_element_too_big_for_the_trash_budget_is_unrecoverable(dataset, monkeypatch):
    monkeypatch.setattr(zarr_safe, "DEFAULT_MAX_TRASH_BYTES", 1)
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "element:images/ext_slide2"]
    assert node.deletable is True
    assert node.recoverable == si.RECOVER_NONE


def test_an_unrecognised_element_is_not_deletable(dataset):
    stray = dataset.cache / "shapes" / "mystery"
    stray.mkdir()
    (stray / "zarr.json").write_text('{"attributes": {"encoding-type": "x"}}')
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "element:shapes/mystery"]
    assert node.deletable is False
    assert node.blocked_reason == si._UNRECOGNISED


def test_every_parent_resolves_and_elements_hang_off_a_group(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    keys = set(_by_key(sections))
    for node in _all_nodes(sections):
        assert node.parent == "" or node.parent in keys, node.key
        if node.kind == si.ELEMENT:
            assert node.parent.startswith("group:")


# ── table contents ───────────────────────────────────────────────────────────

def test_obs_clusterings_are_deletable_and_carry_a_cluster_count(dataset):
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "obs:table/clustering_leiden_r1.0"]
    assert node.deletable is True
    assert node.detail == "2 clusters"
    assert node.parent == "element:tables/table"


def test_the_bare_twin_a_leiden_run_leaves_is_deletable_with_its_clustering(dataset):
    """A clustering is two obs columns, so deleting one must take both.

    The recorded step writes adata.obs[$key]; save_clustering_to_adata writes
    clustering_<key>. Offering only the prefixed one left an identical copy.
    """
    from xenium_viewer.utils import zarr_safe as zs
    adata = dataset.sdata["table"]
    adata.obs["leiden_igraph_r1.0"] = adata.obs["clustering_leiden_r1.0"]
    adata.obs["clustering_leiden_igraph_r1.0"] = adata.obs["clustering_leiden_r1.0"]
    zs.safe_write_element(dataset.sdata, "table", adata)

    sections = si.build_inventory(dataset.data_path, dataset.cache)
    twin = _by_key(sections)["obs:table/leiden_igraph_r1.0"]
    assert twin.deletable is True
    assert "clustering_leiden_igraph_r1.0" in twin.detail

    plan = si.plan_deletion(sections, ["obs:table/clustering_leiden_igraph_r1.0"])
    assert "obs:table/leiden_igraph_r1.0" in {n.key for n in plan.nodes}
    assert "obs:table/leiden_igraph_r1.0" in {n.key for n in plan.added}


def test_a_bare_column_with_no_clustering_twin_stays_blocked(dataset):
    """The pairing must not make ordinary Xenium columns selectable."""
    nodes = _by_key(si.build_inventory(dataset.data_path, dataset.cache))
    assert nodes["obs:table/region"].deletable is False
    assert nodes["obs:table/marker"].deletable is False


def test_a_structural_column_is_never_paired_even_if_named_like_one(dataset):
    from xenium_viewer.utils import zarr_safe as zs
    adata = dataset.sdata["table"]
    # A user who names a clustering "region" must not lose the real column.
    adata.obs["clustering_region"] = adata.obs["clustering_leiden_r1.0"]
    zs.safe_write_element(dataset.sdata, "table", adata)
    nodes = _by_key(si.build_inventory(dataset.data_path, dataset.cache))
    assert nodes["obs:table/region"].deletable is False


def test_a_structural_obs_column_is_not_deletable(dataset):
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "obs:table/region"]
    assert node.deletable is False


def test_table_content_nodes_carry_no_path(dataset):
    """So no executor can ever unlink one — deleting a column is a table rewrite."""
    for node in _all_nodes(si.build_inventory(dataset.data_path, dataset.cache)):
        if node.kind in (si.OBS, si.UNS, si.OBSM):
            assert node.path is None, node.key


def test_obs_is_enumerated_from_disk_not_from_memory(dataset):
    dataset.sdata["table"].obs["clustering_not_saved"] = pd.Categorical(
        ["x"] * dataset.sdata["table"].n_obs)
    keys = set(_by_key(si.build_inventory(dataset.data_path, dataset.cache)))
    assert "obs:table/clustering_not_saved" not in keys
    assert "obs:table/clustering_leiden_r1.0" in keys

    del dataset.sdata["table"].obs["clustering_leiden_r1.0"]
    keys = set(_by_key(si.build_inventory(dataset.data_path, dataset.cache)))
    assert "obs:table/clustering_leiden_r1.0" in keys, (
        "a column still on disk must stay visible — that is what there is to delete")


def test_x_umap_is_deletable_but_spatial_is_not(dataset):
    nodes = _by_key(si.build_inventory(dataset.data_path, dataset.cache))
    assert nodes["obsm:table/X_umap"].deletable is True
    if "obsm:table/spatial" in nodes:
        assert nodes["obsm:table/spatial"].deletable is False


def test_zarr_metadata_is_never_offered_as_content(dataset):
    """A group's zarr.json *is* the group; listing it invites deleting it."""
    keys = set(_by_key(si.build_inventory(dataset.data_path, dataset.cache)))
    assert not [k for k in keys if k.endswith("/zarr.json")], sorted(
        k for k in keys if k.endswith("/zarr.json"))
    assert "session:file/zarr.json" not in keys


def test_rank_genes_groupby_is_a_deletable_uns_key(dataset):
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "uns:table/rank_genes_groupby"]
    assert node.deletable is True


# ── derived caches ───────────────────────────────────────────────────────────

def test_the_transcript_cache_is_deletable_and_says_what_it_costs(dataset):
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "derived:transcript_cache"]
    assert node.deletable is True
    assert "xenium-preprocess" in node.detail
    assert node.size_bytes and node.size_bytes > 0


def test_the_provenance_sidecar_is_blocked(dataset):
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "sidecar:prov_graph.json"]
    assert node.deletable is False
    assert "provenance" in node.blocked_reason


def test_a_dated_provenance_backup_is_blocked_too(dataset):
    """They exist precisely to survive a graph going wrong."""
    (dataset.data_path / "viewer_cache"
     / "prov_graph.backup_20260729_1407.json").write_text("[]")
    node = _by_key(si.build_inventory(dataset.data_path, dataset.cache))[
        "sidecar:prov_graph.backup_20260729_1407.json"]
    assert node.deletable is False


@pytest.mark.parametrize("name", [
    "analysis.py", "analysis_notebook.ipynb", "plots", "xenium_viewer.log"])
def test_viewer_output_in_the_dataset_folder_is_not_called_raw(dataset, name):
    """It is outside every deletable root, but it is not 10x's either — saying
    "the viewer never modifies it" about the viewer's own log is just false."""
    target = dataset.data_path / name
    if name == "plots":
        target.mkdir()
        (target / "dotplot.svg").write_text("<svg/>")
    else:
        target.write_text("x\n")
    nodes = _by_key(si.build_inventory(dataset.data_path, dataset.cache))
    assert f"raw:{name}" not in nodes
    node = nodes[f"derived:{name}"]
    assert node.deletable is False
    assert "by hand" in node.blocked_reason


# ── backups and trash ────────────────────────────────────────────────────────

def test_the_backup_and_trash_are_listed_and_deletable(dataset):
    nodes = _by_key(si.build_inventory(dataset.data_path, dataset.cache))
    backup = nodes["backup:sdata_cached_backup_20260101_000000.zarr"]
    assert backup.deletable is True and backup.recoverable == si.RECOVER_NONE
    assert nodes["trash:all"].deletable is True
    assert any(k.startswith("trash:tables/") for k in nodes)


# ── the safety predicate ─────────────────────────────────────────────────────

def test_deletable_roots_are_the_four_viewer_directories(dataset):
    roots = si.deletable_roots(dataset.data_path, dataset.cache)
    assert {r.name for r in roots} == {
        "sdata_cached.zarr", "viewer_cache", "transcript_cache",
        "sdata_cached_backup_20260101_000000.zarr",
    }
    assert dataset.data_path.resolve() not in roots


def test_deletable_roots_tolerates_no_cache_and_missing_directories(tmp_path):
    assert si.deletable_roots(tmp_path, None) == ()
    assert si.deletable_roots(tmp_path / "nope", None) == ()
    (tmp_path / "sdata_cached_prev_1.zarr").mkdir()
    assert [r.name for r in si.deletable_roots(tmp_path, None)] == [
        "sdata_cached_prev_1.zarr"]


@pytest.mark.parametrize("bad", ["self", "parent", "root"])
def test_deletable_roots_refuses_a_root_that_contains_the_dataset(dataset, bad):
    candidate = {"self": dataset.data_path,
                 "parent": dataset.data_path.parent,
                 "root": Path("/")}[bad]
    roots = si.deletable_roots(dataset.data_path, candidate)
    assert candidate.resolve() not in roots


@pytest.mark.parametrize("relative", [
    ".", "transcripts.parquet", "experiment.xenium", "morphology_focus",
    "analysis/clustering/gene_expression_graphclust/clusters.csv",
    "sdata_cached.zarr/../transcripts.parquet",
])
def test_assert_deletable_refuses_the_raw_output(dataset, relative):
    roots = si.deletable_roots(dataset.data_path, dataset.cache)
    with pytest.raises(si.NotDeletable):
        si.assert_deletable(dataset.data_path / relative, roots)


@pytest.mark.parametrize("absolute", ["/", "/etc/passwd"])
def test_assert_deletable_refuses_paths_outside_the_dataset(dataset, absolute):
    roots = si.deletable_roots(dataset.data_path, dataset.cache)
    with pytest.raises(si.NotDeletable):
        si.assert_deletable(Path(absolute), roots)


def test_assert_deletable_refuses_a_symlink_pointing_at_a_raw_file(dataset):
    link = dataset.cache / "shapes" / "sneaky"
    link.symlink_to(dataset.data_path / "transcripts.parquet")
    roots = si.deletable_roots(dataset.data_path, dataset.cache)
    with pytest.raises(si.NotDeletable):
        si.assert_deletable(link, roots)


def test_assert_deletable_accepts_an_element_a_sidecar_and_a_backup(dataset):
    roots = si.deletable_roots(dataset.data_path, dataset.cache)
    si.assert_deletable(dataset.cache / "images" / "ext_slide2", roots, kind=si.ELEMENT)
    si.assert_deletable(dataset.data_path / "viewer_cache" / "adata_norm_cache.h5ad",
                        roots, kind=si.SIDECAR)
    si.assert_deletable(dataset.data_path / "sdata_cached_backup_20260101_000000.zarr",
                        roots, kind=si.BACKUP)
    si.assert_deletable(dataset.data_path / "transcript_cache", roots, kind=si.DERIVED)


def test_assert_deletable_refuses_a_root_itself_unless_the_kind_allows_it(dataset):
    roots = si.deletable_roots(dataset.data_path, dataset.cache)
    with pytest.raises(si.NotDeletable):
        si.assert_deletable(dataset.cache, roots, kind=si.ELEMENT)
    with pytest.raises(si.NotDeletable):
        si.assert_deletable(dataset.data_path / "viewer_cache", roots, kind=si.SIDECAR)


def test_assert_deletable_with_no_roots_refuses_everything(dataset):
    with pytest.raises(si.NotDeletable):
        si.assert_deletable(dataset.cache / "images" / "ext_slide2", ())


def test_every_deletable_node_passes_assert_deletable(dataset):
    """The property that makes constraint 1 structural, not a promise."""
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    roots = si.deletable_roots(dataset.data_path, dataset.cache)
    checked = 0
    for node in _all_nodes(sections):
        if not node.deletable:
            continue
        checked += 1
        if node.path is None:
            assert (node.kind in si._PATHLESS_KINDS
                    or node.key.startswith("session:attr/")), node.key
            continue
        si.assert_node_deletable(node, roots)
    assert checked > 5, "the fixture must produce deletable nodes to check"


def test_assert_node_deletable_refuses_every_blocked_node(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    roots = si.deletable_roots(dataset.data_path, dataset.cache)
    for node in _all_nodes(sections):
        if node.deletable:
            continue
        with pytest.raises(si.NotDeletable):
            si.assert_node_deletable(node, roots)


def test_the_property_holds_when_the_cache_is_reached_through_a_symlink(dataset):
    link = dataset.data_path / "linked_cache.zarr"
    link.symlink_to(dataset.cache)
    sections = si.build_inventory(dataset.data_path, link)
    roots = si.deletable_roots(dataset.data_path, link)
    for node in _all_nodes(sections):
        if node.deletable and node.path is not None:
            si.assert_node_deletable(node, roots)


# ── planning ─────────────────────────────────────────────────────────────────

def test_plan_expands_the_ext_landmark_cascade(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    plan = si.plan_deletion(sections, ["element:images/ext_slide2"])
    keys = {n.key for n in plan.nodes}
    assert "element:shapes/ext_slide2_xenium_lm" in keys
    assert "element:shapes/ext_slide2_xenium_lm" in {n.key for n in plan.added}


def test_plan_expands_clustering_to_cluster_labels(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    plan = si.plan_deletion(sections, ["obs:table/clustering_leiden_r1.0"])
    assert {n.key for n in plan.nodes} == {
        "obs:table/clustering_leiden_r1.0", "obs:table/cluster_labels_leiden_r1.0"}


def test_plan_does_not_cascade_into_a_node_that_does_not_exist(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    node = _by_key(sections)["element:images/ext_slide2"]
    # Only the _xenium_lm half exists in the fixture.
    assert "element:shapes/ext_slide2_image_lm" not in node.cascade


def test_plan_dedups_a_parent_and_its_child(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    child = next(n.key for n in _all_nodes(sections)
                 if n.key.startswith("trash:tables/"))
    plan = si.plan_deletion(sections, ["trash:all", child])
    assert [n.key for n in plan.nodes] == ["trash:all"]


def test_plan_orders_table_edits_first_and_backups_last(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    plan = si.plan_deletion(sections, [
        "backup:sdata_cached_backup_20260101_000000.zarr",
        "element:images/ext_slide2",
        "obs:table/clustering_leiden_r1.0",
        "sidecar:adata_norm_cache.h5ad",
    ])
    ranks = [si._KIND_ORDER[n.kind] for n in plan.nodes]
    assert ranks == sorted(ranks)
    assert plan.nodes[0].kind == si.OBS
    assert plan.nodes[-1].kind == si.BACKUP


def test_plan_refuses_a_blocked_key(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    with pytest.raises(si.NotDeletable):
        si.plan_deletion(sections, ["element:tables/table"])


def test_plan_drops_an_unknown_key_with_a_warning(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    plan = si.plan_deletion(sections, ["element:images/gone_already"])
    assert plan.is_empty
    assert plan.dropped == ("element:images/gone_already",)
    assert any("no longer present" in w for w in plan.warnings)


def test_plan_totals_and_names_what_is_unrecoverable(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    plan = si.plan_deletion(sections, [
        "sidecar:adata_norm_cache.h5ad", "element:images/ext_slide2"])
    assert plan.total_bytes == sum(n.size_bytes or 0 for n in plan.nodes)
    unrecoverable = {n.key for n in plan.unrecoverable}
    assert "sidecar:adata_norm_cache.h5ad" in unrecoverable
    assert "element:images/ext_slide2" not in unrecoverable


def test_describe_plan_names_every_path_the_total_and_the_warning(dataset):
    sections = si.build_inventory(dataset.data_path, dataset.cache)
    plan = si.plan_deletion(sections, [
        "obs:table/clustering_leiden_r1.0", "sidecar:adata_norm_cache.h5ad"])
    text = si.describe_plan(plan)
    assert str(dataset.data_path / "viewer_cache" / "adata_norm_cache.h5ad") in text
    assert "Not recoverable" in text
    assert "provenance graph keeps its clustering step" in text
    assert si.human_bytes(plan.total_bytes) in text


def test_describe_plan_of_an_empty_plan_says_so():
    assert si.describe_plan(si.Plan()) == "Nothing selected."


# ── robustness ───────────────────────────────────────────────────────────────

def test_inventory_of_a_store_with_a_broken_root_zarr_json_still_lists_elements(dataset):
    (dataset.cache / "zarr.json").write_text("{")
    keys = set(_by_key(si.build_inventory(dataset.data_path, dataset.cache)))
    assert "element:tables/table" in keys
    assert "element:images/ext_slide2" in keys


def test_inventory_survives_a_missing_cache_directory(tmp_path):
    sections = si.build_inventory(tmp_path, tmp_path / "sdata_cached.zarr")
    cache = next(s for s in sections if s.title == "Viewer cache")
    assert cache.nodes == () and cache.note


def test_session_section_names_its_in_memory_mirror():
    """The session trap, pinned at model level: deleting the disk copy of a
    session group without clearing ctx.he_state puts it back on the next save."""
    assert ("he_state", "affine_3x3") in si.session_memory_keys("session:group/he")
    assert ("arms_state", "affine_3x3") in si.session_memory_keys("session:group/arms")
    assert si.session_memory_keys("element:images/ext_slide2") == ()


# ── source guard ─────────────────────────────────────────────────────────────

def test_store_inventory_never_mutates_the_filesystem():
    """An inventory that can write is a defect by definition."""
    source = Path(si.__file__).read_text()
    forbidden = re.compile(
        r"\b(?:os\.remove|os\.rename|os\.rmdir|shutil\.(?:move|rmtree)|"
        r"\.unlink\(|\.mkdir\(|\.write_text\(|\.write_bytes\(|create=True)")
    offenders = [
        f"{i}: {line.strip()}"
        for i, line in enumerate(source.splitlines(), 1)
        if forbidden.search(line) and not line.lstrip().startswith("#")
    ]
    assert offenders == [], (
        "store_inventory must only read: " + ", ".join(offenders))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
