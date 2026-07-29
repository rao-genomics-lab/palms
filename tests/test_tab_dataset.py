"""The Dataset tab: it builds, and it never reaches outside the viewer's dirs.

The model and the safety predicate are covered by ``test_store_inventory.py``.
What is left for here is the Qt shell and the executor — in particular the two
things that have bitten this codebase before: forgetting that a deleted obs
column also lives in the in-memory ``ctx.clusterings`` dict, and deleting a
session group on disk without clearing the state ``save_session`` rebuilds it
from.
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
pytest.importorskip("qtpy")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")

from xenium_viewer.tabs import tab_dataset  # noqa: E402
from xenium_viewer.utils import store_inventory as si  # noqa: E402


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _ctx(sdata=None, data_path=None, no_cache=False):
    return SimpleNamespace(
        sdata=sdata, adata=(sdata["table"] if sdata is not None else None),
        viewer=None, state={}, no_cache=no_cache, segmentation_source="xenium",
        data_path=data_path or (Path(sdata.path).parent if sdata else None),
        clusterings={}, he_state={}, arms_state={},
        external_images_state=[], patch_overlays_state=[],
        refresh_clustering_choices=None, reload_dataset=None,
    )


# ── the widget ───────────────────────────────────────────────────────────────

def test_build_tab_returns_a_widget_and_exports(qapp, tiny_sdata):
    widget, exports = tab_dataset.build_tab(_ctx(tiny_sdata))
    assert widget is not None
    assert callable(exports["restore_session"])


def test_build_tab_survives_no_cache_mode(qapp, tiny_sdata, tmp_path):
    widget, exports = tab_dataset.build_tab(
        _ctx(None, data_path=tmp_path, no_cache=True))
    assert widget is not None


def test_build_tab_survives_a_dataset_with_no_directory(qapp):
    widget, _ = tab_dataset.build_tab(_ctx(None, data_path=None))
    assert widget is not None


def test_cache_path_falls_back_to_the_conventional_name(tiny_sdata, tmp_path):
    ctx = _ctx(tiny_sdata)
    assert tab_dataset._cache_path(ctx) == Path(tiny_sdata.path)
    # The store exists beside tiny_sdata, so the fallback finds it even with no
    # live sdata; a directory without one yields None.
    assert tab_dataset._cache_path(_ctx(None, data_path=tmp_path)) == Path(tiny_sdata.path)
    empty = tmp_path / "elsewhere"
    empty.mkdir()
    assert tab_dataset._cache_path(_ctx(None, data_path=empty)) is None


def test_the_tab_does_not_scan_at_build_time(qapp, tiny_sdata):
    """Walking the whole dataset directory must not be charged to every launch."""
    ctx = _ctx(tiny_sdata)
    tab_dataset.build_tab(ctx)
    assert "_dataset_sections" not in ctx.state
    assert "_dataset_worker" not in ctx.state


# ── the executor ─────────────────────────────────────────────────────────────

@pytest.fixture
def loaded(tiny_sdata):
    """A ctx over a store with a deletable element, obs columns and sidecars."""
    from spatialdata.models import Image2DModel
    from xenium_viewer.utils import zarr_safe

    cache = Path(tiny_sdata.path)
    data_path = cache.parent
    zarr_safe.safe_write_element(
        tiny_sdata, "ext_slide2",
        Image2DModel.parse(np.zeros((3, 4, 4), dtype="uint8")))

    adata = tiny_sdata["table"]
    n = adata.n_obs
    adata.obs["clustering_leiden_r1.0"] = pd.Categorical(["1", "2"] * (n // 2))
    adata.obs["cluster_labels_leiden_r1.0"] = pd.Categorical(["A", "B"] * (n // 2))
    adata.uns["rank_genes_groupby"] = "clustering_leiden_r1.0"
    zarr_safe.safe_write_element(tiny_sdata, "table", adata)

    sidecars = data_path / "viewer_cache"
    sidecars.mkdir(exist_ok=True)
    (sidecars / "adata_norm_cache.h5ad").write_bytes(b"\x89HDF" * 32)

    ctx = _ctx(tiny_sdata)
    ctx.clusterings = {"leiden_r1.0": adata.obs["clustering_leiden_r1.0"]}
    ctx.state["custom_clusterings"] = {"leiden_r1.0": object()}
    ctx.state["cluster_labels"] = {"leiden_r1.0": {"1": "A"}}
    return ctx


def _plan(ctx, keys):
    sections = si.build_inventory(
        Path(ctx.data_path), tab_dataset._cache_path(ctx))
    return si.plan_deletion(sections, keys)


def test_deleting_an_element_removes_it_and_keeps_a_trash_copy(loaded):
    from spatialdata import read_zarr
    plan = _plan(loaded, ["element:images/ext_slide2"])
    result = tab_dataset._apply_deletion(loaded, plan)

    assert result.failed == [], result.failed
    assert result.removed == ["ext_slide2"]
    assert result.needs_reload is True
    assert not (Path(loaded.sdata.path) / "images" / "ext_slide2").exists()
    trash = Path(loaded.sdata.path) / ".xv_trash" / "images"
    assert any(p.name.startswith("ext_slide2") for p in trash.iterdir())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert "ext_slide2" not in read_zarr(loaded.sdata.path).images


def test_deleting_two_obs_columns_persists_the_table_exactly_once(loaded, monkeypatch):
    calls = []
    monkeypatch.setattr(tab_dataset, "_persist_table",
                        lambda ctx: calls.append(ctx))
    plan = _plan(loaded, ["obs:table/clustering_leiden_r1.0"])   # cascades to labels
    result = tab_dataset._apply_deletion(loaded, plan)

    assert len(result.removed) == 2
    assert len(calls) == 1, "one rewrite for the batch, not one per column"
    assert "clustering_leiden_r1.0" not in loaded.adata.obs.columns
    assert "cluster_labels_leiden_r1.0" not in loaded.adata.obs.columns


def test_deleting_a_clustering_also_forgets_it_in_memory(loaded):
    """refresh_clustering_choices reads ctx.clusterings, not adata.obs."""
    refreshed = []
    loaded.refresh_clustering_choices = lambda: refreshed.append(True)
    plan = _plan(loaded, ["obs:table/clustering_leiden_r1.0"])
    tab_dataset._apply_deletion(loaded, plan)

    assert "leiden_r1.0" not in loaded.clusterings
    assert "leiden_r1.0" not in loaded.state["custom_clusterings"]
    assert "leiden_r1.0" not in loaded.state["cluster_labels"]
    assert refreshed == [True]


def test_deleting_a_uns_key_removes_it_from_the_table(loaded):
    plan = _plan(loaded, ["uns:table/rank_genes_groupby"])
    result = tab_dataset._apply_deletion(loaded, plan)
    assert result.failed == [] and result.removed == ["rank_genes_groupby"]
    assert "rank_genes_groupby" not in loaded.adata.uns


def test_deleting_a_sidecar_unlinks_it(loaded):
    plan = _plan(loaded, ["sidecar:adata_norm_cache.h5ad"])
    result = tab_dataset._apply_deletion(loaded, plan)
    assert result.failed == []
    assert result.bytes_freed > 0
    assert not (Path(loaded.data_path) / "viewer_cache"
                / "adata_norm_cache.h5ad").exists()


def test_deleting_a_session_group_also_clears_the_in_memory_mirror(loaded):
    """Otherwise the next save_session writes the affine straight back."""
    import numpy as _np
    from xenium_viewer.utils import zarr_safe

    cache = Path(loaded.sdata.path)
    with zarr_safe.safe_group_update(cache, "viewer_session") as (group, _stage):
        he = group.create_group("he")
        arr = he.create_array("affine_3x3", shape=(3, 3), dtype="float64")
        arr[:] = _np.eye(3)
        group.attrs.update({"he_filename": "slide.tif"})
    loaded.he_state.update({"affine_3x3": _np.eye(3), "he_filename": "slide.tif"})

    plan = _plan(loaded, ["session:group/he"])
    result = tab_dataset._apply_deletion(loaded, plan)

    assert result.failed == [], result.failed
    assert not (cache / "viewer_session" / "he").exists()
    assert "affine_3x3" not in loaded.he_state
    assert "he_filename" not in loaded.he_state


def test_a_blocked_node_never_reaches_the_filesystem(loaded):
    """plan_deletion refuses it, so the executor is never even asked."""
    sections = si.build_inventory(
        Path(loaded.data_path), tab_dataset._cache_path(loaded))
    with pytest.raises(si.NotDeletable):
        si.plan_deletion(sections, ["element:tables/table"])
    assert (Path(loaded.sdata.path) / "tables" / "table").exists()


def test_a_node_pointing_outside_every_root_is_refused_by_the_executor(loaded):
    """The choke point, exercised through _apply_deletion rather than directly."""
    raw = Path(loaded.data_path) / "transcripts.parquet"
    raw.write_bytes(b"PAR1")
    forged = si.Node(
        key="sidecar:forged", kind=si.SIDECAR, name="transcripts.parquet",
        path=raw, deletable=True,
    )
    result = tab_dataset._apply_deletion(loaded, si.Plan(nodes=(forged,)))
    assert result.removed == []
    assert result.failed and "not inside a directory" in result.failed[0][1]
    assert raw.exists(), "the raw output must survive a forged node"


def test_one_failing_node_does_not_abort_the_batch(loaded):
    """A partly applied batch must finish and then describe itself."""
    outside = si.Node(                       # refused by the choke point
        key="sidecar:outside", kind=si.SIDECAR, name="outside.bin",
        path=Path(loaded.data_path).parent / "outside.bin", deletable=True,
    )
    plan = _plan(loaded, ["sidecar:adata_norm_cache.h5ad"])
    result = tab_dataset._apply_deletion(
        loaded, si.Plan(nodes=(outside,) + plan.nodes))

    assert [n for n, _ in result.failed] == ["outside.bin"]
    assert result.removed == ["viewer_cache/adata_norm_cache.h5ad"]
    summary = result.summary()
    assert "outside.bin" in summary and "adata_norm_cache.h5ad" in summary


def test_a_path_that_vanished_between_scan_and_delete_is_not_an_error(loaded):
    gone = si.Node(
        key="sidecar:gone", kind=si.SIDECAR, name="gone.h5ad",
        path=Path(loaded.data_path) / "viewer_cache" / "gone.h5ad",
        deletable=True,
    )
    result = tab_dataset._apply_deletion(loaded, si.Plan(nodes=(gone,)))
    assert result.failed == []
    assert result.removed == ["gone.h5ad"]


def test_the_result_summary_names_what_failed(loaded):
    result = tab_dataset.DeletionResult(
        removed=["a"], failed=[("b", "because")], bytes_freed=1024)
    text = result.summary()
    assert "a" in text and "b: because" in text and "1.0 KB" in text


# ── source guards ────────────────────────────────────────────────────────────

def test_the_tab_never_calls_delete_element_from_disk():
    """The unsafe path: it unlinks before a replacement exists."""
    source = Path(tab_dataset.__file__).read_text()
    assert "delete_element_from_disk" not in source


def _enclosing_function(source: str, line: int) -> str:
    """Name of the innermost def containing *line*, or "" at module level."""
    import ast
    best, best_span = "", None
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if node.lineno <= line <= (node.end_lineno or node.lineno):
            span = (node.end_lineno or node.lineno) - node.lineno
            if best_span is None or span < best_span:
                best, best_span = node.name, span
    return best


def test_only_one_function_removes_files_and_its_callers_vet_the_node():
    """Filesystem removal has exactly one home, and nothing reaches it unchecked."""
    source = Path(tab_dataset.__file__).read_text()
    removal = re.compile(r"(shutil\.rmtree|\.unlink\(|os\.remove|shutil\.move)")
    offenders = [
        f"{i}: {line.strip()}"
        for i, line in enumerate(source.splitlines(), 1)
        if removal.search(line) and not line.lstrip().startswith("#")
        and _enclosing_function(source, i) != "_remove_tree"
    ]
    assert offenders == [], (
        "route filesystem removal through _remove_tree: " + ", ".join(offenders))

    callers = {
        _enclosing_function(source, i)
        for i, line in enumerate(source.splitlines(), 1)
        if "_remove_tree(" in line and "def _remove_tree" not in line
    }
    import ast
    bodies = {
        node.name: ast.get_source_segment(source, node) or ""
        for node in ast.walk(ast.parse(source))
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    for caller in callers:
        assert "assert_node_deletable" in bodies.get(caller, ""), (
            f"{caller} removes files without calling assert_node_deletable first")


def test_kind_order_puts_table_edits_first_and_backups_last():
    """Deleting custom_table then re-persisting obs would recreate it."""
    assert si._KIND_ORDER[si.OBS] < si._KIND_ORDER[si.ELEMENT]
    assert si._KIND_ORDER[si.BACKUP] == max(si._KIND_ORDER.values())


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
