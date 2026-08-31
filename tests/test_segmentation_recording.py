"""Swapping the cell segmentation must reach the recorded notebook.

Tools → Segmentation replaces the Xenium cells with a custom set
(``scripts/extract_seurat_segmentation.R`` → ``build_custom_segmentation.py``),
rebinding ``ctx.adata`` and clearing every derived result in the GUI. It
recorded nothing at all, so the preamble kept saying ``adata = sdata["table"]``
and every node recorded afterwards claimed to be about the Xenium cells.

That is the failure mode worth a test of its own: a replay of such a notebook
does not error. It runs the whole analysis against a different set of cells and
reports success, so ``allow_errors=False`` — and ``verify_notebook.py`` with it
— sees nothing wrong. The numbers come out; they are about something else.

The swap is recorded by upserting ``preamble``, because the preamble is the node
that says where ``adata`` comes from. Flagging the earlier results stale is a
consequence of that, and is correct: they were computed on cells that are no
longer bound.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("numpy")
pytest.importorskip("pandas")
pytest.importorskip("anndata")
pytest.importorskip("qtpy")

TAB = (Path(__file__).resolve().parent.parent
       / "src" / "palms" / "tabs" / "tab_segmentation.py")


def _ctx(data_path: Path, segmentation: str = "xenium"):
    from palms.tabs._helpers import create_shared_helpers
    from palms.utils.viewer_context import ViewerContext

    context = ViewerContext(
        data_path=data_path,
        state={
            "record_code": True, "code_journal": [],
            "prov_graph_restored": True,
            "segmentation_source": segmentation,
        },
    )
    create_shared_helpers(context)
    return context


def _preamble_code(ctx) -> str:
    ctx.record_preamble()
    return ctx.state["prov_graph"].get("preamble").code


def _raw_dataset(tmp_path: Path) -> Path:
    # One marker is enough for has_raw_xenium_source; see test_preamble_recording.
    (tmp_path / "cell_feature_matrix.h5").write_bytes(b"")
    return tmp_path


def test_native_segmentation_binds_the_xenium_table(tmp_path, qapp):
    """The unchanged case: no swap, no mention of a custom table."""
    code = _preamble_code(_ctx(_raw_dataset(tmp_path)))

    assert 'adata = sdata["table"].copy()' in code
    assert "custom_table" not in code


def test_a_custom_segmentation_binds_the_cached_custom_table(tmp_path, qapp):
    """The custom cells are not in the 10x output, so ``xenium()`` cannot reach
    them. They live in the store the tab caches them into."""
    code = _preamble_code(_ctx(_raw_dataset(tmp_path), segmentation="custom"))

    assert 'sdata = xenium(data_path)' in code          # the raw load is intact
    assert ('adata = sd.read_zarr(data_path / "sdata_cached.zarr")'
            '.tables["custom_table"].copy()') in code
    assert 'adata = sdata["table"]' not in code


def test_a_crop_export_reads_the_custom_table_from_the_store_it_opened(
    tmp_path, qapp,
):
    """A crop export already *is* the zarr store, so opening it twice would be
    a second full read of the same thing under a different name."""
    code = _preamble_code(_ctx(tmp_path, segmentation="custom"))

    assert 'adata = sdata.tables["custom_table"].copy()' in code
    assert code.count("read_zarr") == 1


def test_the_swap_flags_results_computed_on_the_other_cells_stale(
    tmp_path, qapp,
):
    """The acceptance case. A clustering recorded before the swap was computed
    on cells that are no longer bound, and the graph has to say so."""
    ctx = _ctx(_raw_dataset(tmp_path))
    ctx.record_preamble()
    ctx.record_node("clustering:leiden_r1.0", "adata.obs['leiden'] = 1",
                    deps=["preamble"])
    graph = ctx.state["prov_graph"]
    assert graph.get("clustering:leiden_r1.0").stale is False

    ctx.state["segmentation_source"] = "custom"
    ctx.record_preamble()

    assert graph.get("clustering:leiden_r1.0").stale is True
    assert "custom_table" in graph.get("preamble").code


def test_reverting_puts_the_xenium_table_back(tmp_path, qapp):
    ctx = _ctx(_raw_dataset(tmp_path), segmentation="custom")
    ctx.record_preamble()

    ctx.state["segmentation_source"] = "xenium"
    ctx.record_preamble()

    code = ctx.state["prov_graph"].get("preamble").code
    assert 'adata = sdata["table"].copy()' in code
    assert "custom_table" not in code


def test_re_recording_an_unchanged_segmentation_flags_nothing(tmp_path, qapp):
    """Every launch re-emits the preamble, and Tools → Segmentation's restore
    handler re-applies the custom segmentation on top of it. If that round trip
    changed the recorded code, opening the dataset would mark the whole notebook
    ⚠ for nothing — which is what a manual dataset rename used to do."""
    ctx = _ctx(_raw_dataset(tmp_path), segmentation="custom")
    ctx.record_preamble()
    ctx.record_node("clustering:leiden_r1.0", "adata.obs['leiden'] = 1",
                    deps=["preamble"])

    ctx.record_preamble()          # the restore handler's re-emit

    assert ctx.state["prov_graph"].get("clustering:leiden_r1.0").stale is False


def _calls_in(func_name: str) -> set[str]:
    """Names of the functions called inside *func_name* in tab_segmentation."""
    tree = ast.parse(TAB.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == func_name:
            return {
                c.func.id for c in ast.walk(node)
                if isinstance(c, ast.Call) and isinstance(c.func, ast.Name)
            }
    raise AssertionError(f"{func_name} is gone from tab_segmentation.py")


@pytest.mark.parametrize("func", ["_apply_custom_segmentation",
                                  "_revert_xenium_segmentation"])
def test_both_swap_paths_record(func):
    """Source guard: a swap that records nothing is the defect this file is
    about, and it is invisible at runtime — both directions have to record."""
    assert "_record_segmentation" in _calls_in(func)


def test_a_failed_restore_puts_the_flag_back():
    """``app.py`` seeds the flag from the session *before* the restore runs, so
    a restore that does not load the custom cells leaves the preamble claiming a
    table that is not bound. Both failure branches have to correct it."""
    assert "_restore_failed" in _calls_in("_restore_session")
    src = TAB.read_text()
    body = src[src.index("def _restore_failed"):]
    assert 'ctx.state["segmentation_source"] = "xenium"' in body
