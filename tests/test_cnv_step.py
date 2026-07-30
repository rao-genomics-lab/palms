"""The inferCNV step: recorded source == executed source, and it matches the
pipeline it replaces.

inferCNV itself needs a gene-position reference this suite has no way to supply,
so these are structural: they pin the parameters the recorded cell carries and
the divergences the migration closed. The behavioural equivalence that *can* be
checked cheaply is against ``run_cnv_pipeline``'s own source, which the CopyKAT
worker still uses — the two must not drift apart again.

Importing the tab module pulls in Qt/napari, so run headless:
    QT_QPA_PLATFORM=offscreen pytest tests/test_cnv_step.py
"""
from __future__ import annotations

import ast
import inspect
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("anndata")
pytest.importorskip("scanpy")

from xenium_viewer.utils.steps import Step, check_step  # noqa: E402
from xenium_viewer.tabs.tab_cnv import _cnv_template  # noqa: E402


def _step(subset=False):
    params = {
        "reference_clustering": "leiden_r1.0",
        "reference_obs_key": "cnv_reference",
        "reference_categories": ["0"],
        "n_neighbors": 15, "smoothing_neighbors": 20,
        "window_size": 60, "step": 10, "lfc_clip": 4.0, "resolution": 0.2,
    }
    if subset:
        params["include"] = ["0", "1", "2"]
    return Step(id="cnv:infercnv", template=_cnv_template(subset), params=params,
                deps=["clustering:leiden_r1.0"],
                outputs=["adata_cnv", "cnv_cluster_keys", "cnv_clusters", "cnv_score"])


@pytest.mark.parametrize("subset", [False, True])
def test_template_is_valid_python_and_self_contained(subset):
    source = _step(subset).render()
    ast.parse(source)
    assert check_step(_step(subset), available={"sc", "pd", "np", "adata"}) == set()


# ── the three divergences this closes ────────────────────────────────────────

def test_normalisation_records_the_target_sum_it_uses():
    """The old cell said `sc.pp.normalize_total(adata_cnv)` — scanpy's *median*
    default — while the viewer normalised to 1e4 via get_normalized_adata."""
    assert "sc.pp.normalize_total(adata_cnv, target_sum=1e4)" in _cnv_template(False)


def test_lfc_clip_reaches_the_recorded_source():
    """run_infercnv's lfc_clip was applied but never recorded, so a replay used
    infercnvpy's default instead of the pipeline's 4.0."""
    assert "lfc_clip=4.0" in _step().render()


def test_cluster_resolutions_records_dendrogram_false():
    assert "dendrogram=False" in _cnv_template(False)


# ── parameters ───────────────────────────────────────────────────────────────

def test_every_widget_parameter_appears_literally():
    source = _step().render()
    for fragment in ("n_neighbors=15", "smoothing_neighbors=20",
                     "window_size=60", "step=10", "[0.2]",
                     "reference_categories=['0']"):
        assert fragment in source, fragment


def test_cell_type_subset_is_recorded_as_code_not_prose():
    source = _step(subset=True).render()
    assert "isin(['0', '1', '2'])" in source
    assert "isin(" not in _step(subset=False).render()


def test_step_binds_the_names_the_tab_reads_back():
    source = _step().render()
    for name in ("adata_cnv", "cnv_cluster_keys", "cnv_clusters", "cnv_score"):
        assert f"\n{name} " in source or f"\n{name} =" in source or f"{name} = " in source


def test_notebook_does_not_import_the_viewer_package():
    assert "xenium_viewer" not in _cnv_template(True)


# ── the two implementations must not drift apart again ───────────────────────

def test_run_cnv_pipeline_still_normalises_the_way_the_template_does():
    """``run_cnv_pipeline`` is still the CopyKAT path. If its normalisation
    changes, the inferCNV template has to change with it."""
    from xenium_viewer.utils import cnv_analysis, gene_analysis

    assert "get_normalized_adata(adata)" in inspect.getsource(
        cnv_analysis.run_cnv_pipeline)
    norm_src = inspect.getsource(gene_analysis.get_normalized_adata)
    assert "target_sum=1e4" in norm_src
    assert "sc.pp.log1p" in norm_src
    assert "sc.pp.pca" in norm_src


def test_run_cnv_pipeline_still_uses_the_parameters_the_template_records():
    from xenium_viewer.utils import cnv_analysis

    src = inspect.getsource(cnv_analysis.run_cnv_pipeline)
    for fragment in ("lfc_clip=lfc_clip", "dendrogram=False",
                     "calculate_gene_values=True", "drop_unmapped_genes=True"):
        assert fragment in src, fragment


# ── the pandas-3 Arrow shim, which the template inlines by hand ──────────────

def test_the_arrow_shim_actually_converts_the_index():
    """Regression: it looked right and did nothing.

    ``_running_mean`` slices the gene list with a **2-D** index array. Under
    pandas 3 the index is an ArrowStringArray, which routes that to pyarrow's
    ``take()`` and raises ``ArrowInvalid: only handle 1-dimensional arrays``.

    The shim converted the index to object dtype and assigned it back — but
    AnnData *re-infers* string dtypes on assignment to ``.obs``/``.var``, so
    with ``future.infer_string`` at its default the Arrow array landed straight
    back where it started. The conversion has to happen with the option off
    across the assignment, and this executes the real shim to prove it does.
    """
    pd = pytest.importorskip("pandas")
    np = pytest.importorskip("numpy")
    anndata = pytest.importorskip("anndata")
    from xenium_viewer.tabs.tab_cnv import _CNV_ARROW_SHIM

    adata_cnv = anndata.AnnData(np.ones((4, 3), dtype="float32"))
    adata_cnv.var_names = ["GENE1", "GENE2", "GENE3"]
    adata_cnv.var["chromosome"] = ["chr1", "chr1", "chr2"]

    exec(compile(_CNV_ARROW_SHIM, "<shim>", "exec"),  # noqa: S102
         {"pd": pd, "adata_cnv": adata_cnv})

    genes = adata_cnv.var.index.values
    assert isinstance(genes, np.ndarray), type(genes).__name__
    # the operation that used to raise: 2-D fancy indexing
    assert genes[np.array([[0, 1], [1, 2]])].shape == (2, 2)


def test_the_shim_restores_the_option_it_changed():
    """It is a global pandas setting; leaving it off would change every
    DataFrame built later in the session."""
    pd = pytest.importorskip("pandas")
    np = pytest.importorskip("numpy")
    anndata = pytest.importorskip("anndata")
    from xenium_viewer.tabs.tab_cnv import _CNV_ARROW_SHIM

    adata_cnv = anndata.AnnData(np.ones((2, 2), dtype="float32"))
    before = pd.options.future.infer_string
    exec(compile(_CNV_ARROW_SHIM, "<shim>", "exec"),  # noqa: S102
         {"pd": pd, "adata_cnv": adata_cnv})

    assert pd.options.future.infer_string == before


def test_the_template_shim_matches_the_helper_run_cnv_pipeline_uses():
    """Two hand-maintained copies of one conversion; they drifted once already.

    ``run_cnv_pipeline`` calls ``_convert_adata_arrow_strings``, which disables
    ``future.infer_string`` around the assignment. The template inlines the same
    logic as plain source — and omitted exactly that line, so the CopyKAT path
    worked while inferCNV died.
    """
    from xenium_viewer.tabs.tab_cnv import _CNV_ARROW_SHIM
    from xenium_viewer.utils.adata_persistence import _convert_adata_arrow_strings

    helper = inspect.getsource(_convert_adata_arrow_strings)
    for fragment in ("pd.options.future.infer_string = False",
                     "is_string_dtype", "rename_categories",
                     "to_numpy(dtype=object)"):
        assert fragment.replace("pd.", "") in helper.replace("pd.", ""), fragment
        assert fragment in _CNV_ARROW_SHIM, f"template is missing: {fragment}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
