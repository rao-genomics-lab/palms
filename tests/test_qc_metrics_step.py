"""The Xenium QC panel: the two places its arithmetic can quietly be wrong.

Both are about ``total_counts``, and neither moves the printed number much —
measured on ``Xenium_V1_human_Pancreas_FFPE_outs``, Xenium's own
``total_counts`` sums 5,512,036 against ``X``'s 5,511,215, a 0.015% difference
that leaves both control percentages identical to four decimals. So no
assertion on the *values* can catch a regression here, and these are source
guards instead:

- ``calculate_qc_metrics`` must run on a copy, because in the viewer ``adata``
  is the table inside the open store and the in-place form would persist a
  recomputed ``total_counts`` over one the viewer treats as structural;
- the control rates must read the *untouched* object, so the denominator stays
  every codeword class rather than gene expression alone.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")
sc = pytest.importorskip("scanpy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TEMPLATE_ID = "qc.metrics"


def _adata(n=60, g=20, areas=True, controls=True, seed=0):
    rng = np.random.default_rng(seed)
    X = sparse.csr_matrix(rng.poisson(0.8, (n, g)).astype("float32"))
    obs = pd.DataFrame(index=[f"c{i}" for i in range(n)])
    obs["cell_id"] = np.arange(1, n + 1)
    if controls:
        obs["control_probe_counts"] = rng.integers(0, 3, n)
        obs["control_codeword_counts"] = rng.integers(0, 2, n)
        # Xenium's total_counts is every codeword class, not X.sum(1).
        obs["total_counts"] = (np.asarray(X.sum(1)).ravel()
                               + obs["control_probe_counts"]
                               + obs["control_codeword_counts"])
    if areas:
        obs["cell_area"] = rng.uniform(20, 90, n)
        obs["nucleus_area"] = rng.uniform(5, 40, n)
    return anndata.AnnData(X=X, obs=obs,
                           var=pd.DataFrame(index=[f"Gene{i}" for i in range(g)]))


def _run(adata, blocks, tmp_path, percent_top=(10,)):
    from palms.utils.step_templates import builtin_assemble
    from palms.utils.steps import Step, StepExecutor

    paths = [str(tmp_path / "qc.png")]
    executor = StepExecutor(
        namespace={"sc": sc, "plt": plt, "Path": Path, "adata": adata})
    step = Step(id="plot:qc_metrics",
                template=builtin_assemble(TEMPLATE_ID, blocks),
                params={"title": "QC", "paths": paths, "percent_top": percent_top},
                outputs=["fig", "cprobes", "cwords"])
    return executor.run(step), paths


def _assemblies():
    from palms.utils.step_templates import builtin_spec
    return list(builtin_spec(TEMPLATE_ID).assemblies)


@pytest.mark.parametrize("blocks", _assemblies(), ids=lambda b: "+".join(b))
def test_every_assembly_draws_and_writes(blocks, tmp_path):
    out, paths = _run(_adata(), blocks, tmp_path)
    assert Path(paths[0]).exists()
    assert len(out["fig"].axes) == (4 if "areas" in blocks else 2)
    assert 0.0 <= out["cprobes"] <= 100.0
    assert 0.0 <= out["cwords"] <= 100.0
    plt.close(out["fig"])


@pytest.mark.parametrize("blocks", _assemblies(), ids=lambda b: "+".join(b))
def test_the_bound_table_is_never_annotated(blocks, tmp_path):
    """The whole reason for ``adata_qc = adata.copy()``."""
    adata = _adata()
    obs_before, var_before = list(adata.obs.columns), list(adata.var.columns)
    totals_before = adata.obs["total_counts"].to_numpy().copy()

    out, _ = _run(adata, blocks, tmp_path)
    plt.close(out["fig"])

    assert list(adata.obs.columns) == obs_before
    assert list(adata.var.columns) == var_before
    np.testing.assert_array_equal(adata.obs["total_counts"].to_numpy(), totals_before)


def test_control_rates_use_xeniums_total_not_the_recomputed_one(tmp_path):
    """The denominator is every codeword class, which is what the rate means."""
    adata = _adata()
    out, _ = _run(adata, ["head", "controls", "plot2", "save"], tmp_path)
    plt.close(out["fig"])

    expected = (adata.obs["control_probe_counts"].sum()
                / adata.obs["total_counts"].sum() * 100)
    assert out["cprobes"] == pytest.approx(expected)

    gene_only = (adata.obs["control_probe_counts"].sum()
                 / np.asarray(adata.X.sum()).ravel()[0] * 100)
    assert out["cprobes"] != pytest.approx(gene_only), (
        "the fixture must make the two denominators differ, or this test "
        "cannot tell them apart"
    )


def test_the_area_panels_need_the_columns_xenium_writes(tmp_path):
    """A custom-segmentation table need not carry them; the 2-panel form is why."""
    bare = _adata(areas=False)
    with pytest.raises(Exception):
        _run(bare, ["head", "controls", "plot4", "areas", "save"], tmp_path)

    out, _ = _run(bare, ["head", "controls", "plot2", "save"], tmp_path)
    assert len(out["fig"].axes) == 2
    plt.close(out["fig"])


def test_percent_top_is_clipped_to_the_panel():
    """scanpy raises "Positions outside range of features" past ``n_vars``."""
    from palms.tabs.tab_qc import PERCENT_TOP, percent_top_for

    assert percent_top_for(377) == PERCENT_TOP
    assert percent_top_for(60) == (10, 20, 50)
    assert percent_top_for(8) == ()


def test_a_small_panel_still_draws(tmp_path):
    """The clipping is not decorative: 8 genes is what the tab tests run on."""
    from palms.tabs.tab_qc import percent_top_for

    adata = _adata(g=8)
    out, _ = _run(adata, ["head", "controls", "plot2", "save"], tmp_path,
                  percent_top=percent_top_for(8))
    plt.close(out["fig"])


def test_the_shipped_text_keeps_the_two_objects_apart():
    """Source guard for both rules, since no value assertion can catch them."""
    from palms.tabs.tab_qc import _qc_metrics_template

    code = "\n".join(line for line in _qc_metrics_template(True).splitlines()
                     if not line.lstrip().startswith("#"))
    assert "adata_qc = adata.copy()" in code
    assert "sc.pp.calculate_qc_metrics(adata_qc," in code
    assert "sc.pp.calculate_qc_metrics(adata," not in code
    for line in code.splitlines():
        if "cprobes" in line or "cwords" in line or "control_" in line:
            assert "adata_qc.obs" not in line, (
                "the control rates must divide by Xenium's own total_counts, "
                f"not the recomputed one: {line!r}"
            )
