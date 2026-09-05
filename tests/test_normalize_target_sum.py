"""The normalisation target, and the memo that has to notice it changing.

``normalize`` had no parameters: it scaled every cell to 1e4, written as a
literal in the template. That is one of two conventions — ``scanpy``'s own
default is the median count across cells, which is what the Celldega tutorial
uses — and the difference is one of the reasons the same analysis run two ways
found different numbers of clusters on identical cells.

The setting is now a whole-block choice rather than a parameter that may be
``None``, for a reason about *reading* rather than about types: ``target_sum=
None`` in a recorded cell is valid scanpy and means the median, but it reads as
an argument someone failed to fill in.

No Qt, no dataset: a synthetic AnnData and the real templates.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
sc = pytest.importorskip("scanpy")
ad = pytest.importorskip("anndata")

from palms.tabs._helpers import (                              # noqa: E402
    DEFAULT_TARGET_SUM, normalize_preview,
)
from palms.utils.step_templates import builtin_spec            # noqa: E402
from palms.utils.steps import Step, StepExecutor               # noqa: E402


@pytest.fixture
def adata():
    rng = np.random.default_rng(0)
    a = ad.AnnData(rng.poisson(5, (40, 6)).astype("float32"))
    a.obs_names = [f"c{i}" for i in range(40)]
    a.var_names = [f"Gene{i}" for i in range(6)]
    return a


def _run(adata, target_sum):
    """Execute the normalize step the way ``ensure_normalized`` would."""
    from palms.utils.step_templates import builtin_assemble

    blocks, params, _ = normalize_preview(target_sum)
    step = Step(id="normalize", template=builtin_assemble("normalize", blocks),
                params=params, kind="setup", outputs=["adata_norm"])
    ex = StepExecutor(namespace={"sc": sc, "np": np, "pd": pd, "adata": adata})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = ex.run(step)
    return step, out["adata_norm"]


# ── the two conventions ──────────────────────────────────────────────────────

def test_a_fixed_target_scales_every_cell_to_it(adata):
    step, norm = _run(adata, 1e4)
    assert "target_sum=10000.0" in step.render()
    # log1p has been applied, so undo it before checking the row sums.
    sums = np.expm1(norm.X).sum(axis=1)
    assert np.allclose(sums, 1e4, rtol=1e-3)


def test_the_median_variant_passes_no_target_and_scanpy_uses_the_median(adata):
    step, norm = _run(adata, None)
    source = step.render()
    assert "sc.pp.normalize_total(adata_norm)" in source
    assert "target_sum=" not in source

    expected = float(np.median(adata.X.sum(axis=1)))
    sums = np.expm1(norm.X).sum(axis=1)
    assert np.allclose(sums, expected, rtol=1e-3)
    assert not np.allclose(expected, 1e4), (
        "the fixture must make the two conventions distinguishable, or this "
        "test would pass against either one"
    )


def test_both_conventions_are_declared_assemblies():
    spec = builtin_spec("normalize")
    assert ("copy", "scale.fixed", "tail") in spec.assemblies
    assert ("copy", "scale.median", "tail") in spec.assemblies


def test_the_default_is_the_viewers_historical_behaviour():
    """A dataset opened after the upgrade must analyse as it did before."""
    assert DEFAULT_TARGET_SUM == 1e4
    assert normalize_preview(DEFAULT_TARGET_SUM).blocks == [
        "copy", "scale.fixed", "tail"]


# ── the memo ─────────────────────────────────────────────────────────────────

def test_changing_the_target_re_runs_the_step_rather_than_reusing_the_copy(adata):
    """The defect this guards: ``_norm_src_id`` keyed on ``id(ctx.adata)`` alone.

    The cell set has not changed when the target changes, so a memo that only
    watches the cells hands back the ``adata_norm`` the *old* setting produced —
    and every later analysis silently uses it while the graph records the new
    code. Mirrors ``_ensure_spatial_neighbors``, which keys on its ``n_neighs``
    for the same reason.
    """
    import ast

    source = (Path(__file__).resolve().parent.parent / "src" / "palms"
              / "tabs" / "_helpers.py").read_text()
    tree = ast.parse(source)
    func = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_ensure_normalized")
    body = ast.unparse(func)
    assert "key = (id(ctx.adata), target_sum)" in body, (
        "the scaling target must be part of the memo key"
    )
    assert "state['_norm_src_id'] = key" in body


def test_the_two_settings_give_different_matrices(adata):
    """The end-to-end reason the memo matters."""
    _, fixed = _run(adata, 1e4)
    _, median = _run(adata, None)
    assert not np.allclose(fixed.X, median.X)


# ── the tab ──────────────────────────────────────────────────────────────────

def test_the_tab_writes_the_state_key_ensure_normalized_reads():
    """Two halves of one setting, in different modules — so pin the name."""
    import ast

    tab = (Path(__file__).resolve().parent.parent / "src" / "palms" / "tabs"
           / "tab_preprocess.py").read_text()
    helpers = (Path(__file__).resolve().parent.parent / "src" / "palms" / "tabs"
               / "_helpers.py").read_text()
    assert '"normalize_target_sum"' in tab
    tree = ast.parse(helpers)
    func = next(n for n in ast.walk(tree)
                if isinstance(n, ast.FunctionDef) and n.name == "_ensure_normalized")
    assert "normalize_target_sum" in ast.unparse(func)


def test_the_setting_round_trips_through_the_session():
    """A stored ``None`` means the median and must not read back as "unset"."""
    from palms.utils.session import _UNSET_TARGET_SUM, _build_session_attrs

    def attrs(**state):
        base = {"segmentation_source": "xenium"}
        base.update(state)
        return _build_session_attrs(state=base, he_state={}, snapshot={},
                                    prev_attrs={})

    assert attrs(normalize_target_sum=1e4)["normalize_target_sum"] == 1e4
    assert attrs(normalize_target_sum=None)["normalize_target_sum"] is None
    assert attrs()["normalize_target_sum"] == _UNSET_TARGET_SUM, (
        "a store that never held the setting must be distinguishable from one "
        "where the median was chosen — None cannot say both"
    )


def test_a_store_written_before_this_feature_opens_on_the_old_behaviour():
    """The opt-in guarantee, read at the source since app.py needs Qt."""
    import ast

    app = (Path(__file__).resolve().parent.parent / "src" / "palms"
           / "app.py").read_text()
    tree = ast.parse(app)
    assigns = [ast.unparse(n) for n in ast.walk(tree)
               if isinstance(n, (ast.Assign, ast.Expr))]
    seeded = [a for a in assigns if "normalize_target_sum" in a]
    assert any("DEFAULT_TARGET_SUM" in a and "_UNSET_TARGET_SUM" in a
               for a in seeded), (
        "the unset sentinel must fall back to the historical 1e4, not to None"
    )


def test_the_setting_is_listed_but_never_deletable():
    from palms.utils.store_inventory import _BLOCKED_SESSION_ATTRS

    assert "normalize_target_sum" in _BLOCKED_SESSION_ATTRS
