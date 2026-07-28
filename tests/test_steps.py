"""Unit tests for the step executor (pure, no Qt/napari/scanpy).

The point of these tests is the exactness guarantee: the string executed and
the string recorded must be the same object's value, for every parameter type
the GUI can produce.

Run standalone:   python tests/test_steps.py
Or with pytest:   pytest tests/test_steps.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from xenium_viewer.utils.prov_graph import SETUP, ProvGraph  # noqa: E402
from xenium_viewer.utils.steps import (  # noqa: E402
    ParamError, Step, StepError, StepExecutor,
    check_step, coerce, free_names, render, validate_param,
)


# ── the core guarantee ───────────────────────────────────────────────────────

def test_executed_source_is_recorded_source():
    """The string handed to exec is the string stored on the node."""
    ex = StepExecutor(namespace={"total": 0})
    step = Step(
        id="add",
        template="total = total + $amount",
        params={"amount": 7},
        outputs=["total"],
    )
    out = ex.run(step)

    assert out["total"] == 7
    assert ex.graph.get("add").code == step.render()
    assert ex.graph.get("add").code == "total = total + 7"


def test_params_reach_both_sides_identically():
    """A param appears in the recorded source exactly as it was executed."""
    ex = StepExecutor(namespace={})
    step = Step(
        id="clustering:leiden_r0.8",
        template=(
            "labels = leiden(resolution=$resolution, key_added=$key, "
            "flavor=$flavor, n_iterations=$n_iterations)"
        ),
        params={
            "resolution": 0.8,
            "key": "leiden_r0.8",
            "flavor": "igraph",
            "n_iterations": 2,
        },
        outputs=["labels"],
    )
    ex.ns["leiden"] = lambda **kwargs: kwargs
    out = ex.run(step)

    recorded = ex.graph.get("clustering:leiden_r0.8").code
    assert "resolution=0.8" in recorded
    assert "key_added='leiden_r0.8'" in recorded
    assert "n_iterations=2" in recorded
    # what ran carried the same values the source shows
    assert out["labels"] == {
        "resolution": 0.8, "key_added": "leiden_r0.8",
        "flavor": "igraph", "n_iterations": 2,
    }


def test_params_are_also_recorded_machine_readably():
    ex = StepExecutor()
    ex.run(Step(id="s", template="x = $n", params={"n": 3}, outputs=["x"]))
    assert ex.graph.get("s").params == {"n": 3}


# ── parameter validation ─────────────────────────────────────────────────────

@pytest.mark.parametrize("value", [
    1, 1.5, "text", True, False, None, [1, 2], ("a", "b"), {"k": 1}, {1, 2},
    "quote's and \"double\"", "C:\\path\\to\\file", 0.1, 1e-10,
])
def test_literal_params_round_trip(value):
    validate_param("p", value)  # must not raise
    assert eval(repr(value)) == value  # noqa: S307


def test_float_noise_still_round_trips_exactly():
    noisy = 1.0000000000000002
    validate_param("resolution", noisy)
    assert eval(repr(noisy)) == noisy  # noqa: S307


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_non_finite_floats_rejected(value):
    with pytest.raises(ParamError, match="non-finite"):
        validate_param("p", value)


def test_object_without_literal_form_rejected():
    class Opaque:
        pass

    with pytest.raises(ParamError, match="no literal form"):
        validate_param("p", Opaque())


def test_numpy_scalar_rejected_but_coercible():
    np = pytest.importorskip("numpy")
    with pytest.raises(ParamError):
        validate_param("resolution", np.float64(0.8))
    # the widget boundary is where this gets fixed
    coerced = coerce(np.float64(0.8))
    assert type(coerced) is float
    validate_param("resolution", coerced)


def test_coerce_handles_arrays_and_containers():
    np = pytest.importorskip("numpy")
    assert coerce(np.array([1, 2, 3])) == [1, 2, 3]
    assert coerce({"a": np.int64(2)}) == {"a": 2}
    assert coerce([np.float64(1.5), "x"]) == [1.5, "x"]
    assert type(coerce(np.bool_(True))) is bool


def test_bool_not_silently_widened_to_int():
    validate_param("flag", True)
    assert repr(True) == "True"


# ── template handling ────────────────────────────────────────────────────────

def test_undeclared_placeholder_fails_loudly():
    with pytest.raises(StepError, match=r"\$missing"):
        render("x = $missing", {})


def test_braces_in_template_survive():
    """string.Template leaves dict literals and f-strings alone."""
    out = render("d = {'k': $v}", {"v": 2})
    assert out == "d = {'k': 2}"


def test_free_names_finds_undeclared_dependency():
    assert free_names("y = adata.X") == {"adata"}
    assert free_names("import scanpy as sc\ny = sc.pp.pca(1)") == set()
    assert free_names("x = 1\ny = x + len('a')") == set()


def test_check_step_against_namespace():
    step = Step(id="s", template="result = sc.tl.leiden(adata, resolution=$r)",
                params={"r": 1.0})
    assert check_step(step, available={"sc", "adata"}) == set()
    assert check_step(step, available={"sc"}) == {"adata"}


# ── execution semantics ──────────────────────────────────────────────────────

def test_namespace_is_shared_across_steps():
    ex = StepExecutor()
    ex.run(Step(id="a", template="x = $v", params={"v": 2}, kind=SETUP, outputs=["x"]))
    ex.run(Step(id="b", template="y = x * $f", params={"f": 3},
                deps=["a"], outputs=["y"]))
    assert ex.ns["y"] == 6


def test_progress_reports_each_statement():
    ex = StepExecutor()
    seen = []
    step = Step(
        id="multi",
        template="a = $one\nb = a + $one\nc = b + $one",
        params={"one": 1},
        outputs=["c"],
    )
    ex.run(step, progress=lambda i, n, label: seen.append((i, n, label)))
    assert [(i, n) for i, n, _ in seen] == [(1, 3), (2, 3), (3, 3)]
    assert seen[0][2] == "a = 1"
    assert ex.ns["c"] == 3


def test_failure_names_the_step_and_records_nothing():
    ex = StepExecutor()
    step = Step(id="boom", template="x = $a\nraise ValueError('nope')",
                params={"a": 1})
    with pytest.raises(StepError, match=r"step 'boom' failed at statement 2/2"):
        ex.run(step)
    assert ex.graph.get("boom") is None  # no node claiming a nonexistent artifact


def test_syntax_error_names_the_step():
    ex = StepExecutor()
    with pytest.raises(StepError, match="not valid Python"):
        ex.run(Step(id="bad", template="x = = $v", params={"v": 1}))


def test_undeclared_output_is_an_error():
    ex = StepExecutor()
    with pytest.raises(StepError, match="did not bind"):
        ex.run(Step(id="s", template="x = $v", params={"v": 1}, outputs=["missing"]))


def test_rerun_revises_node_and_marks_descendants_stale():
    """The executor inherits prov_graph's upsert semantics for changed params."""
    ex = StepExecutor(graph=ProvGraph())
    ex.run(Step(id="clustering:leiden", template="labels = $r",
                params={"r": 1.0}, outputs=["labels"]))
    ex.run(Step(id="nhood:leiden", template="nh = labels",
                deps=["clustering:leiden"], outputs=["nh"]))
    assert ex.graph.get("nhood:leiden").stale is False

    ex.run(Step(id="clustering:leiden", template="labels = $r",
                params={"r": 0.5}, outputs=["labels"]))
    assert ex.graph.get("clustering:leiden").code == "labels = 0.5"
    assert ex.graph.get("nhood:leiden").stale is True


def test_missing_dependency_is_an_error_at_record_time():
    ex = StepExecutor()
    with pytest.raises(KeyError):
        ex.run(Step(id="s", template="x = $v", params={"v": 1}, deps=["nope"]))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
