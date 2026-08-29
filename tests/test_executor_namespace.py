"""The declared executor namespace and the real one must not drift.

``EXECUTOR_BASE_NAMES`` is what template validation checks a rendered template
against. If the executor is seeded with a name the declaration omits, a template
can come to rely on it, validation will still pass, and the notebook — whose
preamble binds only the declared set — replays with a ``NameError``. That is a
reproducibility failure that only shows up in a clean kernel, so it is worth a
cheap static guard here.

The call site is parsed rather than executed: building a real ``StepExecutor``
needs a ``ViewerContext`` with widgets, which means Qt, napari and scanpy.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from palms.utils.step_templates.namespace import (
    EXECUTOR_BASE_NAMES,
    NamespaceMismatch,
    check_base_namespace,
)

_HELPERS = Path(__file__).resolve().parent.parent / "src" / "palms" / "tabs" / "_helpers.py"


def _seeded_namespace_keys() -> set[str]:
    """The keys of the ``base_namespace`` dict literal in ``_get_executor``."""
    tree = ast.parse(_HELPERS.read_text())
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
        if "base_namespace" not in targets:
            continue
        assert isinstance(node.value, ast.Dict), (
            "base_namespace is no longer a dict literal; this test can no "
            "longer read it statically"
        )
        keys = set()
        for key in node.value.keys:
            assert isinstance(key, ast.Constant) and isinstance(key.value, str), (
                f"non-literal key in base_namespace: {ast.dump(key)}"
            )
            keys.add(key.value)
        return keys
    raise AssertionError("no `base_namespace = {...}` assignment found in _helpers.py")


def test_the_executor_is_seeded_with_exactly_the_declared_names():
    assert _seeded_namespace_keys() == set(EXECUTOR_BASE_NAMES)


def test_the_call_site_verifies_itself_at_runtime():
    """The static check above can be defeated; the runtime one cannot.

    A future refactor that builds the namespace some other way would slip past
    ``_seeded_namespace_keys``, so the call site must also call
    ``check_base_namespace`` on whatever it actually built.
    """
    source = _HELPERS.read_text()
    assert "check_base_namespace(" in source, (
        "_get_executor must validate its namespace against the declared set"
    )


def test_a_missing_name_is_rejected():
    with pytest.raises(NamespaceMismatch, match="missing.*adata"):
        check_base_namespace(EXECUTOR_BASE_NAMES - {"adata"})


def test_an_extra_name_is_rejected():
    """Both directions matter — see the module docstring."""
    with pytest.raises(NamespaceMismatch, match="unexpected.*adata_norm"):
        check_base_namespace(set(EXECUTOR_BASE_NAMES) | {"adata_norm"})


def test_the_exact_set_passes():
    check_base_namespace(EXECUTOR_BASE_NAMES)
    check_base_namespace({name: None for name in EXECUTOR_BASE_NAMES})


def test_step_produced_names_are_not_in_the_base_set():
    """``adata_norm`` and friends are bound by a step, not by the executor.

    A step that needs one must declare the producing step in ``deps``; putting
    it in the base set would let a template read it with no edge in the graph,
    which is exactly the lie the provenance DAG exists to prevent.
    """
    produced_by_steps = {
        "adata_norm", "adata_leiden", "adata_cnv", "rank_df", "rank_results",
        "roi_polygons", "roi_deg_df", "ligrec_res", "fig",
    }
    assert not (EXECUTOR_BASE_NAMES & produced_by_steps)
