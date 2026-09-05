"""The ``environment`` node: what the analysis was run with.

A notebook only reproduces a result against the same software. The recorded code
named the functions but never the versions that answered the call, so a replay
that disagreed left no way to separate a real difference from a scanpy upgrade.

These tests pin the three properties that make the node useful rather than
decorative: it records the versions *and* executes, it sorts first without
becoming a dependency of everything, and re-opening a dataset in an unchanged
environment does not rewrite it.

Qt is imported by the tab helpers, so run with ``QT_QPA_PLATFORM=offscreen``.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("numpy")
pytest.importorskip("qtpy")

from palms.utils.environment import (  # noqa: E402
    RECORDED_PACKAGES, environment_code, package_versions, pins, same_environment,
)


@pytest.fixture
def ctx(tmp_path, qapp):
    from palms.tabs._helpers import create_shared_helpers
    from palms.utils.viewer_context import ViewerContext

    context = ViewerContext(
        data_path=tmp_path,
        state={"record_code": True, "code_journal": [],
               "prov_graph_restored": True},
    )
    create_shared_helpers(context)
    return context


# ── the cell itself ──────────────────────────────────────────────────────────

def test_the_recorded_versions_are_the_ones_actually_installed():
    versions = package_versions()
    assert versions["python"] == sys.version.split()[0]
    assert versions["numpy"] is not None            # a hard dependency

    code = environment_code(versions)
    recorded = {line.split()[1]: line.split()[2] for line in pins(code)}
    assert recorded["numpy"] == versions["numpy"]
    assert recorded["python"] == versions["python"]


def test_the_recorded_palms_version_is_the_sources_not_the_installs():
    """xv-0ob: an editable install's dist-info is frozen at install time.

    ``importlib.metadata.version('palms')`` answers with whatever the checkout
    was installed as, which on a dev box that predates a release is a version
    older than most of the code it is describing — and this value goes into the
    provenance node and into ``report_*.json``, which is the evidence behind a
    reproducibility claim. The module attribute tracks the source, so it is what
    the recorded pin must say. Third-party names are the opposite case and stay
    on the metadata, which this asserts too.
    """
    from importlib.metadata import version as installed

    import palms
    from palms.utils.environment import THIS_DISTRIBUTION

    versions = package_versions()
    assert versions[THIS_DISTRIBUTION] == palms.__version__
    assert versions["numpy"] == installed("numpy")


def test_the_cell_is_executable_python_that_seeds_the_rngs():
    """A comment block would tell a replay nothing it could act on."""
    import random

    import numpy as np

    code = environment_code({"python": "3.12.0", "numpy": "2.1.0"},
                            timestamp="2026-07-29 10:00:00 BST", seed=7)
    ast.parse(code)

    exec(compile(code, "<environment>", "exec"), {})  # noqa: S102
    after_cell = (random.random(), float(np.random.random()))

    random.seed(7)
    np.random.seed(7)
    assert (random.random(), float(np.random.random())) == after_cell


def test_a_package_that_is_not_installed_is_named_rather_than_omitted():
    """Absence is information: infercnvpy missing explains a missing CNV cell."""
    code = environment_code({"python": "3.12.0", "infercnvpy": None})
    assert "not installed: infercnvpy" in code
    assert "#   infercnvpy  None" not in code


def test_the_recorded_set_covers_what_a_step_can_depend_on():
    for name in ("scanpy", "anndata", "squidpy", "numpy", "pandas",
                 "scikit-learn", "leidenalg", "palms"):
        assert name in RECORDED_PACKAGES, name


def test_two_cells_recorded_in_the_same_environment_compare_equal():
    """Only the pins matter — the timestamp differs on every call."""
    versions = {"python": "3.12.0", "scanpy": "1.11.0"}
    first = environment_code(versions, timestamp="2026-07-29 10:00:00 BST")
    second = environment_code(versions, timestamp="2026-07-30 22:15:00 BST")

    assert first != second
    assert same_environment(first, second)
    assert pins(first) == ["#   python  3.12.0", "#   scanpy  1.11.0"]


def test_an_upgrade_makes_the_cells_differ():
    a = environment_code({"scanpy": "1.11.0"})
    b = environment_code({"scanpy": "1.12.0"})
    assert not same_environment(a, b)


# ── how it sits in the graph ─────────────────────────────────────────────────

def test_recording_the_preamble_records_the_environment_with_it(ctx):
    ctx.record_preamble()
    graph = ctx.state["prov_graph"]

    node = graph.get("environment")
    assert node is not None
    assert node.kind == "setup"
    assert "sc.logging.print_header()" in node.code


def test_the_environment_sorts_first(ctx):
    ctx.record_preamble()
    ctx.record_node("clustering:k", "sc.tl.leiden(adata)", deps=["preamble"])

    assert ctx.state["prov_graph"].topo_sort()[0] == "environment"


def test_nothing_depends_on_the_environment(ctx):
    """A version change is something to read, not a reason to grey out results.

    If ``preamble`` depended on it, re-recording in an upgraded environment
    would flag every descendant stale — every clustering, every DEG table —
    which is a warning about the whole session where none is warranted.
    """
    ctx.record_preamble()
    ctx.record_node("clustering:k", "sc.tl.leiden(adata)", deps=["preamble"])
    graph = ctx.state["prov_graph"]

    assert graph.get("preamble").deps == []
    assert all("environment" not in node.deps
               for node in (graph.get(nid) for nid in graph.topo_sort()))


def test_re_recording_in_the_same_environment_leaves_the_node_alone(ctx):
    """Tabs call ``record_preamble`` freely; each call must not churn the cell.

    Re-recording would rewrite the timestamp, re-emit the cell into the flat
    journal and rewrite the sidecar — on every launch, for no new information.
    """
    ctx.record_preamble()
    original = ctx.state["prov_graph"].get("environment").code
    journal_length = len(ctx.state["code_journal"])

    ctx.record_preamble()
    ctx.record_environment()

    assert ctx.state["prov_graph"].get("environment").code == original
    assert len(ctx.state["code_journal"]) == journal_length


def test_a_restored_session_keeps_its_original_stamp_when_versions_match(ctx):
    """The graph comes back from disk; the node must survive the next launch."""
    from palms.utils.prov_graph import ProvGraph

    ctx.record_preamble()
    items = ctx.state["prov_graph"].to_list()
    recorded = ctx.state["prov_graph"].get("environment").code

    # next launch: same environment, graph restored from the sidecar
    ctx.state["prov_graph"] = ProvGraph.from_list(items)
    ctx.state.pop("_environment_code", None)
    ctx.record_preamble()

    assert ctx.state["prov_graph"].get("environment").code == recorded


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
