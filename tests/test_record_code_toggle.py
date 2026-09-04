"""Every ``run_step`` path must work with Preferences -> "Record reproducible code" off.

The failure this pins (xv-24o) was uniform across all 31 templates and nothing in
the suite turned the flag off, so it survived every green run::

    KeyError: node 'normalize' depends on unknown node 'preamble'
              (dependencies must be recorded first)

Recording was gated in exactly one layer and ``run_step`` was not in it:
``_record_node`` returns early with the toggle off, ``app.py`` seeds the preamble
only when it is on — but ``StepExecutor.run`` upserts unconditionally, against a
graph that was deliberately left empty. Two things made it worse than a plain
crash: the step's code had *already executed* when the upsert raised (so the user
saw a failure for a computation that succeeded, and its output binding was
discarded), and ``run_step`` has none of ``_record_node``'s degrade-and-report
handling, so the exception propagated out of the tab callback.

The second half of the fix is the one with no crash to point at: with the toggle
off the session's own graph must be left *untouched*. ``_record_node`` still
returns early, so letting the migrated steps write into it would persist a graph
at exit with ``normalize`` and ``clustering:*`` present but ``environment``,
``clustering:<key>`` and every terminal missing — an analysis nobody ran that way.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("numpy")
pytest.importorskip("pandas")
pytest.importorskip("anndata")
pytest.importorskip("qtpy")


def _ctx(tmp_path: Path, *, record: bool):
    from palms.tabs._helpers import create_shared_helpers
    from palms.utils.viewer_context import ViewerContext

    context = ViewerContext(
        data_path=tmp_path,
        state={
            "record_code": record, "code_journal": [],
            "prov_graph_restored": True,
        },
    )
    create_shared_helpers(context)
    return context


def _step(node_id="normalize", deps=("preamble",), out="marker", value=1):
    from palms.utils.steps import Step

    return Step(id=node_id, template=f"{out} = {value}\n", deps=list(deps),
                outputs=[out])


def _ids(graph):
    return sorted(d["id"] for d in graph.to_list())


# ── the crash ────────────────────────────────────────────────────────────────

def test_a_step_runs_with_recording_off(tmp_path, qapp):
    """The acceptance case: no KeyError, and the output comes back."""
    ctx = _ctx(tmp_path, record=False)

    assert ctx.run_step(_step()) == {"marker": 1}


def test_a_step_runs_with_recording_off_and_no_preamble_call(tmp_path, qapp):
    """``app.py`` skips ``record_preamble()`` entirely when recording is off.

    So ``run_step`` cannot assume anything seeded the root for it. Note the
    session graph is not absent, it is *empty* — ``create_shared_helpers``
    seeds one either way — which is why the failure was a missing dep rather
    than a missing graph.
    """
    ctx = _ctx(tmp_path, record=False)
    assert _ids(ctx.state["prov_graph"]) == []

    assert ctx.run_step(_step()) == {"marker": 1}


def test_a_chain_of_steps_runs_with_recording_off(tmp_path, qapp):
    """A dependent step still resolves its dep — the scratch graph is a real graph."""
    ctx = _ctx(tmp_path, record=False)
    ctx.run_step(_step())

    assert ctx.run_step(
        _step("clustering:x", deps=("normalize",), out="labels", value=2)
    ) == {"labels": 2}


# ── the session graph is left alone ──────────────────────────────────────────

def test_recording_off_does_not_touch_the_session_graph(tmp_path, qapp):
    ctx = _ctx(tmp_path, record=False)
    ctx.run_step(_step())

    assert _ids(ctx.state["prov_graph"]) == [], (
        "a step recorded into the session graph with recording off; "
        "save_session would persist it at exit"
    )
    assert ctx.state["code_journal"] == []


def test_recording_off_keeps_a_restored_graph_exactly_as_it_was(tmp_path, qapp):
    """The realistic case: a dataset with previous work, opened with recording off.

    Here the graph is *not* empty — ``app.py`` restores it unconditionally, only
    the preamble re-emit is gated — so nothing crashes and the damage is silent.
    """
    from palms.utils.prov_graph import ProvGraph, SETUP

    ctx = _ctx(tmp_path, record=False)
    restored = ProvGraph()
    restored.upsert("preamble", "data_path = 1\n", kind=SETUP)
    ctx.state["prov_graph"] = restored
    before = restored.to_list()

    ctx.run_step(_step())

    assert restored.to_list() == before


def test_the_scratch_graph_is_seeded_with_a_real_preamble(tmp_path, qapp):
    """Not a stub: the dep has to resolve, and it is the same text either way."""
    ctx = _ctx(tmp_path, record=False)
    ctx.run_step(_step())
    scratch = ctx.state["_unrecorded_prov_graph"]

    assert _ids(scratch) == ["normalize", "preamble"]
    assert "data_path = Path(" in scratch.get("preamble").code


# ── recording on is unchanged ────────────────────────────────────────────────

def test_recording_on_is_unchanged(tmp_path, qapp):
    ctx = _ctx(tmp_path, record=True)
    ctx.record_preamble()
    ctx.run_step(_step())

    assert _ids(ctx.state["prov_graph"]) == ["environment", "normalize", "preamble"]
    assert "_unrecorded_prov_graph" not in ctx.state
    assert ctx.state["code_journal"], "the flat journal is still written"


def test_run_step_seeds_its_own_root_when_recording_is_on(tmp_path, qapp):
    """``run_step`` establishes its own precondition rather than inheriting one.

    ``app.py`` seeds the preamble at launch, but that is a side effect a call
    site cannot see; before this, missing it raised *after* the code had run.
    """
    ctx = _ctx(tmp_path, record=True)

    assert ctx.run_step(_step()) == {"marker": 1}
    assert ctx.state["prov_graph"].get("preamble") is not None


def test_turning_recording_back_on_records_into_the_new_graph(tmp_path, qapp):
    """``_on_record_toggled`` replaces ``state["prov_graph"]`` wholesale.

    The executor is built once and held on ``ctx``, so it has to follow that
    swap — otherwise steps after a re-enable land in the discarded graph.
    """
    from palms.utils.prov_graph import ProvGraph

    ctx = _ctx(tmp_path, record=False)
    ctx.run_step(_step())

    ctx.state["record_code"] = True
    ctx.state["prov_graph"] = ProvGraph()      # what the toggle handler does
    ctx.record_preamble()
    ctx.run_step(_step("clustering:x", deps=("preamble",), out="labels", value=2))

    assert "clustering:x" in _ids(ctx.state["prov_graph"])
