"""Launch -> operate -> persist -> relaunch, through the real widgets.

The lifecycle these tests cover has never been exercised end to end. Its two
halves each had a reason to be missed: the persistence block used to be
straight-line code after ``napari.run()``, unreachable without an event loop,
and every tab was only ever driven with a stubbed context, where nothing is
written and nothing comes back.

Each test re-reads the store with plain ``zarr`` / ``spatialdata`` as well as
through a fresh context. A reload that agreed with itself but not with the bytes
on disk would pass an in-memory-only check.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

RESOLUTION = 0.3
KEY = f"leiden_igraph_r{RESOLUTION}"
NODE = f"clustering:{KEY}"


def _session_attrs(dataset: Path) -> dict:
    """The persisted session, read without the viewer."""
    import zarr
    group = zarr.open_group(
        str(dataset / "sdata_cached.zarr" / "viewer_session"),
        mode="r", use_consolidated=False,
    )
    return dict(group.attrs)


def _table_on_disk(dataset: Path):
    """The table as ``spatialdata`` reads it back, independent of any context."""
    import warnings
    from spatialdata import read_zarr
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read_zarr(str(dataset / "sdata_cached.zarr"))["table"]


# ── a recorded analysis step ─────────────────────────────────────────────────

def test_a_leiden_run_survives_a_close_and_reload(e2e_dataset, launch):
    app = launch(e2e_dataset)

    assert NODE not in app.graph, "pick a resolution the fixture has not used"
    assert KEY not in app.obs_columns()

    app.rig.tab("Cells", "Clustering")
    app.rig.mg("Resolution").value = RESOLUTION
    app.rig.click("Run Leiden Clustering")
    app.rig.wait_idle("Run Leiden Clustering", timeout=300)
    app.settle()

    # A clustering is more than one obs column: the recorded step writes the
    # bare key (so the notebook reproduces it) and save_clustering_to_adata
    # writes the prefixed twin the viewer colours by.
    assert {KEY, f"clustering_{KEY}"} <= app.obs_columns()
    assert NODE in app.graph
    assert app.graph.get(NODE).deps == ["normalize"]
    assert not app.graph.get(NODE).stale

    recorded = app.graph.get(NODE).code
    app.shutdown()
    app.close()

    # ── on disk, without the viewer ──────────────────────────────────────────
    table = _table_on_disk(e2e_dataset)
    assert f"clustering_{KEY}" in table.obs.columns, "the column never reached the store"
    ids = {node["id"] for node in _session_attrs(e2e_dataset)["prov_graph"]}
    assert NODE in ids, "the step was persisted without its code"

    sidecar = json.loads((e2e_dataset / "viewer_cache" / "prov_graph.json").read_text())
    assert NODE in {node["id"] for node in sidecar}

    # ── and through a fresh launch ───────────────────────────────────────────
    again = launch(e2e_dataset)
    assert {KEY, f"clustering_{KEY}"} <= again.obs_columns()
    assert NODE in again.graph
    assert again.graph.get(NODE).code == recorded, "the reloaded step is not the one that ran"
    assert KEY in again.ctx.clusterings, "restored, but not offered to the tabs"


def test_a_leiden_run_reaches_the_derived_outputs(e2e_dataset, launch):
    """analysis.py and the notebook are rebuilt, and name the step that ran.

    The notebook write at the end of ``shutdown_viewer`` is inside a bare
    ``except Exception: pass``, so nothing else in the suite would notice it
    failing — the file would simply be the one from the previous session, or
    absent.
    """
    import nbformat

    assert not (e2e_dataset / "analysis.py").exists(), "the fixture ships derived outputs"
    assert not (e2e_dataset / "analysis_notebook.ipynb").exists()

    app = launch(e2e_dataset)
    app.rig.tab("Cells", "Clustering")
    app.rig.mg("Resolution").value = RESOLUTION
    app.rig.click("Run Leiden Clustering")
    app.rig.wait_idle("Run Leiden Clustering", timeout=300)
    app.settle()
    n_nodes = len(app.graph)
    app.shutdown()

    script = e2e_dataset / "analysis.py"
    assert script.exists(), "analysis.py was not written"
    assert KEY in script.read_text()

    notebook = e2e_dataset / "analysis_notebook.ipynb"
    assert notebook.exists(), "analysis_notebook.ipynb was not written"
    nb = nbformat.read(notebook, as_version=4)
    assert len(nb.cells) >= n_nodes, (
        f"{len(nb.cells)} cells for {n_nodes} nodes — the notebook is behind the graph")
    assert any(KEY in cell.source for cell in nb.cells if cell.cell_type == "code")


# ── the QC filter: the one action that rebinds the cell set ──────────────────

def test_the_qc_filter_survives_a_close_and_reload(e2e_dataset, launch):
    app = launch(e2e_dataset)
    all_cells = app.ctx.adata.n_obs
    assert app.ctx.state.get("qc_filter") is None

    app.rig.tab("Tools", "QC")
    app.rig.mg("Filter cells").value = True
    app.rig.mg("Min transcripts per cell").value = 60
    app.rig.mg("Filter genes").value = True
    app.rig.mg("Min cells per gene").value = 3
    app.rig.click("Apply filter")
    app.rig.wait_idle("Apply filter", timeout=300)
    app.settle()

    stored = app.ctx.state["qc_filter"]
    assert stored == {"min_counts": 60, "min_cells": 3}
    assert app.ctx.adata.n_obs < all_cells, "the filter kept every cell — pick a harder cutoff"
    kept = app.ctx.adata.n_obs

    # The filter starts a second lineage; it does not revise the first.
    assert app.ctx.cell_root() == "qc_filter"
    assert app.ctx.cell_scoped_id("normalize") == "normalize:qc"
    assert app.graph.get("qc_filter").barrier is True
    assert app.graph.get("normalize").deps == ["preamble"], "the first lineage was re-rooted"

    # Applying a *different* cutoff than the one the fixture last used revises
    # `qc_filter`, so the results computed under the old filter are stale — that
    # is upsert being honest. What must not happen is the unfiltered lineage
    # going with them: those results were computed on every cell, and still were.
    stale = {node.id for node in app.graph.nodes() if node.stale}
    assert "normalize" not in stale
    assert "clustering:leiden_igraph_r1.0" not in stale
    assert "normalize:qc" in stale, "the filtered lineage should have been revised"

    # Written on every recorded step, not only at exit.
    sidecar = json.loads((e2e_dataset / "viewer_cache" / "prov_graph.json").read_text())
    assert "qc_filter" in {node["id"] for node in sidecar}

    app.shutdown()
    app.close()

    assert _session_attrs(e2e_dataset)["qc_filter"] == stored

    # The filter narrows the analysis, never the store. Persisting the subset
    # would drop cells the raw output does not contain, irreversibly.
    assert _table_on_disk(e2e_dataset).n_obs == all_cells

    again = launch(e2e_dataset)
    assert again.ctx.state["qc_filter"] == stored
    assert again.ctx.adata.n_obs == kept, "the filter was stored but not re-applied"
    assert again.ctx.full_adata.n_obs == all_cells
    assert again.ctx.cell_root() == "qc_filter"
    # ensure_qc_filter runs before any tab restore for exactly this reason: a
    # handler reaching ensure_normalized first would normalise the unfiltered
    # table and upsert `normalize` back onto `preamble`.
    assert again.graph.get("normalize").deps == ["preamble"]
    # Re-applying the *stored* filter is not a revision of it, so the relaunch
    # must not stale anything the apply did not. It would, if the re-emitted
    # `qc_filter` code differed by so much as a literal's formatting.
    assert {node.id for node in again.graph.nodes() if node.stale} == stale, \
        "a relaunch under a stored filter changed which results are stale"


def test_the_qc_filter_repoints_the_label_raster(e2e_dataset, launch):
    """Dropped cells render transparent rather than borrowing another's colour.

    ``label_to_obs`` is indexed by the raster's pixel value and holds an obs row
    *position*, so a map left over from the unfiltered table paints each cell
    with some other cell's value and raises nothing.
    """
    app = launch(e2e_dataset)
    before = len(app.ctx.label_to_obs)

    app.rig.tab("Tools", "QC")
    app.rig.mg("Filter cells").value = True
    app.rig.mg("Min transcripts per cell").value = 60
    app.rig.click("Apply filter")
    app.rig.wait_idle("Apply filter", timeout=300)
    app.settle()

    assert len(app.ctx.label_to_obs) == before, "the map must keep its raster indexing"
    assert (app.ctx.label_to_obs == -1).any(), "no cell was marked dropped"
    live = app.ctx.label_to_obs[app.ctx.label_to_obs >= 0]
    assert live.max() < app.ctx.adata.n_obs, "the map still points past the filtered table"


# ── relaunching an untouched dataset must change nothing ─────────────────────

def test_relaunching_an_untouched_dataset_changes_nothing(e2e_dataset, launch):
    """No node moves, and nothing goes stale.

    The gate on a whole class of bug: `app.py` re-emits `preamble` for the
    current `data_path` on every launch, so anything that makes the re-emitted
    code differ — a moved dataset, an unseeded `segmentation_source`, a QC
    filter applied in the wrong order — flags every descendant stale and marks
    the entire notebook ⚠ for nothing. Each of those has happened.
    """
    first = launch(e2e_dataset)
    before = [node.__dict__.copy() for node in first.graph.nodes()]
    assert not any(node["stale"] for node in before)
    first.shutdown()
    first.close()

    second = launch(e2e_dataset)
    after = [node.__dict__.copy() for node in second.graph.nodes()]

    assert [n["id"] for n in after] == [n["id"] for n in before]
    stale = [n["id"] for n in after if n["stale"]]
    assert not stale, f"a relaunch flagged {len(stale)} node(s) stale: {stale}"
    assert after == before


def test_the_source_fixture_is_never_written_to(app, e2e_source, e2e_dataset):
    """The staged copy is the one that gets used.

    Asserted here as well as enforced in the fixture's teardown, so the property
    has a test that names it: a launch writes `palms.log`, repairs the cache and
    rewrites `analysis.py` beside whatever it is pointed at.
    """
    assert app.ctx.data_path == e2e_dataset
    assert e2e_source not in e2e_dataset.parents
    assert (e2e_dataset / "palms.log").exists(), "the launch logged somewhere else"
    assert not (e2e_source / "palms.log").exists()


@pytest.mark.parametrize("name", ["analysis.py", "analysis_notebook.ipynb", "palms.log"])
def test_the_committed_fixture_carries_no_derived_output(e2e_source, name):
    """Otherwise the tests above would pass on a stale file from a past session."""
    if e2e_source.name != "crop_7":
        pytest.skip("only the committed fixture is stripped")
    assert not (e2e_source / name).exists()
