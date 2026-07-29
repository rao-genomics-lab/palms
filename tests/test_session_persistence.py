"""Session save must not destroy the session it is trying to write.

``save_session`` used to call ``create_group("viewer_session", overwrite=True)``
— which unlinks the group — and only write the replacement ~110 lines later.
Any exception in between (most plausibly a non-JSON-serialisable value reaching
``attrs.update``) left an empty group and one printed warning, at a point where
the user had already closed the window.

It also rebuilt attrs from scratch, so the four ``migrated_*`` markers were
wiped on every clean exit and their migrations re-ran at the next launch —
including two that themselves rewrote the whole cell table.

Most of these are pure-function tests over ``_build_session_attrs`` and need no
zarr at all.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_session_persistence.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")

from xenium_viewer.utils.session import (  # noqa: E402
    _PRESERVED_ON_SAVE, _build_session_attrs, _json_safe, _session_summary,
)


def _args(**overrides):
    base = {
        "state": {"segmentation_source": "xenium"},
        "he_state": {},
        "snapshot": {},
        "prev_attrs": {},
    }
    base.update(overrides)
    return base


# ── the migration-marker regression ──────────────────────────────────────────

def test_migration_markers_survive_a_save():
    """Wiping these re-armed the startup migrations on every single launch."""
    prev = {marker: True for marker in _PRESERVED_ON_SAVE}
    attrs = _build_session_attrs(**_args(prev_attrs=prev))
    for marker in _PRESERVED_ON_SAVE:
        assert attrs[marker] is True, marker


# ── the provenance graph must never shrink on save ───────────────────────────

def _graph(*ids):
    from xenium_viewer.utils.prov_graph import ProvGraph
    graph = ProvGraph()
    for node_id in ids:
        graph.upsert(node_id, f"# {node_id}", kind="setup")
    return graph


def test_a_saved_graph_is_never_replaced_by_a_smaller_one():
    """Regression, and it cost a real 13-node analysis.

    A launch that failed to restore came up holding only the preamble the tabs
    seed during construction. Its exit then wrote that single node over the
    stored graph — the last copy — because save simply serialised whatever was
    in memory. Nothing in the GUI removes nodes, so a shrink is always a bug.
    """
    stored = [{"id": f"n{i}", "code": "pass", "deps": []} for i in range(13)]
    attrs = _build_session_attrs(**_args(
        state={"prov_graph": _graph("preamble")},
        prev_attrs={"prov_graph": stored},
    ))
    assert len(attrs["prov_graph"]) == 13


def test_a_session_with_no_graph_does_not_wipe_the_stored_one():
    stored = [{"id": "preamble", "code": "pass", "deps": []}]
    attrs = _build_session_attrs(**_args(
        state={}, prev_attrs={"prov_graph": stored},
    ))
    assert attrs["prov_graph"] == stored


def test_a_grown_graph_is_saved():
    """The guard must not freeze the graph — growth is the normal case."""
    stored = [{"id": "preamble", "code": "pass", "deps": []}]
    attrs = _build_session_attrs(**_args(
        state={"prov_graph": _graph("preamble", "clustering:k")},
        prev_attrs={"prov_graph": stored},
    ))
    assert [item["id"] for item in attrs["prov_graph"]] == ["preamble", "clustering:k"]


def test_a_revised_graph_of_the_same_size_is_saved():
    """Same node count, different code — re-running a step must persist."""
    stored = [{"id": "preamble", "code": "OLD", "deps": []}]
    attrs = _build_session_attrs(**_args(
        state={"prov_graph": _graph("preamble")},
        prev_attrs={"prov_graph": stored},
    ))
    assert attrs["prov_graph"][0]["code"] == "# preamble"


def test_unknown_previous_keys_are_carried_forward():
    """Preserve-by-default, so a key added elsewhere isn't silently dropped."""
    attrs = _build_session_attrs(**_args(prev_attrs={"some_future_key": 42}))
    assert attrs["some_future_key"] == 42


def test_computed_values_win_over_previous_ones():
    """Clearing the H&E image must actually clear it, not fall back to prev."""
    attrs = _build_session_attrs(**_args(
        prev_attrs={"he_filename": "old.tif"}, he_state={},
    ))
    assert attrs["he_filename"] is None


def test_arms_attrs_fall_back_to_previous_when_absent():
    """tab_arms writes these in real time, so 'absent' means 'unchanged'."""
    attrs = _build_session_attrs(**_args(
        prev_attrs={"arms_he_filename": "arms.tif", "arms_csv_path": "/tmp/a.csv"},
    ))
    assert attrs["arms_he_filename"] == "arms.tif"
    assert attrs["arms_csv_path"] == "/tmp/a.csv"


# ── serialisability ──────────────────────────────────────────────────────────

def test_built_attrs_are_json_serializable():
    state = {
        "segmentation_source": "custom",
        "cluster_labels": {"leiden_r1.0": {0: "Tumour", 1: "Stroma"}},
        "marker_genes_json": '{"A": ["g1"]}',
        "rank_genes_groupby": "leiden_r1.0",
    }
    snapshot = {
        "roi_data": [np.zeros((4, 2))],
        "arms_state": {"affine_3x3": np.eye(3), "he_filename": "a.tif"},
        "external_images_ui": [{"name": "x"}],
    }
    attrs = _build_session_attrs(**_args(state=state, snapshot=snapshot))
    json.dumps(attrs)          # must not raise


def test_cluster_label_keys_are_stringified():
    state = {"cluster_labels": {"leiden": {0: "A", 1: "B"}}}
    attrs = _build_session_attrs(**_args(state=state))
    assert attrs["cluster_labels"] == {"leiden": {"0": "A", "1": "B"}}


def test_numpy_affine_becomes_a_plain_list():
    snapshot = {"arms_state": {"affine_3x3": np.eye(3)}}
    attrs = _build_session_attrs(**_args(snapshot=snapshot))
    assert attrs["arms_affine_3x3"] == [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]


def test_json_safe_drops_only_the_offender():
    """One bad value used to abort the whole save — after the group was gone."""
    safe, dropped = _json_safe({"good": 1, "bad": object(), "also_good": "x"})
    assert safe == {"good": 1, "also_good": "x"}
    assert dropped == ["bad"]


def test_a_non_serializable_prov_graph_does_not_abort_the_save():
    class Unserializable:
        def __len__(self):
            return 1

        def to_list(self):
            return [object()]

    attrs = _build_session_attrs(**_args(state={"prov_graph": Unserializable()}))
    safe, dropped = _json_safe(attrs)
    assert dropped == ["prov_graph"]
    assert safe["segmentation_source"] == "xenium"


def test_summary_reports_what_was_saved():
    attrs = {"roi_count": 3, "he_filename": "slide.tif", "has_rank_genes": True}
    summary = _session_summary(attrs)
    assert "3 ROIs" in summary and "slide.tif" in summary and "rank genes" in summary
    assert _session_summary({}) == "empty session"


def test_ui_residuals_survive_a_save_with_nothing_loaded():
    """An empty snapshot list means "none loaded now", not "forget them".

    It is empty before restore runs and right after a cache recovery, so
    overwriting on empty blanked recovered contrast/affine settings.
    """
    prev = {"external_images_ui": [{"element_name": "ext_a"}],
            "patch_overlays_ui": [{"element_name": "patch_b"}]}
    attrs = _build_session_attrs(**_args(prev_attrs=prev, snapshot={
        "external_images_ui": [], "patch_overlays_ui": [],
    }))
    assert attrs["external_images_ui"] == [{"element_name": "ext_a"}]
    assert attrs["patch_overlays_ui"] == [{"element_name": "patch_b"}]


def test_a_real_ui_snapshot_still_wins():
    prev = {"external_images_ui": [{"element_name": "old"}]}
    attrs = _build_session_attrs(**_args(prev_attrs=prev, snapshot={
        "external_images_ui": [{"element_name": "new"}],
    }))
    assert attrs["external_images_ui"] == [{"element_name": "new"}]


# ── on-disk behaviour ────────────────────────────────────────────────────────

def test_save_session_roundtrips(tiny_sdata):
    from xenium_viewer.utils.session import save_session

    cache = Path(tiny_sdata.path)
    save_session(cache, {"marker_genes_json": '{"A": ["g1"]}'},
                 {"he_filename": "slide.tif"}, {})

    import zarr
    session = zarr.open_group(str(cache / "viewer_session"), mode="r",
                              use_consolidated=False)
    assert session.attrs["he_filename"] == "slide.tif"
    assert "he" in session and "arms" in session


def test_a_failed_save_leaves_the_previous_session_intact(tiny_sdata, monkeypatch):
    """The reported failure shape: the group must survive a mid-save error."""
    import zarr
    from xenium_viewer.utils import session as session_mod

    cache = Path(tiny_sdata.path)
    session_mod.save_session(cache, {}, {"he_filename": "original.tif"}, {})

    monkeypatch.setattr(session_mod, "_write_array", lambda *a, **k: (
        _ for _ in ()).throw(RuntimeError("boom")))
    session_mod.save_session(cache, {}, {"he_filename": "clobbered.tif",
                                         "affine_3x3": np.eye(3)}, {})
    monkeypatch.undo()

    session = zarr.open_group(str(cache / "viewer_session"), mode="r",
                              use_consolidated=False)
    assert session.attrs["he_filename"] == "original.tif"


def test_markers_written_by_migrations_survive_a_real_save(tiny_sdata):
    """End-to-end version of the marker regression, through zarr."""
    import zarr
    from xenium_viewer.utils.session import save_session

    cache = Path(tiny_sdata.path)
    save_session(cache, {}, {}, {})
    store = zarr.open_group(str(cache / "viewer_session"), mode="r+",
                            use_consolidated=False)
    store.attrs["migrated_to_adata"] = True

    save_session(cache, {}, {}, {})

    reopened = zarr.open_group(str(cache / "viewer_session"), mode="r",
                               use_consolidated=False)
    assert reopened.attrs["migrated_to_adata"] is True


def test_clearing_a_registration_removes_its_stored_affine(tiny_sdata):
    """The seeded copy must not leave a stale affine behind."""
    import zarr
    from xenium_viewer.utils.session import save_session

    cache = Path(tiny_sdata.path)
    save_session(cache, {}, {"affine_3x3": np.eye(3)}, {})
    session = zarr.open_group(str(cache / "viewer_session"), mode="r",
                              use_consolidated=False)
    assert "affine_3x3" in session["he"]

    save_session(cache, {}, {}, {})
    session = zarr.open_group(str(cache / "viewer_session"), mode="r",
                              use_consolidated=False)
    assert "affine_3x3" not in session["he"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
