"""The `viewer_session` a cropped export carries, and what it deliberately omits.

Round-tripped through a real zarr store and read back with the same
`load_session` the viewer uses, because the point of this group is that the
*viewer* can read it — an in-memory dict proves nothing about that.
"""

from __future__ import annotations

import numpy as np
import pytest

from palms.utils import crop_session


class _Ctx:
    def __init__(self, **kw):
        self.he_state = kw.get("he_state", {})
        self.arms_state = kw.get("arms_state", {})
        self.external_images_state = kw.get("external_images_state", [])
        self.patch_overlays_state = kw.get("patch_overlays_state", [])
        self.state = kw.get("state", {})
        self.data_path = kw.get("data_path", "/src/dataset")


_REG = [[-2.3535, 0.0, 34396.65], [0.0, -2.3535, 53689.95], [0.0, 0.0, 1.0]]


def _attrs(**kw):
    ctx = kw.pop("ctx", None) or _Ctx()
    return crop_session.build_session_attrs(
        ctx,
        carried_elements=kw.pop("carried", set()),
        cluster_labels=kw.pop("cluster_labels", {}),
        graph_items=kw.pop("graph_items", []),
        **kw,
    )


# ── what must NOT be there ───────────────────────────────────────────────────

def test_no_affine_is_written_into_the_export_session():
    """One writer for geometry: the element.

    A second copy drifts the moment someone re-registers inside the export —
    that writes the element immediately but the session only at exit.
    """
    ctx = _Ctx(he_state={"affine_3x3": np.array(_REG), "he_filename": "slide.tif"},
               arms_state={"affine_3x3": np.array(_REG)})

    attrs = _attrs(ctx=ctx)

    for key in ("affine_3x3", "coarse_affine", "arms_affine_3x3"):
        assert key not in attrs, f"{key} must come from the element, not the session"


def test_the_flips_are_false_because_the_element_transform_contains_them():
    ctx = _Ctx(he_state={"flip_v": True, "flip_h": True},
               arms_state={"flip_v": True, "flip_h": True})

    attrs = _attrs(ctx=ctx)

    assert attrs["flip_v"] is False and attrs["flip_h"] is False
    assert attrs["arms_flip_v"] is False and attrs["arms_flip_h"] is False


def test_the_arms_tile_source_paths_are_nulled():
    """Otherwise the export re-reads the whole slide's tiles into a crop.

    `tab_arms._on_arms_restored` falls back to these files whenever the sdata
    tiles are empty, and they are absolute paths that are still valid — so the
    crop would come up covered in source-frame tiles it does not contain.
    """
    ctx = _Ctx(arms_state={"geojson_path": "/data/tiles.geojson",
                           "csv_path": "/data/clusters.csv",
                           "he_path": "/data/slide.svs"})

    attrs = _attrs(ctx=ctx)

    assert attrs["arms_geojson_path"] is None
    assert attrs["arms_csv_path"] is None
    assert attrs["arms_he_path"] is None
    assert attrs["he_path"] is None


def test_a_ui_row_keeps_its_settings_but_loses_its_affine_copy():
    ctx = _Ctx(patch_overlays_state=[{
        "element_name": "patch_a", "palette_name": "tab20",
        "affine_matrix": _REG, "opacity": 0.8,
    }])

    row = _attrs(ctx=ctx, carried={"patch_a"})["patch_overlays_ui"][0]

    assert "affine_matrix" not in row, "the element carries the matrix now"
    assert row["palette_name"] == "tab20" and row["opacity"] == 0.8


def test_a_ui_row_for_an_overlay_that_did_not_travel_is_dropped():
    """A row pointing at an absent element is a broken entry in the tab."""
    ctx = _Ctx(patch_overlays_state=[{"element_name": "patch_dropped"}])

    assert _attrs(ctx=ctx, carried=set())["patch_overlays_ui"] == []


# ── what must be there ───────────────────────────────────────────────────────

def test_the_filename_survives_because_the_patch_link_is_by_layer_name():
    ctx = _Ctx(he_state={"he_filename": "slide.tif"})
    assert _attrs(ctx=ctx)["he_filename"] == "slide.tif"


def test_migration_flags_are_set_so_the_first_launch_does_not_re_migrate():
    attrs = _attrs()
    for key in ("migrated_to_adata", "migrated_landmarks_to_sdata",
                "migrated_rank_genes_to_adata", "migrated_deg_to_sdata"):
        assert attrs[key] is True


def test_cluster_labels_travel():
    labels = {"leiden_r1.0": {0: "Pr LumEp", 1: "CD8+ T"}}
    assert _attrs(cluster_labels=labels)["cluster_labels"] == labels


# ── the provenance graph ─────────────────────────────────────────────────────

def test_graph_paths_are_rewritten_to_the_export():
    """Otherwise the first launch re-emits `preamble` for its own path, `upsert`
    sees a changed node, and every descendant comes up flagged stale."""
    items = [{"id": "preamble",
              "code": 'data_path = Path(r"/src/dataset")\nsdata = xenium(data_path)'}]

    out = crop_session.rewrite_graph_paths(items, "/src/dataset", "/out/crop_3")

    assert "/out/crop_3" in out[0]["code"]
    assert "/src/dataset" not in out[0]["code"]


def test_rewriting_is_prefix_only():
    """A sibling directory did not move, and must not be rewritten."""
    items = [{"id": "n", "code": 'a = r"/src/dataset_backup/x"\nb = r"/src/dataset/y"'}]

    out = crop_session.rewrite_graph_paths(items, "/src/dataset", "/out/crop_3")

    assert "/src/dataset_backup/x" in out[0]["code"], (
        "/src/dataset must not match /src/dataset_backup"
    )
    assert "/out/crop_3/y" in out[0]["code"]


def test_the_crop_export_note_names_what_was_not_carried():
    from palms.utils.crop_export import crop_export_note

    node_id, code, label = crop_export_note(
        "crop_3", "/out", True, ["patch_HE_R2: nothing inside the crop"])

    assert node_id == "crop_export:crop_3"
    assert "not carried: patch_HE_R2" in code
    assert "crop_3" in label


def test_the_note_is_built_in_exactly_one_place():
    """It is written into the carried graph *and* the source's own graph.

    Two constructions of the same string in two files is the drift `run_step`
    exists to prevent, so the tab must call the shared function.
    """
    from pathlib import Path

    tab = (Path(__file__).resolve().parent.parent / "src" / "palms"
           / "tabs" / "tab_crop_dataset.py").read_text()

    assert "crop_export_note" in tab, (
        "tab_crop_dataset must build the provenance note through "
        "crop_export.crop_export_note rather than formatting its own copy"
    )


# ── the round trip that matters ──────────────────────────────────────────────

def test_the_written_session_is_readable_by_load_session(tmp_path):
    """`load_session` is what the viewer calls; anything else is a proxy for it."""
    import spatialdata as sd
    from spatialdata.models import Image2DModel
    from palms.utils.session import load_session

    staging = tmp_path / "export"
    staging.mkdir()
    sdata = sd.SpatialData(images={"morphology_focus": Image2DModel.parse(
        np.zeros((3, 32, 32), dtype=np.uint8), dims=("c", "y", "x"))})
    sdata.write(str(staging / "sdata_cached.zarr"))

    graph = [{"id": "preamble", "code": "x = 1", "deps": [], "kind": "setup",
              "label": "Setup", "params": {}, "stale": False, "seq": 1}]
    attrs = _attrs(ctx=_Ctx(he_state={"he_filename": "slide.tif"}),
                   cluster_labels={"leiden_r1.0": {0: "Pr LumEp"}},
                   graph_items=graph, he_shape_yx=(2254, 16371))

    crop_session.write_export_session(staging, attrs, graph)

    got = load_session(staging / "sdata_cached.zarr")
    assert got is not None, "the export must not open as a sessionless store"
    assert got["he_filename"] == "slide.tif"
    assert got["affine_3x3"] is None, "no affine copy — the element owns placement"
    # Keys survive as written — `load_cluster_labels_from_sdata` produces int
    # cluster ids, and nothing in this path stringifies them.
    assert got["cluster_labels"] == {"leiden_r1.0": {0: "Pr LumEp"}}

    sidecar = staging / "viewer_cache" / "prov_graph.json"
    assert sidecar.exists(), "the graph must also land where app.py prefers to read it"
    import json
    assert json.loads(sidecar.read_text())[0]["id"] == "preamble"
