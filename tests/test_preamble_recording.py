"""The recorded preamble must load the data the way *this* dataset can be read.

A Crop Dataset export has no raw 10x output — the zarr store written by the crop
*is* the data. ``spatialdata_io.xenium()`` cannot read one, so a notebook whose
preamble calls it fails on its very first cell, before any analysis runs.

That was true of every notebook recorded on a crop export, including the demo
dataset shipped in the release bundle: replaying it stopped at ``preamble`` with
``FileNotFoundError: .../cells.zarr.zip``. Nothing caught it, because recording
and persisting the graph both succeed — only *executing* it fails, and the only
thing that executes it is ``scripts/verify_notebook.py``.

These tests pin the branch rather than the wording, except where the wording is
the thing that has to be right: the call itself.
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


def _ctx(data_path: Path):
    from palms.tabs._helpers import create_shared_helpers
    from palms.utils.viewer_context import ViewerContext

    context = ViewerContext(
        data_path=data_path,
        state={
            "record_code": True, "code_journal": [],
            "prov_graph_restored": True,
        },
    )
    create_shared_helpers(context)
    return context


def _preamble_code(ctx) -> str:
    ctx.record_preamble()
    return ctx.state["prov_graph"].get("preamble").code


def test_raw_output_is_read_with_spatialdata_io(tmp_path, qapp):
    """With raw 10x output present, the preamble is unchanged: ``xenium()``."""
    # One marker is enough — has_raw_xenium_source is deliberately generous, so
    # that a *partially* broken raw dataset keeps raising its own error rather
    # than being reclassified as a crop export.
    (tmp_path / "cell_feature_matrix.h5").write_bytes(b"")

    code = _preamble_code(_ctx(tmp_path))

    assert "from spatialdata_io import xenium" in code
    assert "sdata = xenium(data_path)" in code
    assert "read_zarr" not in code


def test_a_crop_export_is_read_from_its_own_zarr(tmp_path, qapp):
    """With no raw output, the preamble reads the store the crop wrote."""
    code = _preamble_code(_ctx(tmp_path))

    assert 'sdata = sd.read_zarr(data_path / "sdata_cached.zarr")' in code
    assert "xenium(data_path)" not in code
    # ``sd`` has to be bound for that line to run.
    assert "import spatialdata as sd" in code


def test_a_declared_cache_only_export_is_read_from_zarr_even_beside_raw_files(
    tmp_path, qapp,
):
    """The acceptance case for xv-iy9.

    A raw-format crop export would write raw-shaped files *and* stamp
    ``cache_only`` in its manifest. Branching on file presence would then record
    ``xenium(data_path)``, which reads the raw half and silently drops every
    derived layer the crop carried — a worse failure than the loud one this
    branch was added to fix, because the notebook runs and produces less.
    """
    import json

    from palms.utils.zarr_safe import MANIFEST_FILE

    (tmp_path / "cell_feature_matrix.h5").write_bytes(b"")
    cache = tmp_path / "sdata_cached.zarr"
    cache.mkdir()
    (cache / MANIFEST_FILE).write_text(json.dumps({"cache_only": True}))

    code = _preamble_code(_ctx(tmp_path))

    assert 'sdata = sd.read_zarr(data_path / "sdata_cached.zarr")' in code
    assert "xenium(data_path)" not in code


def test_the_crop_export_preamble_still_derives_the_store_from_data_path(tmp_path, qapp):
    """``palms-rename-dataset`` rewrites exactly one line; keep it sufficient.

    The tool substitutes the recorded ``data_path = Path(r"…")``. If the store
    were recorded as its own absolute path, a renamed dataset would keep a
    notebook pointing at the old location — silently, since the first line would
    still be rewritten and look repaired.
    """
    code = _preamble_code(_ctx(tmp_path))

    assert f'data_path = Path(r"{tmp_path}")' in code
    assert str(tmp_path / "sdata_cached.zarr") not in code


def test_both_branches_bind_adata_the_same_way(tmp_path, qapp):
    """Everything downstream reads ``adata``; the load path must not change it."""
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "cell_feature_matrix.h5").write_bytes(b"")
    crop = tmp_path / "crop"
    crop.mkdir()

    for path in (raw, crop):
        assert 'adata = sdata["table"].copy()' in _preamble_code(_ctx(path))
