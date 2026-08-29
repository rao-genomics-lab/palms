"""The H&E / ARMS image paths must stay lazy.

Loading an H&E used to build RAM fast enough to kill the session. Measured on a
16384x12288 RGB TIFF, load + write to the cache peaked at **2.61 GB** and the
session restore at **1.71 GB**; the cause in both cases was a dense
materialisation of the full-resolution image, not anything napari did.

Three properties are pinned here, each against the defect it describes:

* ``parse_rgb_image_for_store`` never densifies (2.61 -> 0.91 GB), and produces
  a *byte-identical* element to the ``np.asarray`` version it replaced.
* the session restore hands napari dask, not numpy (1.71 -> 0.49 GB). The
  identical ARMS code was fixed in ``9cad210`` and the H&E copy beside it was
  missed for five months, which is exactly the kind of drift a source guard is
  for.
* ``raster_io.level_is_computed`` still recognises a chained pyramid, since
  ``app._warn_if_pyramid_is_not_stored`` is the only thing standing between a
  user and a silent OOM on an untiled file.

Run with:  pytest tests/test_he_image_memory.py
"""
from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pytest

da = pytest.importorskip("dask.array")
tifffile = pytest.importorskip("tifffile")

from palms.utils.raster_io import level_is_computed  # noqa: E402
from palms.utils.registration import (  # noqa: E402
    describe_pyramid, load_he_pyramid, parse_rgb_image_for_store,
)

SRC = Path(__file__).resolve().parent.parent / "src" / "palms"


def _write_flat_rgb(path, shape=(512, 384, 3), tile=(128, 128)):
    """A TIFF with no internal pyramid — the case that synthesises levels."""
    rng = np.random.default_rng(0)
    img = rng.integers(0, 255, shape, dtype=np.uint8)
    tifffile.imwrite(path, img, tile=tile, photometric="rgb")
    return img


# ── the reader stays lazy ────────────────────────────────────────────────────

def test_load_he_pyramid_returns_dask_and_reads_nothing(tmp_path):
    path = tmp_path / "flat.tif"
    expected = _write_flat_rgb(path)

    pyramid, tif = load_he_pyramid(path)
    try:
        assert len(pyramid) == 5, "four levels are synthesised above the base"
        assert all(isinstance(level, da.Array) for level in pyramid)
        assert np.array_equal(np.asarray(pyramid[0]), expected)
        # Halving each time, spatial dims only.
        assert [tuple(l.shape) for l in pyramid] == [
            (512, 384, 3), (256, 192, 3), (128, 96, 3), (64, 48, 3), (32, 24, 3)
        ]
    finally:
        tif.close()


def test_describe_pyramid_names_a_chain_and_an_untiled_base(tmp_path):
    """The diagnostic has to distinguish the two failure modes, not just log."""
    path = tmp_path / "flat.tif"
    _write_flat_rgb(path)
    pyramid, tif = load_he_pyramid(path)
    try:
        text = describe_pyramid(pyramid, "flat.tif")
        assert "CHAINED" in text, "synthesised levels are a chain and must say so"
        assert "512, 384, 3" in text.replace("(", "").replace(")", "")
    finally:
        tif.close()

    # A stored level is one task per chunk; a chained one carries the levels below.
    stored = da.from_array(np.zeros((64, 64)), chunks=(32, 32))
    assert not level_is_computed(stored)
    assert level_is_computed(da.coarsen(np.mean, stored, {0: 2, 1: 2}))


# ── the writer stays lazy, and writes the same bytes ─────────────────────────

def test_parse_rgb_image_for_store_does_not_materialise(tmp_path):
    path = tmp_path / "flat.tif"
    _write_flat_rgb(path)
    pyramid, tif = load_he_pyramid(path)
    try:
        parsed, shape_yx = parse_rgb_image_for_store(pyramid[0])
    finally:
        tif.close()

    assert shape_yx == (512, 384)
    scale0 = parsed["scale0"]["image"].data
    assert isinstance(scale0, da.Array), "the base must reach the store lazily"
    assert scale0.shape == (3, 512, 384), "written as (c, y, x)"
    assert scale0.dtype == np.uint8
    # scale_factors=[2, 2, 2, 2] -> five levels, unchanged from the eager version.
    assert len(parsed.children) == 5


def test_parse_rgb_image_for_store_matches_the_eager_version(tmp_path):
    """Byte-identical to `np.asarray(...).astype(np.uint8)`, or caches diverge."""
    from spatialdata.models import Image2DModel

    path = tmp_path / "flat.tif"
    _write_flat_rgb(path)
    pyramid, tif = load_he_pyramid(path)
    # The tif handle must outlive the lazy element — map_blocks hides the store
    # from spatialdata's introspection, it does not detach the file.
    lazy, _ = parse_rgb_image_for_store(pyramid[0])
    base = np.asarray(pyramid[0])

    eager = Image2DModel.parse(
        np.transpose(base, (2, 0, 1)).astype(np.uint8), dims=("c", "y", "x"),
        scale_factors=[2, 2, 2, 2], chunks=(3, 1024, 1024),
    )
    for name in eager.children:
        assert np.array_equal(
            np.asarray(lazy[name]["image"].data),
            np.asarray(eager[name]["image"].data),
        ), f"{name} differs from what the eager writer produced"
    tif.close()


def test_parse_rgb_image_for_store_accepts_a_cyx_base():
    """A restored H&E is already (c, y, x); it must not be transposed again."""
    base = da.from_array(np.zeros((3, 40, 20), dtype=np.uint8), chunks=(3, 20, 20))
    parsed, shape_yx = parse_rgb_image_for_store(base)
    assert shape_yx == (40, 20)
    assert parsed["scale0"]["image"].data.shape == (3, 40, 20)


# ── source guards: the eager idioms must not come back ───────────────────────

def _function_source(path: Path, name: str) -> str:
    """Source of the innermost function called *name*, closures included."""
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
            return ast.unparse(node)
    raise AssertionError(f"{name} not found in {path.name}")


@pytest.mark.parametrize("module, func", [
    ("tabs/tab_he_registration.py", "_save_he_to_sdata"),
    ("tabs/tab_arms.py", "_save_arms_he_to_sdata"),
])
def test_image_writers_never_densify(module, func):
    src = _function_source(SRC / module, func)
    assert "np.asarray" not in src, (
        f"{func} must not pull the full-resolution image into RAM; "
        "use registration.parse_rgb_image_for_store"
    )
    assert "Image2DModel" not in src, f"{func} should go through the shared helper"


@pytest.mark.parametrize("module, func", [
    ("tabs/tab_he_registration.py", "_load_he_from_sdata"),
    ("tabs/tab_arms.py", "_load_arms_from_sdata"),
])
def test_session_restore_hands_napari_dask(module, func):
    src = _function_source(SRC / module, func)
    assert ".compute()" not in src, (
        f"{func} must stay lazy — napari fetches only the tiles it draws"
    )
    assert "np.transpose" not in src, f"{func} must use da.transpose"


def test_pyramid_warning_covers_every_image_element():
    """Scoping it to morphology_focus is why an H&E could die quietly."""
    src = _function_source(SRC / "app.py", "_warn_if_pyramid_is_not_stored")
    assert "sdata.images.get('morphology_focus')" not in src.replace('"', "'")
    assert "for name in list(sdata.images" in src


# ── replacing an H&E that the cache already holds ────────────────────────────

@pytest.mark.parametrize("module, func, element", [
    ("tabs/tab_he_registration.py", "_save_he_to_sdata", "he_image"),
    ("tabs/tab_arms.py", "_save_arms_he_to_sdata", "arms_he_image"),
])
def test_image_writers_can_replace_a_stored_element(module, func, element):
    """Loading a second H&E over a restored one must actually persist.

    `safe_write_element`'s guard inspects `sdata`'s element dicts, so a restored
    `he_image` makes the rewrite look unsafe even after the old napari layer is
    gone. Without `replace_backed=True` the write is refused, the new image
    displays, and the next launch brings back the old one.
    """
    src = _function_source(SRC / module, func)
    assert "replace_backed=True" in src, (
        f"{func} must tell safe_write_element it has torn down the old layer"
    )
    assert element in src


def test_custom_segmentation_only_opts_in_where_the_layer_is_gone():
    """`_on_update_sdata` writes the on-screen layer's own arrays; it must not."""
    load_path = SRC / "tabs" / "tab_segmentation.py"
    tree = ast.parse(load_path.read_text())
    calls = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "save_custom_seg_to_sdata"
    ]
    assert len(calls) == 2, "two callers, and they differ in what is safe"
    opted_in = [
        any(kw.arg == "replace_backed" for kw in c.keywords) for c in calls
    ]
    assert sorted(opted_in) == [False, True], (
        "exactly one caller may opt in: the one loading a *different* "
        "segmentation after _apply_custom_segmentation removed the old layer"
    )
