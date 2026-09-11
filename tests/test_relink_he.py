"""``palms-relink-he`` puts back the H&E path and pixel size a restore dropped.

The bug it repairs: ``tab_he_registration._restore_session`` did not pass
``he_path`` / ``he_pixel_size_um`` through to the restore handler, and the H&E
session attrs are written "computed value wins, ``None`` included" (which is how
clearing an H&E clears it), so one restore erased the path from the store
permanently. Stores written before that fix carry ``he_path = None`` and cannot
repair themselves -- nothing left on disk says where the image came from.

These tests are filesystem-only: a fake store is a directory holding a
``sdata_cached.zarr`` with a ``viewer_session`` group, which is all
``_read_prev_attrs`` and ``safe_group_update`` need.
"""

from __future__ import annotations

import numpy as np
import pytest
import tifffile
import zarr

from palms.scripts.relink_he import (
    find_datasets, he_file_facts, index_he_files, main, plan_dataset,
)


def _make_store(root, name, **attrs):
    data_path = root / name
    session = zarr.open_group(str(data_path / "sdata_cached.zarr"),
                              mode="a").require_group("viewer_session")
    for key, value in attrs.items():
        session.attrs[key] = value
    return data_path


def _make_he(path, height, width, px_um=0.25):
    """A tiny OME-TIFF that declares a physical pixel size, like a real slide."""
    tifffile.imwrite(
        str(path), np.zeros((height, width, 3), np.uint8), photometric="rgb",
        metadata={"PhysicalSizeX": px_um, "PhysicalSizeXUnit": "µm",
                  "PhysicalSizeY": px_um, "PhysicalSizeYUnit": "µm"},
    )
    return path


@pytest.fixture
def he_dir(tmp_path):
    directory = tmp_path / "slides"
    directory.mkdir()
    _make_he(directory / "section.ome.tif", 40, 64)
    _make_he(directory / "other.ome.tif", 24, 24)
    return directory


def test_he_file_facts_reads_shape_and_pixel_size(he_dir):
    assert he_file_facts(he_dir / "section.ome.tif") == (40, 64, pytest.approx(0.25))


def test_a_dropped_path_is_relinked(tmp_path, he_dir):
    data_path = _make_store(tmp_path, "ds", he_filename="section.ome.tif",
                            he_path=None, he_shape_yx=[40, 64])

    result = plan_dataset(data_path, index_he_files(he_dir), relink_all=False)

    assert result.status == "relinked"
    assert result.he_path == str(he_dir / "section.ome.tif")
    assert result.px_um == pytest.approx(0.25)


def test_a_name_that_matches_a_different_image_is_refused(tmp_path, he_dir):
    """The guard that makes this safe at all.

    A filename match is not evidence that it is the same section, and every
    transform the dataset carries is expressed in the original's pixels -- so a
    wrong relink would silently re-point a registration at the wrong image.
    """
    data_path = _make_store(tmp_path, "ds", he_filename="section.ome.tif",
                            he_path=None, he_shape_yx=[999, 999])

    result = plan_dataset(data_path, index_he_files(he_dir), relink_all=False)

    assert result.status == "shape-mismatch"
    assert result.he_path is None


def test_a_dataset_that_is_already_linked_is_left_alone(tmp_path, he_dir):
    data_path = _make_store(tmp_path, "ds", he_filename="section.ome.tif",
                            he_path=str(he_dir / "section.ome.tif"),
                            he_pixel_size_um=0.25, he_shape_yx=[40, 64])

    assert plan_dataset(data_path, index_he_files(he_dir),
                        relink_all=False).status == "linked"


def test_a_dataset_with_no_he_is_not_a_failure(tmp_path, he_dir):
    data_path = _make_store(tmp_path, "ds", roi_count=0)
    assert plan_dataset(data_path, index_he_files(he_dir),
                        relink_all=False).status == "no-he"


def test_a_missing_candidate_is_reported_not_guessed(tmp_path, he_dir):
    data_path = _make_store(tmp_path, "ds", he_filename="absent.ome.tif",
                            he_path=None, he_shape_yx=[40, 64])
    assert plan_dataset(data_path, index_he_files(he_dir),
                        relink_all=False).status == "not-found"


def test_find_datasets_walks_a_root(tmp_path, he_dir):
    _make_store(tmp_path / "run", "a", he_filename="section.ome.tif")
    _make_store(tmp_path / "run", "b", he_filename="other.ome.tif")

    assert {p.name for p in find_datasets([tmp_path / "run"], recursive=True)} == {"a", "b"}
    assert find_datasets([tmp_path / "run" / "a"], recursive=False) == \
        [(tmp_path / "run" / "a").resolve()]


def test_dry_run_writes_nothing(tmp_path, he_dir, capsys):
    data_path = _make_store(tmp_path, "ds", he_filename="section.ome.tif",
                            he_path=None, he_shape_yx=[40, 64])

    assert main([str(data_path), "--he-dir", str(he_dir), "--dry-run"]) == 0

    from palms.utils.session import _read_prev_attrs
    assert _read_prev_attrs(data_path / "sdata_cached.zarr").get("he_path") is None
    assert "--dry-run" in capsys.readouterr().out


def test_the_repair_survives_a_reread(tmp_path, he_dir):
    """The end the tool exists for: a later launch picks `from_file`, not the cache."""
    from palms.utils.session import _read_prev_attrs

    data_path = _make_store(tmp_path, "ds", he_filename="section.ome.tif",
                            he_path=None, he_shape_yx=[40, 64], flip_v=True)

    assert main([str(data_path), "--he-dir", str(he_dir)]) == 0

    attrs = _read_prev_attrs(data_path / "sdata_cached.zarr")
    assert attrs["he_path"] == str(he_dir / "section.ome.tif")
    assert attrs["he_pixel_size_um"] == pytest.approx(0.25)
    # Every other attr survives the group swap.
    assert attrs["flip_v"] is True
    assert attrs["he_filename"] == "section.ome.tif"
