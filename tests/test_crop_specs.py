"""A crop's ``experiment.xenium`` must describe the crop, not its parent.

The export used to copy the parent's file verbatim, so every quantity in it was
the parent's — and since ``loader``'s cache-freshness fingerprint is a sha256 of
that file, a crop and its parent hashed identically as well.

The tests here assert the *partition* rather than a list of remembered keys:
every top-level key of the result is carried, recomputed, or the ``palms_crop``
block, with no fourth possibility. That is what makes an unrecognised key from a
future XOA version safe by default — it is demoted into ``run_specs``, where
it cannot be read as a statement about the crop.

The formulas the recompute uses were verified against a real 10x bundle rather
than derived; :func:`palms.utils.xenium_specs.crop_table_stats` records that
measurement. What is pinned here is that the code implements those formulas.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from palms.utils import xenium_specs
from palms.utils.xenium_specs import (
    CARRIED_KEYS,
    CROP_KEY,
    RECOMPUTED_KEYS,
    crop_specs,
    crop_table_stats,
)


#: A real ``experiment.xenium`` (XOA 3.1.0), trimmed to one key of each kind.
PARENT = {
    # carried — the run, the panel, the instrument, the software
    "major_version": 5, "minor_version": 1, "patch_version": 0,
    "run_name": "COMBAT_PROSTATE_CRC_04DEC24",
    "run_start_time": "2024-12-04T16:01:48Z",
    "region_name": "Region_4",
    "preservation_method": "ffpe",
    "cassette_name": "Prostate", "slide_id": "0037793",
    "panel_name": "hMulti_100g", "panel_organism": "Human",
    "chemistry_version": "v2",
    "pixel_size": 0.2125,
    "instrument_sn": "XETG00283", "instrument_sw_version": "3.1.0.0",
    "analysis_sw_version": "xenium-3.1.0.4",
    "analysis_uuid": "2bb6a123-dc79-4aee-9824-c65e2881fa88",
    "experiment_uuid": "6c33acf8-cd22-4578-a6dc-987ff82c1394",
    "z_step_size": 3.0,
    "segmentation_stain": "Xenium Multi-Tissue Stain",
    # recomputed — every one of these describes the parent's extent
    "num_cells": 310003,
    "transcripts_per_cell": 209,
    "total_cell_area": 22178312.051457267,
    "region_area": 68050616.9528125,
    # demoted — no restatement exists for a crop
    "num_transcripts": 132030595,
    "num_transcripts_high_quality": 117566440,
    "fraction_transcripts_assigned": 0.7555220180180671,
    "transcripts_per_100um": 400.4995231102975,
    "nuclear_transcripts_per_100um": 452.4533489151673,
    "thickness_of_high_quality_decoded_transcripts": 5.213768371076243,
    "non_zero_matrix_entries": 63557814,
    "segmented_cell_stain_frac": 0.977126027812634,
    "imported_cell_frac": 0.0,
    "images": {"morphology_focus_filepath": "morphology_focus/"},
    "xenium_explorer_files": {"cells_zarr_filepath": "cells.zarr.zip"},
}

#: The parent quantities that must never survive at the top level: each one
#: measures the parent's extent and has no definition this code can restate.
EXTENT_KEYS = (
    "transcripts_per_100um",
    "num_transcripts", "num_transcripts_high_quality",
    "fraction_transcripts_assigned", "nuclear_transcripts_per_100um",
    "thickness_of_high_quality_decoded_transcripts",
    "non_zero_matrix_entries", "segmented_cell_stain_frac",
    "imported_cell_frac",
)


def _obs(n=4, counts=(10, 20, 30, 40), areas=(100.0, 100.0, 200.0, 200.0)):
    """A cropped table's obs, with the two columns the stats read."""
    return pd.DataFrame(
        {
            "transcript_counts": np.asarray(counts[:n], dtype=np.int32),
            "cell_area": np.asarray(areas[:n], dtype=np.float32),
            "cell_labels": np.arange(n, dtype=np.int64),
        }
    )


def _crop(obs=None, region_area=1234.5, parent=None, source="/data/parent"):
    return crop_specs(
        parent if parent is not None else PARENT,
        stats=crop_table_stats(obs if obs is not None else _obs()),
        region_area_um2=region_area,
        source_path=source,
    )


# ── the partition ───────────────────────────────────────────────────────────

def test_every_top_level_key_is_carried_recomputed_or_the_crop_block():
    """No fourth category: a key is provenance, a restatement, or demoted."""
    specs = _crop()
    unexplained = set(specs) - CARRIED_KEYS - set(RECOMPUTED_KEYS) - {CROP_KEY}
    assert unexplained == set()


def test_no_parent_extent_quantity_survives_at_the_top_level():
    specs = _crop()
    assert [k for k in EXTENT_KEYS if k in specs] == []


def test_a_demoted_key_is_kept_verbatim_under_run_specs():
    """Demotion is not deletion — nothing in the parent's file is lost."""
    specs = _crop()
    parent_copy = specs[CROP_KEY]["run_specs"]
    for key in EXTENT_KEYS:
        assert parent_copy[key] == PARENT[key]
    assert set(specs[CROP_KEY]["demoted"]) >= set(EXTENT_KEYS)


def test_an_unrecognised_key_is_demoted_rather_than_carried():
    """The default for a key from a future XOA version is 'not a fact about
    this crop', because the alternative failure is a silent lie."""
    parent = dict(PARENT, some_future_metric=17.0)
    specs = _crop(parent=parent)
    assert "some_future_metric" not in specs
    assert "some_future_metric" in specs[CROP_KEY]["demoted"]
    assert specs[CROP_KEY]["run_specs"]["some_future_metric"] == 17.0


def test_the_keys_spatialdata_io_reads_survive():
    """``pixel_size`` and the analyser version are the only two the reader
    itself touches; losing either makes the export unopenable."""
    specs = _crop()
    assert specs["pixel_size"] == PARENT["pixel_size"]
    assert specs["analysis_sw_version"] == PARENT["analysis_sw_version"]


def test_a_xenium_ranger_block_is_carried():
    """Resegmented data carries its version in ``xenium_ranger``, which the
    reader prefers over ``analysis_sw_version``."""
    parent = dict(PARENT, xenium_ranger={"version": "xenium-3.0.1.1"})
    assert _crop(parent=parent)["xenium_ranger"] == {"version": "xenium-3.0.1.1"}


# ── the restated quantities ─────────────────────────────────────────────────

def test_the_recomputed_quantities_describe_the_crop():
    obs = _obs()
    specs = _crop(obs=obs, region_area=1234.5)
    assert specs["num_cells"] == 4
    assert specs["transcripts_per_cell"] == 25.0            # median(10,20,30,40)
    assert specs["total_cell_area"] == pytest.approx(600.0)
    assert specs["region_area"] == 1234.5
    assert set(specs[CROP_KEY]["recomputed"]) == set(RECOMPUTED_KEYS)


def test_transcripts_per_100um_is_demoted_not_restated():
    """The obvious formula reproduces a 3.1.0 bundle exactly and is 30% out on a
    2.0.0 one, so there is no single definition to state — see
    ``crop_table_stats``. Pinned because the formula looks verified if you only
    ever check one dataset."""
    specs = _crop()
    assert "transcripts_per_100um" not in specs
    assert "transcripts_per_100um" in specs[CROP_KEY]["demoted"]
    assert "transcripts_per_100um" not in crop_table_stats(_obs())


def test_transcripts_per_cell_is_the_median_not_the_mean():
    """Pinned because 10x's own value is the median, verified exactly against a
    real bundle — a mean would be wrong by 37% on that dataset."""
    obs = _obs(counts=(1, 1, 1, 1000), areas=(1.0, 1.0, 1.0, 1.0))
    assert crop_table_stats(obs)["transcripts_per_cell"] == 1.0


def test_a_missing_source_column_omits_the_key_rather_than_guessing():
    obs = _obs().drop(columns=["cell_area"])
    specs = _crop(obs=obs)
    assert "total_cell_area" not in specs
    assert "num_cells" in specs
    # And the parent's value does not sneak back in through the carried set.
    assert "total_cell_area" in specs[CROP_KEY]["demoted"]


def test_an_absent_region_area_is_omitted_not_inherited():
    specs = _crop(region_area=None)
    assert "region_area" not in specs
    assert "region_area" in specs[CROP_KEY]["demoted"]


# ── the file on disk ────────────────────────────────────────────────────────

def test_the_result_is_json_serialisable_from_numpy_backed_obs():
    """obs columns hand back numpy scalars, which json.dumps refuses; the stats
    coerce them, so this is what catches a coercion going missing."""
    text = json.dumps(_crop())
    assert json.loads(text)["num_cells"] == 4


def test_a_crop_no_longer_hashes_the_same_as_its_parent(tmp_path):
    """loader's freshness fingerprint is a sha256 of this file, so a verbatim
    copy gave a crop the parent's identity."""
    parent_file = tmp_path / "parent.xenium"
    crop_file = tmp_path / "crop.xenium"
    parent_file.write_text(json.dumps(PARENT, indent=2))
    crop_file.write_text(json.dumps(_crop(), indent=2))
    assert (hashlib.sha256(parent_file.read_bytes()).hexdigest()
            != hashlib.sha256(crop_file.read_bytes()).hexdigest())


def test_read_specs_round_trips(tmp_path):
    path = tmp_path / "experiment.xenium"
    path.write_text(json.dumps(PARENT))
    assert xenium_specs.read_specs(path) == PARENT
    assert not xenium_specs.is_crop(PARENT)
    assert xenium_specs.is_crop(_crop())


# ── cropping a crop ─────────────────────────────────────────────────────────

def test_cropping_a_crop_does_not_nest_and_keeps_the_chain():
    first = _crop(source="/data/parent")
    second = crop_specs(
        first, stats=crop_table_stats(_obs(n=2)),
        region_area_um2=10.0, source_path="/data/crop_1",
    )
    assert second[CROP_KEY]["derived_from"] == ["/data/parent", "/data/crop_1"]
    # The embedded copy is the original run's file, not a crop block inside a
    # crop block — the intermediate is named in derived_from instead.
    assert CROP_KEY not in second[CROP_KEY]["run_specs"]
    assert second[CROP_KEY]["run_specs"]["num_transcripts"] == PARENT["num_transcripts"]
    assert second["num_cells"] == 2


def test_a_crop_of_a_crop_still_carries_the_run_provenance():
    first = _crop()
    second = crop_specs(first, stats=crop_table_stats(_obs()),
                        region_area_um2=1.0, source_path="/data/crop_1")
    assert second["run_name"] == PARENT["run_name"]
    assert second["pixel_size"] == PARENT["pixel_size"]


# ── source guard ────────────────────────────────────────────────────────────

def test_crop_export_never_copies_the_parent_experiment_file():
    """The defect was a ``shutil.copy`` of the parent's specs. Parsed rather
    than grepped so a rename of the constant cannot hide it."""
    src = Path(__file__).resolve().parents[1] / "src/palms/utils/crop_export.py"
    tree = ast.parse(src.read_text())
    copies = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr in {"copy", "copy2", "copyfile"}
        and any("experiment.xenium" in ast.dump(arg) for arg in node.args)
    ]
    assert copies == [], "experiment.xenium must be restated, not copied"
