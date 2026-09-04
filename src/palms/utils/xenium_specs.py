"""Restate ``experiment.xenium`` for a dataset that is a *crop* of another.

A Crop Dataset export used to copy the parent's ``experiment.xenium`` verbatim,
which made every quantity in it describe the wrong dataset — measured on
``demo_data/crop_6``: the file said ``num_cells`` 299769 for a table holding
76,577 cells, with ``region_area``, ``total_cell_area`` and
``transcripts_per_cell`` likewise the parent's. It is the one Xenium-format
file a crop writes, so it is the one a reader would trust; and because
``loader``'s cache-freshness hash is a sha256 of it, a crop and its parent
hashed identically.

Three rules decide what happens to a key, and the third is the one that keeps
this honest as 10x's format moves:

1. **Carried** (:data:`CARRIED_KEYS`) — facts about the *run*, which stay true
   of any subset of it: the instrument, the panel, the software versions, the
   uuids identifying the acquisition, ``pixel_size``. These identify the run the
   crop came from, not the crop.
2. **Recomputed** (:data:`RECOMPUTED_KEYS`) — restated from the crop itself by
   :func:`crop_table_stats` and the caller's polygon area. Only quantities whose
   definition was *verified* against real 10x output are here; see that
   function's docstring for the measurement.
3. **Everything else is demoted, not deleted.** An unrecognised key is more
   likely a measured quantity (which would then be a lie) than a provenance
   string, so the default is to drop it from the top level — but the run's own
   file is preserved verbatim under ``palms_crop.run_specs``, so nothing is
   lost, it is only moved somewhere it cannot be mistaken for a statement about
   the crop. This mirrors ``store_inventory``'s "unrecognised defaults to not
   deletable": the allow-list decides only the yes cases.

The transcript-derived quantities (``num_transcripts``,
``num_transcripts_high_quality``, ``fraction_transcripts_assigned``,
``transcripts_per_100um``, ``nuclear_transcripts_per_100um``,
``thickness_of_high_quality_decoded_transcripts``), the segmentation fractions
and ``non_zero_matrix_entries`` are all demoted rather than recomputed,
deliberately — see :func:`crop_table_stats` for what each measurement said.

A recomputed key is written even when the parent's XOA version never defined it
(``total_cell_area`` and ``region_area`` are absent from a 2.0.0 file), since it
is a true statement about the crop; ``palms_crop.recomputed`` is what tells a
reader which keys PALMS wrote rather than 10x.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Optional

#: The block this module adds, naming the parent and what was done to the file.
CROP_KEY = "palms_crop"

#: Facts about the run that stay true of a subset of it. Everything here is
#: copied verbatim. The uuids and ``region_name`` identify the *parent*
#: acquisition — which is the point: they say where the crop came from.
#: ``pixel_size`` and ``analysis_sw_version`` (or ``xenium_ranger['version']``)
#: are the only keys ``spatialdata_io.xenium()`` itself reads, so both must stay.
CARRIED_KEYS = frozenset({
    "major_version", "minor_version", "patch_version",
    "run_name", "run_start_time", "region_name", "preservation_method",
    "cassette_name", "cassette_uuid", "slide_id", "well_uuid",
    "panel_type", "panel_design_id", "panel_predesigned_id",
    "panel_designer_version", "panel_name", "panel_organism",
    "panel_tissue_type", "panel_num_targets_predesigned",
    "panel_num_targets_custom",
    "chemistry_version", "segmentation_stain", "z_step_size", "pixel_size",
    "instrument_sn", "instrument_sw_version",
    "analysis_sw_version", "analysis_uuid", "xenium_ranger",
    "experiment_uuid", "roi_uuid", "calibration_uuid",
})

#: Keys this module restates from the crop. ``region_area`` comes from the
#: caller (the drawn polygon); the rest from :func:`crop_table_stats`.
RECOMPUTED_KEYS = (
    "num_cells",
    "transcripts_per_cell",
    "total_cell_area",
    "region_area",
)


def crop_table_stats(obs) -> dict[str, Any]:
    """The spec quantities derivable from a cropped table's ``obs``.

    Each formula was checked against a real 10x bundle
    (``output-XETG00283__0037793__Region_4``, XOA 3.1.0, 310,003 cells) by
    recomputing it from ``cells.parquet`` and comparing with what 10x itself
    wrote in ``experiment.xenium``. All four reproduce it *exactly*:

    ========================  =================================================
    ``num_cells``             ``len(obs)`` — 310003
    ``transcripts_per_cell``  **median** of ``transcript_counts`` — 209
    ``total_cell_area``       ``cell_area.sum()`` — 22178312.051457267
    ========================  =================================================

    ``transcripts_per_cell`` is reported as a float, which a *median* is: 10x
    writes an int in both bundles checked, but only because a median over that
    many integer counts lands on one. A crop of an even number of cells can sit
    between two, and rounding there would be a second definition.

    A key whose source column is absent is omitted rather than guessed, so a
    table that has been through something unusual loses the field instead of
    gaining a wrong one.

    What is deliberately *not* here, because the measurement said the obvious
    formula is wrong:

    - ``num_transcripts`` is **not** the row count of ``transcripts.parquet``
      (132,030,595 against 132,310,899 rows on that dataset — a 280,304-row
      gap this code has no definition for), so every quantity derived from it
      is demoted too.
    - ``non_zero_matrix_entries`` counts the **all-feature** matrix (63,557,814
      over 9,687 features), while the viewer's table is read with
      ``gex_only=True`` and holds only the 5,101 Gene Expression rows. Counting
      the crop's ``X`` would answer a different question under 10x's name —
      the same trap as the QC control-rate denominator.
    - ``transcripts_per_100um`` is the one that argues for checking a formula
      against **more than one XOA version**. ``transcript_counts.sum() /
      cell_area.sum() * 100`` reproduces the 3.1.0 bundle to the last digit
      (400.4995231102975) and is *wrong* on a 2.0.0 one — 63.006 against the
      81.975 10x wrote, a ratio of 1.301 that none of the obvious variants
      (``total_counts``, ``nucleus_area``, per-cell mean or median) closes. So
      the quantity has no single definition this code can state, and it is
      demoted like the rest.
    """
    stats: dict[str, Any] = {"num_cells": len(obs)}

    counts = obs["transcript_counts"] if "transcript_counts" in obs else None
    area = obs["cell_area"] if "cell_area" in obs else None

    if counts is not None and len(counts):
        # Median, not mean: verified against 10x's own value.
        stats["transcripts_per_cell"] = float(counts.median())
    if area is not None and len(area):
        total_area = float(area.sum())
        stats["total_cell_area"] = total_area
    return stats


def crop_specs(
    parent_specs: Mapping[str, Any],
    *,
    stats: Mapping[str, Any],
    region_area_um2: Optional[float],
    source_path: Any,
    tool: str = "palms",
) -> dict[str, Any]:
    """Build the ``experiment.xenium`` contents for a crop of *parent_specs*.

    *stats* is :func:`crop_table_stats` over the exported table;
    *region_area_um2* is the area of the drawn crop polygon, which is this
    dataset's analogue of ``region_area`` — it is PALMS's restatement, not
    10x's quantity, which is why ``palms_crop.recomputed`` names it.

    Cropping a crop does not nest: the parent's own ``palms_crop`` block is
    stripped from the embedded copy and its source path is prepended to
    ``derived_from``, so the chain stays flat and readable however many times a
    dataset is cropped.
    """
    parent_specs = dict(parent_specs)
    parent_block = parent_specs.pop(CROP_KEY, None)

    # ``run_specs`` is always the *acquisition's* own file, never an
    # intermediate crop's: cropping a crop inherits it rather than embedding
    # the parent whole. Otherwise each generation would carry a copy of the
    # one before it, and the measured facts a reader actually wants from it —
    # the segmentation fractions, the transcript counts — exist only in the
    # original anyway. ``derived_from`` is what names the intermediates.
    run_specs = parent_specs
    if isinstance(parent_block, Mapping):
        run_specs = dict(parent_block.get("run_specs") or parent_specs)

    specs: dict[str, Any] = {
        k: v for k, v in parent_specs.items() if k in CARRIED_KEYS
    }
    recomputed = []
    for key in RECOMPUTED_KEYS:
        value = region_area_um2 if key == "region_area" else stats.get(key)
        if value is None:
            continue
        specs[key] = float(value) if key != "num_cells" else int(value)
        recomputed.append(key)

    chain = []
    if isinstance(parent_block, Mapping):
        chain = [str(p) for p in parent_block.get("derived_from", [])]
    chain.append(str(source_path))

    specs[CROP_KEY] = {
        "tool": tool,
        "derived_from": chain,
        "recomputed": recomputed,
        # Read off the run's file rather than RECOMPUTED_KEYS, so a key whose
        # source column was missing is reported as demoted (it is) while a
        # recomputed one is not — and so the list stays the same at every
        # generation, since it describes ``run_specs``.
        "demoted": sorted(
            k for k in run_specs
            if k not in CARRIED_KEYS and k not in recomputed
        ),
        "note": (
            "This dataset is a crop of another. The keys listed in 'recomputed' "
            "describe this crop; every other key describes the run it was cut "
            "from. Keys in 'demoted' measured that run's extent and have no "
            "restatement here — they are kept verbatim under 'run_specs', "
            "which is always the original acquisition's own experiment.xenium."
        ),
        "run_specs": run_specs,
    }
    return specs


def read_specs(path: Path) -> dict[str, Any]:
    """Parse an ``experiment.xenium``. Raises ``OSError``/``ValueError``."""
    with open(Path(path)) as f:
        return json.load(f)


def is_crop(specs: Mapping[str, Any]) -> bool:
    """Whether *specs* was written by :func:`crop_specs`."""
    return isinstance(specs.get(CROP_KEY), Mapping)
