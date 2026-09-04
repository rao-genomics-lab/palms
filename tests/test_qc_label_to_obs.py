"""Re-pointing the label map is the one thing here that fails without raising.

``label_to_obs`` is indexed by the raster's *pixel value* and holds an obs *row
position*. ``CellColorManager`` then does

    color_arr[valid_labels] = rgba_obs[obs_indices]

so if the table loses rows and the map is not re-pointed, every cell after the
first dropped one is painted with some other cell's value — and nothing raises,
as long as the largest surviving row index is still inside the shorter table.

The property that matters is therefore not "the numbers changed" but "a label
still names the same cell", which is what this asserts.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")
anndata = pytest.importorskip("anndata")


def _table(n=50, seed=1):
    rng = np.random.default_rng(seed)
    return anndata.AnnData(
        np.zeros((n, 4), dtype="float32"),
        obs=pd.DataFrame({"cell_id": rng.permutation(np.arange(1, n + 1))},
                         index=[f"c{i}" for i in range(n)]),
    )


@pytest.mark.parametrize("seed", range(6))
def test_a_label_still_names_the_same_cell(seed):
    from palms.loader import label_to_obs_for
    from palms.utils.rebind_cells import kept_mask, repoint_label_to_obs

    full = _table(seed=seed)
    keep = np.random.default_rng(seed).random(full.n_obs) > 0.4
    subset = full[keep].copy()

    full_map = label_to_obs_for(full)
    new_map = repoint_label_to_obs(full_map, kept_mask(full, subset))

    survivors = set(subset.obs_names)
    checked = 0
    for label in range(len(full_map)):
        if full_map[label] < 0:
            assert new_map[label] == -1
            continue
        name = full.obs_names[full_map[label]]
        if name in survivors:
            assert subset.obs_names[new_map[label]] == name
            checked += 1
        else:
            assert new_map[label] == -1
    assert checked == subset.n_obs


def test_the_map_keeps_its_length():
    """Its length is the raster's maximum label, not the table's row count.

    Shortening it to the surviving maximum would make ``CellColorManager`` build
    a colormap the raster runs off the end of.
    """
    from palms.loader import label_to_obs_for
    from palms.utils.rebind_cells import kept_mask, repoint_label_to_obs

    full = _table()
    full_map = label_to_obs_for(full)
    # Drop the cell holding the highest label, so a naive rebuild would shrink.
    highest = int(np.argmax(full.obs["cell_id"].to_numpy()))
    keep = np.ones(full.n_obs, dtype=bool)
    keep[highest] = False
    subset = full[keep].copy()

    new_map = repoint_label_to_obs(full_map, kept_mask(full, subset))
    assert len(new_map) == len(full_map)
    assert len(label_to_obs_for(subset)) < len(full_map), (
        "the fixture must actually drop the top label, or this proves nothing"
    )


def test_dropped_cells_become_minus_one():
    from palms.loader import label_to_obs_for
    from palms.utils.rebind_cells import kept_mask, repoint_label_to_obs

    full = _table()
    keep = np.zeros(full.n_obs, dtype=bool)
    keep[:5] = True
    subset = full[keep].copy()
    new_map = repoint_label_to_obs(label_to_obs_for(full), kept_mask(full, subset))
    assert int((new_map >= 0).sum()) == 5


def test_a_mismatched_map_is_refused():
    """A map for a different table would silently gather the wrong rows."""
    from palms.utils.rebind_cells import repoint_label_to_obs

    with pytest.raises(ValueError):
        repoint_label_to_obs(np.array([-1, 0, 1, 2]), np.array([True, True]))


def test_the_colour_array_matches_the_filtered_rows():
    """End to end through the real manager: right length, right values."""
    pytest.importorskip("matplotlib")
    from palms.loader import label_to_obs_for
    from palms.utils.coloring import CellColorManager
    from palms.utils.rebind_cells import kept_mask, repoint_label_to_obs

    full = _table(n=30)
    full.obs["group"] = pd.Categorical(
        [str(i % 3) for i in range(full.n_obs)])
    keep = np.arange(full.n_obs) % 2 == 0
    subset = full[keep].copy()

    full_map = label_to_obs_for(full)
    new_map = repoint_label_to_obs(full_map, kept_mask(full, subset))
    manager = CellColorManager(subset, new_map)

    colors, _ = manager.get_cluster_colors(
        pd.Series(subset.obs["group"].astype(int).values,
                  index=subset.obs["cell_id"].values, name="group"))
    assert len(colors) == len(full_map), (
        "the colour array is indexed by label value, so it must cover the "
        "whole raster even when the table does not"
    )
    dropped = np.flatnonzero(new_map < 0)
    assert np.all(colors[dropped] == 0), "a filtered-out cell renders transparent"
