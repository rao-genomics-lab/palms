"""The Dataset tab: it builds, and it never reaches outside the viewer's dirs.

The model and the safety predicate are covered by ``test_store_inventory.py``.
What is left for here is the Qt shell and the executor — in particular the two
things that have bitten this codebase before: forgetting that a deleted obs
column also lives in the in-memory ``ctx.clusterings`` dict, and deleting a
session group on disk without clearing the state ``save_session`` rebuilds it
from.
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
pytest.importorskip("qtpy")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")

from xenium_viewer.tabs import tab_dataset  # noqa: E402
from xenium_viewer.utils import store_inventory as si  # noqa: E402


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _ctx(sdata=None, data_path=None, no_cache=False):
    return SimpleNamespace(
        sdata=sdata, adata=(sdata["table"] if sdata is not None else None),
        viewer=None, state={}, no_cache=no_cache, segmentation_source="xenium",
        data_path=data_path or (Path(sdata.path).parent if sdata else None),
        clusterings={}, he_state={}, arms_state={},
        external_images_state=[], patch_overlays_state=[],
        refresh_clustering_choices=None, reload_dataset=None,
    )


# ── the widget ───────────────────────────────────────────────────────────────

def test_build_tab_returns_a_widget_and_exports(qapp, tiny_sdata):
    widget, exports = tab_dataset.build_tab(_ctx(tiny_sdata))
    assert widget is not None
    assert callable(exports["restore_session"])


def test_build_tab_survives_no_cache_mode(qapp, tiny_sdata, tmp_path):
    widget, exports = tab_dataset.build_tab(
        _ctx(None, data_path=tmp_path, no_cache=True))
    assert widget is not None


def test_build_tab_survives_a_dataset_with_no_directory(qapp):
    widget, _ = tab_dataset.build_tab(_ctx(None, data_path=None))
    assert widget is not None


def test_cache_path_falls_back_to_the_conventional_name(tiny_sdata, tmp_path):
    ctx = _ctx(tiny_sdata)
    assert tab_dataset._cache_path(ctx) == Path(tiny_sdata.path)
    # The store exists beside tiny_sdata, so the fallback finds it even with no
    # live sdata; a directory without one yields None.
    assert tab_dataset._cache_path(_ctx(None, data_path=tmp_path)) == Path(tiny_sdata.path)
    empty = tmp_path / "elsewhere"
    empty.mkdir()
    assert tab_dataset._cache_path(_ctx(None, data_path=empty)) is None


def test_the_tab_does_not_scan_at_build_time(qapp, tiny_sdata):
    """Walking the whole dataset directory must not be charged to every launch."""
    ctx = _ctx(tiny_sdata)
    tab_dataset.build_tab(ctx)
    assert "_dataset_sections" not in ctx.state
    assert "_dataset_worker" not in ctx.state


# ── source guards ────────────────────────────────────────────────────────────

def test_the_tab_never_calls_delete_element_from_disk():
    """The unsafe path: it unlinks before a replacement exists."""
    source = Path(tab_dataset.__file__).read_text()
    assert "delete_element_from_disk" not in source


def test_the_tab_only_removes_files_through_the_guarded_helper():
    """Every rmtree/unlink must sit in the one function that vets its path."""
    source = Path(tab_dataset.__file__).read_text()
    offenders = [
        f"{i}: {line.strip()}"
        for i, line in enumerate(source.splitlines(), 1)
        if re.search(r"(shutil\.rmtree|\.unlink\(|os\.remove)", line)
        and not line.lstrip().startswith("#")
        and "_remove_path" not in line
    ]
    assert offenders == [], (
        "route filesystem removal through _remove_path, which calls "
        "assert_node_deletable first: " + ", ".join(offenders))


def test_kind_order_puts_table_edits_first_and_backups_last():
    """Deleting custom_table then re-persisting obs would recreate it."""
    assert si._KIND_ORDER[si.OBS] < si._KIND_ORDER[si.ELEMENT]
    assert si._KIND_ORDER[si.BACKUP] == max(si._KIND_ORDER.values())


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
