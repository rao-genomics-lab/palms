"""The real persistence functions must not leave a broken store behind.

``tests/test_zarr_safe.py`` covers the primitive. These exercise the call sites
that used delete-then-write, against a real on-disk store, so a future edit that
reintroduces the pattern at one of them fails here.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_persistence_safety.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

pytest.importorskip("spatialdata")
np = pytest.importorskip("numpy")
pytest.importorskip("geopandas")
pytest.importorskip("shapely")

from palms.utils import adata_persistence as ap  # noqa: E402
from palms.utils.zarr_safe import JOURNAL_DIR, STAGING_DIR  # noqa: E402


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _ctx(sdata):
    """The slice of ViewerContext the persistence functions actually touch."""
    return SimpleNamespace(
        sdata=sdata, adata=sdata["table"], no_cache=False,
        segmentation_source="xenium", annotation_layer=None,
    )


def _reread(cache):
    from spatialdata import read_zarr
    return read_zarr(str(cache))


def _clean(cache):
    """No journal and no staging left behind by a successful write."""
    return (
        list((cache / JOURNAL_DIR).glob("*.json")) == []
        and list((cache / STAGING_DIR).iterdir()) == []
    )


# ── the hot path ─────────────────────────────────────────────────────────────

def test_persist_table_roundtrips_and_leaves_no_debris(tiny_sdata):
    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    tiny_sdata["table"].obs["clustering_leiden_r1.0"] = ["a"] * tiny_sdata["table"].n_obs

    ap._persist_table(ctx)

    assert "clustering_leiden_r1.0" in _reread(cache)["table"].obs
    assert _clean(cache)


def test_persist_table_keeps_the_previous_version_recoverable(tiny_sdata):
    from palms.utils.zarr_safe import list_trash

    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    tiny_sdata["table"].obs["clustering_x"] = ["a"] * tiny_sdata["table"].n_obs
    ap._persist_table(ctx)

    assert list_trash(cache).get("tables/table")


def test_a_failed_persist_leaves_the_table_readable(tiny_sdata, monkeypatch):
    """The reported failure, at the real call site.

    Under delete-then-write this left no table on disk and the loader then
    discarded the whole cache.
    """
    from spatialdata import SpatialData

    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    monkeypatch.setattr(SpatialData, "write", lambda *a, **k: (_ for _ in ()).throw(
        OSError("No space left on device")))

    ap._persist_table(ctx)          # swallows, as it always has

    assert "table" in _reread(cache).tables
    assert list((cache / STAGING_DIR).iterdir()) == []


# ── shapes: the delete-then-early-return class ───────────────────────────────

def test_saving_rois_roundtrips(tiny_sdata):
    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    rois = [np.array([[0.0, 0.0], [0.0, 4.0], [4.0, 4.0], [4.0, 0.0]])]

    ap.save_rois_to_sdata(ctx, rois)
    assert "rois" in _reread(cache).shapes
    assert len(ap.load_rois_from_sdata(_reread(cache))) == 1
    assert _clean(cache)


def test_an_empty_roi_list_clears_rois_deliberately(tiny_sdata):
    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    ap.save_rois_to_sdata(ctx, [np.array([[0.0, 0.0], [0.0, 4.0], [4.0, 4.0]])])
    ap.save_rois_to_sdata(ctx, [])
    assert "rois" not in _reread(cache).shapes


def test_a_failed_roi_save_does_not_erase_the_stored_rois(tiny_sdata, monkeypatch):
    """Regression: the old order deleted the element and *then* built the new
    one, so a failure — or a transiently empty layer — lost the ROIs."""
    from spatialdata import SpatialData

    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    ap.save_rois_to_sdata(ctx, [np.array([[0.0, 0.0], [0.0, 4.0], [4.0, 4.0]])])
    assert "rois" in _reread(cache).shapes

    monkeypatch.setattr(SpatialData, "write", lambda *a, **k: (_ for _ in ()).throw(
        OSError("disk full")))
    ap.save_rois_to_sdata(ctx, [np.array([[1.0, 1.0], [1.0, 5.0], [5.0, 5.0]])])
    monkeypatch.undo()

    assert "rois" in _reread(cache).shapes


def test_saving_landmarks_roundtrips(tiny_sdata):
    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    pts = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])

    ap.save_landmarks_to_sdata(ctx, "he_xenium_landmarks", pts)
    loaded = ap.load_landmarks_from_sdata(_reread(cache), "he_xenium_landmarks")
    np.testing.assert_allclose(loaded, pts)
    assert _clean(cache)


def test_clearing_landmarks_removes_the_element(tiny_sdata):
    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    ap.save_landmarks_to_sdata(ctx, "he_xenium_landmarks", np.array([[1.0, 2.0]]))
    ap.save_landmarks_to_sdata(ctx, "he_xenium_landmarks", None)
    assert "he_xenium_landmarks" not in _reread(cache).shapes


def test_saving_arms_tiles_roundtrips(tiny_sdata):
    cache = Path(tiny_sdata.path)
    ctx = _ctx(tiny_sdata)
    polys = [np.array([[0.0, 0.0], [0.0, 2.0], [2.0, 2.0], [2.0, 0.0]])]

    ap.save_arms_tiles_to_sdata(ctx, polys, ["tile_1"], [3])
    stored = _reread(cache)["arms_tiles"]
    assert list(stored["tile_name"]) == ["tile_1"]
    assert list(stored["cluster_id"]) == [3]
    assert _clean(cache)


# ── the guard that keeps this from regressing ────────────────────────────────

def test_no_source_file_still_does_delete_then_write():
    """delete_element_from_disk immediately followed by write_element is the
    bug. Nothing in src/ may reintroduce it."""
    import re

    src = Path(__file__).resolve().parent.parent / "src" / "palms"
    pattern = re.compile(r"delete_element_from_disk\([^)]*\)", re.S)
    offenders = []
    for path in sorted(src.rglob("*.py")):
        if path.name == "zarr_safe.py":
            continue
        text = path.read_text()
        for match in pattern.finditer(text):
            line = text[: match.start()].count("\n") + 1
            offenders.append(f"{path.relative_to(src)}:{line}")
    assert offenders == [], (
        "delete_element_from_disk must not be called outside zarr_safe.py — "
        "use safe_write_element / safe_delete_element: " + ", ".join(offenders)
    )


def test_restore_is_not_gated_on_a_stored_session():
    """`restore_fn` must be reachable for a store that has no `viewer_session`.

    The defect: `restore_fn(session)` sat inside `if session is not None:`, under
    `if not no_cache and zarr_path.exists():`, with an `elif` that could therefore
    never run while the zarr existed. A cropped export — elements on disk, no
    session — opened with no H&E, no ROIs and no cluster names, and every
    `load_*_from_sdata` result computed just above was silently discarded.

    Parse the call rather than grep it: the property is *reachability*, and a
    regression would look like ordinary code.
    """
    import ast

    src = Path(__file__).resolve().parent.parent / "src" / "palms" / "app.py"
    tree = ast.parse(src.read_text(), str(src))

    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Name) and n.func.id == "restore_fn"]
    assert calls, "restore_fn is no longer called in app.py at all"

    # Every enclosing `if` test for each call site, by line span.
    def _guards_of(call):
        out = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.If):
                continue
            body_lines = [n.lineno for b in node.body for n in ast.walk(b)
                          if hasattr(n, "lineno")]
            if body_lines and min(body_lines) <= call.lineno <= max(body_lines):
                out.append(ast.unparse(node.test))
        return out

    for call in calls:
        guards = _guards_of(call)
        offending = [g for g in guards
                     if "load_session" in g or "have_session" in g
                     or "session is not None" in g or "stored is not None" in g]
        assert not offending, (
            f"restore_fn at app.py:{call.lineno} is gated on whether a stored session "
            f"exists ({offending}). A store with no viewer_session still holds the "
            "elements worth restoring."
        )


def test_an_empty_session_never_overwrites_a_stored_registration():
    """Restoring must not downgrade the disk.

    `_on_arms_restored` calls `_save_arms_affine_to_sdata()` on every restore. With
    no session the affine is None, the composed transform is the identity, and that
    identity used to be written over a real registration — losing it during a launch
    that only meant to read. Both overlay tabs must consult the element first.
    """
    import re

    src = Path(__file__).resolve().parent.parent / "src" / "palms" / "tabs"
    for name, element in (("tab_arms.py", "arms_he_image"),
                          ("tab_he_registration.py", "he_image")):
        text = (src / name).read_text()
        fn = re.search(r"def _save_\w*affine_to_sdata\(\):(.*?)\n    def ", text, re.S)
        assert fn, f"{name}: could not find the affine-saving function"
        body = fn.group(1)
        assert "_load_affine_from_sdata_element" in body, (
            f"{name}: _save_*_affine_to_sdata writes a transform without first checking "
            f"whether {element} already carries one. An identity composed from an empty "
            "state must not overwrite a stored registration."
        )
        # The *call*, not the import — `set_transformation` is imported at the top
        # of the same function, which would otherwise always sort first.
        guard_at = body.index("_load_affine_from_sdata_element(sdata")
        write_at = body.index("set_transformation(sdata")
        assert guard_at < write_at, (
            f"{name}: the stored-transform check must happen before set_transformation"
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
