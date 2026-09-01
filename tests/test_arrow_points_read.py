"""The pyarrow-filesystem points patch, and the tripwire that retires it.

The patch is the only compatibility shim in this repo that reaches into a
*private* third-party module path, so it is tested from two directions: that it
does what it claims, and that CI says so the moment it stops being needed.
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from palms.utils.arrow_points_read import arrow_points_reader, _is_local  # noqa: E402


def _io_points():
    spatialdata = pytest.importorskip("spatialdata")
    try:
        import spatialdata._io.io_points as io_points
    except Exception:  # pragma: no cover - layout changed; the next test says so
        pytest.skip("spatialdata._io.io_points has moved")
    assert spatialdata  # silence the unused-import linter
    return io_points


# ── the tripwire ─────────────────────────────────────────────────────────────

def test_delete_this_patch_once_spatialdata_passes_a_filesystem():
    """The whole reason the patch exists is that ``_read_points`` names no
    filesystem, so dask gives it the fsspec reader, which never pushes a filter
    down. When upstream fixes that, this test fails — and the fix is to delete
    ``utils/arrow_points_read.py``, its call in ``loader._open_cache``, the
    ``== True`` note in ``transcripts.gene.tmpl``, and this file."""
    io_points = _io_points()
    source = inspect.getsource(io_points._read_points)

    assert "filesystem" not in source, (
        "spatialdata._read_points now names a filesystem itself — the PALMS "
        "patch is obsolete. Delete utils/arrow_points_read.py, its call in "
        "loader._open_cache, the '== True' rationale in transcripts.gene.tmpl, "
        "and this test file. See docs/transcript-read-pushdown.md."
    )


def test_the_private_path_the_patch_needs_still_exists():
    """A silent degrade is the right *failure* mode and the wrong thing to
    discover months later: if spatialdata moves this, the viewer keeps working
    and quietly loses the 4.1x."""
    io_points = _io_points()
    assert callable(getattr(io_points, "read_parquet", None))


# ── what the patch does ──────────────────────────────────────────────────────

def test_the_patch_asks_for_the_pyarrow_filesystem():
    io_points = _io_points()
    seen = {}
    stock = io_points.read_parquet
    io_points.read_parquet = lambda path, *a, **k: seen.update(k) or "frame"
    try:
        with arrow_points_reader() as active:
            assert active
            io_points.read_parquet("/tmp/store/points.parquet")
    finally:
        io_points.read_parquet = stock

    assert seen.get("filesystem") == "arrow"


def test_the_stock_reader_is_restored_afterwards():
    io_points = _io_points()
    stock = io_points.read_parquet
    with arrow_points_reader():
        assert io_points.read_parquet is not stock
    assert io_points.read_parquet is stock


def test_the_stock_reader_is_restored_even_if_the_read_raises():
    io_points = _io_points()
    stock = io_points.read_parquet
    with pytest.raises(ValueError):
        with arrow_points_reader():
            raise ValueError("read failed")
    assert io_points.read_parquet is stock


def test_a_remote_store_is_left_on_fsspec():
    """``_read_points`` prefixes remote stores with ``simplecache::``, a chained
    fsspec URL the pyarrow filesystem cannot parse. Reading one that way would
    turn a working remote dataset into a crash."""
    io_points = _io_points()
    calls = []
    stock = io_points.read_parquet
    io_points.read_parquet = lambda path, *a, **k: calls.append((path, k))
    try:
        with arrow_points_reader():
            io_points.read_parquet("simplecache::https://host/points.parquet")
            io_points.read_parquet("s3://bucket/points.parquet")
    finally:
        io_points.read_parquet = stock

    assert [k.get("filesystem") for _p, k in calls] == [None, None]


@pytest.mark.parametrize("path, local", [
    ("/data/x/points.parquet", True),
    (Path("/data/x/points.parquet"), True),
    ("simplecache::https://host/p.parquet", False),
    ("s3://bucket/p.parquet", False),
    ("https://host/p.parquet", False),
])
def test_only_local_stores_are_redirected(path, local):
    assert _is_local(path) is local


def test_an_explicit_filesystem_outranks_the_patch():
    """Including a future spatialdata that passes one itself, which is what
    makes the wrapper inert rather than harmful before anyone deletes it."""
    io_points = _io_points()
    seen = {}
    stock = io_points.read_parquet
    io_points.read_parquet = lambda path, *a, **k: seen.update(k)
    try:
        with arrow_points_reader():
            io_points.read_parquet("/tmp/p.parquet", filesystem="fsspec")
    finally:
        io_points.read_parquet = stock

    assert seen["filesystem"] == "fsspec"


def test_a_failed_pyarrow_read_falls_back_rather_than_raising():
    io_points = _io_points()
    attempts = []

    def _flaky(path, *a, **k):
        attempts.append(k.get("filesystem"))
        if k.get("filesystem") == "arrow":
            raise RuntimeError("no pyarrow filesystem here")
        return "frame"

    stock = io_points.read_parquet
    io_points.read_parquet = _flaky
    try:
        with arrow_points_reader():
            assert io_points.read_parquet("/tmp/p.parquet") == "frame"
    finally:
        io_points.read_parquet = stock

    assert attempts == ["arrow", None]


# ── the other half of the fix ────────────────────────────────────────────────

def test_the_gene_template_states_is_gene_as_a_comparison():
    """A bare boolean term makes dask's ``_DNF.extract_pq_filters`` return None,
    which collapses the whole conjunction — so the patch above buys nothing
    unless the template is written this way. The two are one fix."""
    from palms.utils.step_templates import builtin_text

    text = builtin_text("transcripts.gene")
    assert "(_transcripts['is_gene'] == True)" in text
    assert "& _transcripts['is_gene']\n" not in text
