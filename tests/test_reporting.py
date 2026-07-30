"""Write failures must be recorded somewhere a user can find them.

Every failure used to go to stdout, which a GUI user never reads, and the only
dialog was for permission errors — shown once per *process* via a module-level
flag, so the second failure and everything after it was invisible. When the
cache was being corrupted, the warnings that would have explained it were lost.

Run headless:  QT_QPA_PLATFORM=offscreen pytest tests/test_reporting.py
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from xenium_viewer.utils import reporting  # noqa: E402
from xenium_viewer.utils.reporting import (  # noqa: E402
    LOG_FILENAME, classify, failure_summary, failures, get_logger,
    report_write_failure, reset_failures, setup_logging,
)


@pytest.fixture(autouse=True)
def _clean():
    reset_failures()
    reporting._configured_path = None
    for handler in list(reporting._log.handlers):
        reporting._log.removeHandler(handler)
        handler.close()
    yield
    reset_failures()


# ── the log file ─────────────────────────────────────────────────────────────

def test_setup_creates_a_log_beside_the_dataset(tmp_path):
    path = setup_logging(tmp_path)
    assert path == tmp_path / LOG_FILENAME
    assert path.exists()


def test_warnings_reach_the_file(tmp_path):
    setup_logging(tmp_path)
    get_logger("xenium_viewer.tests").warning("could not persist adata table: disk full")
    assert "disk full" in (tmp_path / LOG_FILENAME).read_text()


def test_failures_are_logged_with_a_traceback(tmp_path):
    setup_logging(tmp_path)
    try:
        raise OSError("No space left on device")
    except OSError as exc:
        report_write_failure(exc, "clustering data")

    text = (tmp_path / LOG_FILENAME).read_text()
    assert "clustering data" in text
    assert "Traceback" in text          # the part stdout never carried


def test_setup_is_idempotent(tmp_path):
    setup_logging(tmp_path)
    before = len(reporting._log.handlers)
    setup_logging(tmp_path)
    assert len(reporting._log.handlers) == before


def test_switching_datasets_moves_the_log(tmp_path):
    first, second = tmp_path / "a", tmp_path / "b"
    first.mkdir()
    second.mkdir()
    setup_logging(first)
    setup_logging(second)

    get_logger("xenium_viewer.tests").warning("after the switch")
    assert "after the switch" in (second / LOG_FILENAME).read_text()
    assert "after the switch" not in (first / LOG_FILENAME).read_text()


def test_a_read_only_dataset_directory_is_not_fatal(tmp_path, monkeypatch):
    """A read-only mount is a legitimate configuration."""
    def boom(*args, **kwargs):
        raise OSError("Read-only file system")

    monkeypatch.setattr(reporting, "RotatingFileHandler", boom, raising=False)
    monkeypatch.setattr(
        "logging.handlers.RotatingFileHandler",
        lambda *a, **k: (_ for _ in ()).throw(OSError("Read-only file system")),
    )
    assert setup_logging(tmp_path) is None
    # console logging still works
    get_logger("xenium_viewer.tests").warning("still reported")


def test_the_log_lives_outside_the_zarr_store(tmp_path):
    """Anything inside the store makes zarr's hierarchy walk warn."""
    setup_logging(tmp_path)
    assert (tmp_path / LOG_FILENAME).parent.name != "sdata_cached.zarr"


# ── classification drives what the user is told ──────────────────────────────

@pytest.mark.parametrize("exc,expected", [
    (PermissionError("nope"), "permission"),
    (OSError("Permission denied"), "permission"),
    (OSError("attempt to write a read-only database"), "permission"),
    (OSError("No space left on device"), "space"),
    (OSError("not enough free space to rebuild the cache"), "space"),
    (ValueError("something else entirely"), "other"),
])
def test_classify(exc, expected):
    assert classify(exc) == expected


def test_errno_based_classification():
    exc = OSError("failed")
    exc.errno = 28
    assert classify(exc) == "space"
    exc.errno = 13
    assert classify(exc) == "permission"


# ── the running tally ────────────────────────────────────────────────────────

def test_failures_accumulate(tmp_path):
    setup_logging(tmp_path)
    report_write_failure(OSError("a"), "ROIs")
    report_write_failure(PermissionError("b"), "annotations")

    recorded = failures()
    assert [f["operation"] for f in recorded] == ["ROIs", "annotations"]
    assert [f["kind"] for f in recorded] == ["other", "permission"]


def test_summary_reads_cleanly(tmp_path):
    setup_logging(tmp_path)
    assert "No write failures" in failure_summary()
    report_write_failure(OSError("x"), "ROIs")
    report_write_failure(PermissionError("y"), "table")
    summary = failure_summary()
    assert "2 write failure(s)" in summary
    assert "1 other" in summary and "1 permission" in summary


def test_every_failure_is_recorded_not_just_the_first(tmp_path):
    """The old dialog suppressed itself after one, and nothing else surfaced."""
    setup_logging(tmp_path)
    for i in range(5):
        report_write_failure(PermissionError(f"attempt {i}"), "table")

    assert len(failures()) == 5
    text = (tmp_path / LOG_FILENAME).read_text()
    for i in range(5):
        assert f"attempt {i}" in text


def test_reporting_never_raises_without_a_gui(tmp_path):
    setup_logging(tmp_path)
    report_write_failure(PermissionError("no qt here"), "table")   # must not raise


# ── the modal must never block a process with nobody at the keyboard ─────────

@pytest.mark.parametrize("platform,headless", [
    ("offscreen", True), ("minimal", True), ("vnc", True),
    ("xcb", False), ("wayland", False), ("", False),
])
def test_headless_platforms_are_recognised(monkeypatch, platform, headless):
    from xenium_viewer.utils.reporting import _headless
    monkeypatch.setenv("QT_QPA_PLATFORM", platform)
    assert _headless() is headless


def test_no_modal_is_raised_when_qt_is_headless(tmp_path, monkeypatch, qapp):
    """Regression: this hung the suite for an hour.

    ``QApplication.instance() is not None`` was the only guard, so once a
    fixture created one before the tests that inject write failures, the
    permission dialog's ``exec_()`` blocked with no event loop and no human able
    to dismiss it. The failure must still be logged and notified — only the
    modal is suppressed.
    """
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    import qtpy.QtWidgets as qtw
    import xenium_viewer.utils.reporting as reporting

    constructed = []

    class _Tripwire:
        def __init__(self, *a, **k):
            constructed.append(True)
            raise AssertionError("a modal dialog would have blocked here")

    monkeypatch.setattr(qtw, "QMessageBox", _Tripwire)
    setup_logging(tmp_path)

    reporting._surface("permission", "table", "read-only fs", show_modal=True)

    assert constructed == []
    # the failure is still recorded — suppressing the modal must not hide it
    reporting.report_write_failure(PermissionError("read-only fs"), "table")
    assert failures()[-1]["operation"] == "table"


# ── a step that ran but could not be recorded ────────────────────────────────

def test_a_recording_failure_reaches_the_log_with_its_node_id(tmp_path):
    """It used to go to ``warnings.warn`` — a terminal a GUI user never reads,
    and only once per unique message under the default ``once`` filter."""
    setup_logging(tmp_path)
    try:
        raise KeyError("unknown dependency: normalize")
    except KeyError as exc:
        reporting.report_recording_failure("nhood:leiden_r1.0", exc)

    text = (tmp_path / LOG_FILENAME).read_text()
    assert "nhood:leiden_r1.0" in text
    assert "normalize" in text
    assert "Traceback" in text


def test_recording_failures_accumulate_separately_from_write_failures(tmp_path):
    """Different things: one lost the code, the other lost the data."""
    setup_logging(tmp_path)
    reporting.report_write_failure(OSError("disk"), "ROIs")
    reporting.report_recording_failure("rank_genes:k", KeyError("missing dep"))

    assert [f["node_id"] for f in reporting.recording_failures()] == ["rank_genes:k"]
    assert [f["operation"] for f in failures()] == ["ROIs"]


def test_every_recording_failure_is_kept(tmp_path):
    setup_logging(tmp_path)
    for i in range(4):
        reporting.report_recording_failure(f"node:{i}", KeyError("missing dep"))

    assert len(reporting.recording_failures()) == 4


def test_the_health_line_mentions_unrecorded_steps(tmp_path):
    """The Cache tab's health section is where a user looks for "did anything
    go wrong?" — a notebook missing a step belongs there."""
    setup_logging(tmp_path)
    assert "could not be recorded" not in failure_summary()

    reporting.report_recording_failure("clustering:k", KeyError("missing dep"))
    summary = failure_summary()
    assert "No write failures" in summary
    assert "1 step(s) could not be recorded" in summary


def test_a_recording_failure_never_raises_a_modal(tmp_path, monkeypatch, qapp):
    """The analysis succeeded — there is nothing for the user to decide, and a
    dialog here would block a worker thread's callback."""
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    import qtpy.QtWidgets as qtw

    class _Tripwire:
        def __init__(self, *a, **k):
            raise AssertionError("a modal dialog would have blocked here")

    monkeypatch.setattr(qtw, "QMessageBox", _Tripwire)
    setup_logging(tmp_path)

    reporting.report_recording_failure("clustering:k", KeyError("missing dep"))
    assert reporting.recording_failures()[-1]["node_id"] == "clustering:k"


def test_the_recorder_does_not_fall_back_to_warnings(tmp_path):
    """Source guard: ``warnings.warn`` in the recorder is what this replaces."""
    src = Path(__file__).resolve().parent.parent / "src" / "xenium_viewer"
    text = (src / "tabs" / "_helpers.py").read_text()
    assert "warnings.warn" not in text


# ── the shim keeps old call sites working ────────────────────────────────────

def test_permission_dialog_shim_routes_to_the_reporter(tmp_path):
    from xenium_viewer.utils.adata_persistence import _maybe_show_permission_dialog

    setup_logging(tmp_path)
    _maybe_show_permission_dialog(PermissionError("legacy path"), "ROI shapes")

    assert failures()[-1]["operation"] == "ROI shapes"
    assert "legacy path" in (tmp_path / LOG_FILENAME).read_text()


def test_no_persistence_module_still_prints_warnings_instead_of_logging():
    """Source guard: cache write paths must log, so the file captures them."""
    import re

    src = Path(__file__).resolve().parent.parent / "src" / "xenium_viewer"
    offenders = []
    for name in ("utils/adata_persistence.py", "utils/session.py"):
        text = (src / name).read_text()
        for match in re.finditer(r'print\(f?"[^"]*Warning[^"]*"', text):
            line = text[: match.start()].count("\n") + 1
            offenders.append(f"{name}:{line}")
    assert offenders == [], (
        "cache write paths must use log.warning / report_write_failure so the "
        "log file captures them: " + ", ".join(offenders)
    )


def test_log_level_keeps_debug_out_of_the_console(tmp_path):
    setup_logging(tmp_path, level=logging.INFO)
    stream = [h for h in reporting._log.handlers
              if isinstance(h, logging.StreamHandler)
              and not hasattr(h, "baseFilename")]
    assert stream and stream[0].level == logging.INFO


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
