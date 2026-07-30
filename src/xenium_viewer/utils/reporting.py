"""Logging and user-visible reporting of failures.

Every write failure in the viewer used to go to stdout, which a GUI user never
sees, and the only dialog was for permission errors — shown **once per process**
via a module-level flag, so the second failure and everything after it was
invisible. When the zarr cache was being corrupted, the warnings that would have
explained it went to a terminal nobody was reading.

Three pieces, deliberately small:

* a rotating log file per dataset, so the next occurrence is captured with a
  traceback whether or not anyone was watching;
* :func:`report_write_failure`, which always logs and shows a **non-modal**
  napari notification, marshalled to the GUI thread — the old dialog could be
  constructed from a ``thread_worker``, which is a real Qt violation;
* a running tally the Cache tab can display, so failures are visible in
  aggregate without a popup per event.
"""

from __future__ import annotations

import logging
import os
import sys
import threading
from pathlib import Path
from typing import Optional

LOG_FILENAME = "xenium_viewer.log"
LOGGER_NAME = "xenium_viewer"

_log = logging.getLogger(LOGGER_NAME)

_configured_path: Optional[Path] = None
_lock = threading.Lock()

# Failures since launch, for the Cache tab's health line.
_failures: list[dict] = []

# Steps whose *recording* failed since launch. Separate from _failures: the
# analysis ran and its result is on screen; what was lost is the code for it.
_recording_failures: list[dict] = []

# Modal dialogs already shown, keyed by (dataset, error class). The previous
# implementation used a single process-wide bool, so one permission error
# suppressed every subsequent dialog for the life of the session.
_modal_shown: set[tuple[str, str]] = set()


def get_logger(name: str = LOGGER_NAME) -> logging.Logger:
    return logging.getLogger(name)


def log_path(data_path) -> Path:
    return Path(data_path) / LOG_FILENAME


def setup_logging(data_path, level: int = logging.INFO) -> Optional[Path]:
    """Attach a rotating file handler for *data_path*. Idempotent.

    Returns the log path, or ``None`` if the directory is not writable (a
    read-only dataset mount is a legitimate configuration — console output still
    works, and the write itself will fail with a clearer message than we could
    produce here).
    """
    global _configured_path
    from logging.handlers import RotatingFileHandler

    target = log_path(data_path)
    with _lock:
        if _configured_path == target:
            return target

        for handler in list(_log.handlers):
            _log.removeHandler(handler)
            handler.close()

        _log.setLevel(level)
        # Root already prints elsewhere in this app; don't double up.
        _log.propagate = False

        formatter = logging.Formatter(
            "%(asctime)s %(levelname)-7s %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        # Keep the existing stdout behaviour — the terminal is still useful when
        # someone *is* watching, and scripts depend on it.
        stream = logging.StreamHandler(sys.stdout)
        stream.setFormatter(logging.Formatter("%(message)s"))
        stream.setLevel(level)
        _log.addHandler(stream)

        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            file_handler = RotatingFileHandler(
                target, maxBytes=2 * 1024 * 1024, backupCount=3, encoding="utf-8",
            )
            file_handler.setFormatter(formatter)
            file_handler.setLevel(logging.DEBUG)
            _log.addHandler(file_handler)
        except OSError as exc:
            _log.warning("could not open the log file at %s: %s", target, exc)
            _configured_path = None
            return None

        _configured_path = target
        _log.info("Xenium Viewer log started (%s)", target)
        return target


def configured_log_path() -> Optional[Path]:
    return _configured_path


# ── failure reporting ────────────────────────────────────────────────────────

def classify(exc: BaseException) -> str:
    """'permission', 'space' or 'other' — what the user can actually do about it."""
    text = str(exc)
    if (isinstance(exc, PermissionError)
            or (isinstance(exc, OSError) and getattr(exc, "errno", None) in (13, 30))
            or "Permission denied" in text
            or "read-only" in text.lower()):
        return "permission"
    if (isinstance(exc, OSError) and getattr(exc, "errno", None) == 28) or (
            "No space left" in text or "not enough free space" in text):
        return "space"
    return "other"


_GUIDANCE = {
    "permission": (
        "The dataset folder may be read-only (e.g. a shared or mounted drive).\n\n"
        "To enable saving, copy the dataset to a writable location and reopen it, "
        "or launch with --no-cache to skip zarr persistence entirely."
    ),
    "space": (
        "The disk holding the dataset is full.\n\n"
        "Free some space and retry — your existing cache has not been modified."
    ),
}


def report_write_failure(exc: BaseException, operation: str = "data",
                         *, dataset=None) -> None:
    """Log a failed save and surface it without blocking the user.

    Always logs with a traceback. Shows a non-modal napari notification for every
    failure, and a modal dialog only for the classes where the user must act —
    once per (dataset, class), not once per process.
    """
    kind = classify(exc)
    _log.warning("could not save %s: %s", operation, exc,
                 exc_info=(type(exc), exc, exc.__traceback__))
    _failures.append({"operation": operation, "error": str(exc), "kind": kind})

    key = (str(dataset or _configured_path or ""), kind)
    show_modal = kind in _GUIDANCE and key not in _modal_shown
    if show_modal:
        _modal_shown.add(key)
    _surface(kind, operation, str(exc), show_modal)


def report_recording_failure(node_id: str, exc: BaseException) -> None:
    """Surface a step that ran but could not be recorded properly.

    The recorder degrades rather than aborting — a bug in provenance bookkeeping
    must never lose the user's analysis — but it used to degrade through
    ``warnings.warn``, which in a GUI process goes to a terminal nobody is
    reading (and only once per unique message, since the default filter is
    ``once``). The result was the failure mode this whole phase exists to
    prevent: a result on screen whose cell is missing from the notebook, with
    nothing said about it.

    Logs with a traceback, keeps the failure for the Cache tab's health line,
    and shows a non-modal napari warning. No dialog: the action succeeded, so
    there is nothing for the user to stop and decide.
    """
    _log.warning("could not record step %r: %s", node_id, exc,
                 exc_info=(type(exc), exc, exc.__traceback__))
    _recording_failures.append({"node_id": node_id, "error": str(exc)})
    _notify(
        f"Code recording degraded for '{node_id}': {exc}. The step ran and its "
        f"code was appended, but without its place in the provenance graph — "
        f"the exported notebook may not replay it in the right order."
    )


def recording_failures() -> list[dict]:
    """Steps whose recording failed since launch (most recent last)."""
    return list(_recording_failures)


def _notify(message: str) -> None:
    """Non-modal napari warning, on the GUI thread; a no-op with no GUI."""
    try:
        from qtpy.QtWidgets import QApplication
        if QApplication.instance() is None:
            return
        from superqt.utils import ensure_main_thread
    except Exception:
        return

    @ensure_main_thread
    def _show():
        try:
            from napari.utils.notifications import show_warning
            show_warning(message)
        except Exception:  # pragma: no cover - napari not always present
            pass

    _show()


def _event_loop_running() -> bool:
    """True when the main thread is actually inside a Qt event loop.

    ``QMessageBox.exec_()`` starts a *nested* event loop and returns only when
    something dismisses the dialog. In a process where nobody ever called
    ``app.exec_()`` there is no window manager, no user and no other loop to
    deliver that click, so it never returns.

    ``QThread.loopLevel()`` is the direct answer — 0 outside any event loop, ≥1
    inside one — and it is readable from a worker thread about the main thread,
    which matters because failures are reported from napari workers.
    """
    try:
        from qtpy.QtWidgets import QApplication
        app = QApplication.instance()
        if app is None:
            return False
        return app.thread().loopLevel() > 0
    except Exception:  # pragma: no cover - old Qt without loopLevel
        # Unknown → assume a loop is running. A missed suppression shows a
        # dialog nobody asked for; a wrong suppression is silence about a
        # failed save, which is worse.
        return True


def _headless() -> bool:
    """True when Qt is up but nobody can answer a dialog.

    ``QApplication.instance() is not None`` is not enough: a test run, a
    headless script or CI creates one with no event loop and no human, and
    ``QMessageBox.exec_()`` then blocks forever with nothing able to dismiss it.
    That is not hypothetical — it hung the test suite for an hour once a fixture
    started creating the QApplication before the tests that inject write
    failures. The notification still fires; only the modal is suppressed.

    Two independent signals, because the first one alone left the bug live on
    exactly the machine where it hurts most. ``QT_QPA_PLATFORM`` catches CI and
    explicitly headless runs — but a developer's desktop *has* a display, so
    ``conftest.py`` deliberately does not set it, and every such signal stayed
    unset while `pytest` still had no event loop. `pytest` then hung, for real,
    in ``test_persistence_safety.py``. The loop-level check is the one that
    actually describes the precondition: not "is there a screen" but "is anyone
    processing events".
    """
    platform = os.environ.get("QT_QPA_PLATFORM", "").split(":")[0].strip().lower()
    if platform in {"offscreen", "minimal", "minimalegl", "vnc"}:
        return True
    return not _event_loop_running()


def _surface(kind: str, operation: str, message: str, show_modal: bool) -> None:
    """Hand off to the GUI thread, or do nothing when there is no GUI."""
    try:
        from qtpy.QtWidgets import QApplication
        if QApplication.instance() is None:
            return
        from superqt.utils import ensure_main_thread
    except Exception:
        return
    if _headless():
        show_modal = False

    @ensure_main_thread
    def _show():
        try:
            from napari.utils.notifications import show_error
            show_error(f"Could not save {operation}: {message}")
        except Exception:  # pragma: no cover - napari not always present
            pass
        if not show_modal:
            return
        try:
            from qtpy.QtWidgets import QMessageBox
            box = QMessageBox()
            box.setIcon(QMessageBox.Warning)
            box.setWindowTitle("Could Not Save")
            box.setText(f"Could not save {operation}.")
            box.setInformativeText(_GUIDANCE[kind])
            box.setDetailedText(message)
            box.exec_()
        except Exception:  # pragma: no cover - defensive
            pass

    _show()


def failures() -> list[dict]:
    """Write failures recorded since launch (most recent last)."""
    return list(_failures)


def failure_summary() -> str:
    """One line for the Cache tab's health section."""
    if _failures:
        kinds: dict[str, int] = {}
        for failure in _failures:
            kinds[failure["kind"]] = kinds.get(failure["kind"], 0) + 1
        detail = ", ".join(f"{count} {kind}" for kind, count in sorted(kinds.items()))
        summary = f"{len(_failures)} write failure(s) this session ({detail})."
    else:
        summary = "No write failures this session."
    if _recording_failures:
        # A notebook missing a step is a health problem too, and this line is
        # the one place a user goes looking for "did anything go wrong?".
        summary += (f" {len(_recording_failures)} step(s) could not be recorded "
                    f"— see the log.")
    return summary


def reset_failures() -> None:
    """Test helper — also used when switching datasets."""
    _failures.clear()
    _recording_failures.clear()
    _modal_shown.clear()


def open_log_file() -> bool:
    """Open the log in the platform's default viewer. Returns False if there isn't one."""
    if _configured_path is None or not _configured_path.exists():
        return False
    try:
        if sys.platform == "darwin":
            os.system(f'open "{_configured_path}"')
        elif os.name == "nt":  # pragma: no cover - not a supported platform
            os.startfile(str(_configured_path))  # type: ignore[attr-defined]
        else:
            os.system(f'xdg-open "{_configured_path}" &')
        return True
    except Exception:  # pragma: no cover - defensive
        return False
