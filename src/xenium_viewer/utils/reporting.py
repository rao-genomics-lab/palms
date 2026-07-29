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


def _surface(kind: str, operation: str, message: str, show_modal: bool) -> None:
    """Hand off to the GUI thread, or do nothing when there is no GUI."""
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
    if not _failures:
        return "No write failures this session."
    kinds: dict[str, int] = {}
    for failure in _failures:
        kinds[failure["kind"]] = kinds.get(failure["kind"], 0) + 1
    detail = ", ".join(f"{count} {kind}" for kind, count in sorted(kinds.items()))
    return f"{len(_failures)} write failure(s) this session ({detail})."


def reset_failures() -> None:
    """Test helper — also used when switching datasets."""
    _failures.clear()
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
