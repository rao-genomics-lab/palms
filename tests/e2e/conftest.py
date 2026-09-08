"""Fixtures for the end-to-end suite: stage a dataset, launch the real app.

These are the only tests that build a ``napari.Viewer`` and run the application's
own startup path. Everything else in ``tests/`` either stubs the context or reads
source. See ``tests/data/README.md`` for what the fixture dataset is.

Three guards are autouse, and each closes a way this suite could waste a CI run
rather than fail one: an unexpected modal dialog blocks forever, a hung Qt call
blocks forever, and a launch that writes into the checked-in fixture corrupts
every later run.
"""
from __future__ import annotations

import faulthandler
import hashlib
import os
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parent.parent.parent
DEFAULT_DATASET = REPO / "tests" / "data" / "crop_7"

#: What the committed fixture records as its own location, and therefore what
#: staging has to rewrite. Kept in step with scripts/prepare_e2e_fixture.py.
PLACEHOLDER_PATH = "/palms-e2e-fixture/crop_7"

#: Seconds before faulthandler dumps every thread's stack and aborts. Generous:
#: a cold launch is ~20 s here and a GitHub runner is slower, so this is a
#: hang detector, not a performance budget.
DEADLINE = 900


def pytest_collection_modifyitems(items):
    """Mark everything under ``tests/e2e/`` ``e2e``, and skip it on macOS.

    A hook rather than a module-level ``pytestmark``: that attribute is only
    honoured in test modules and classes, so setting it in a conftest looks
    right and does nothing — the marks would silently not exist, and a new test
    module would silently not get them.

    macOS runs the suite under the cocoa plugin against a real window server,
    and this is the one suite that opens windows. Get the Linux leg stable
    first; the fixture is committed, so enabling it later is deleting a line.
    """
    darwin = pytest.mark.skipif(
        sys.platform == "darwin",
        reason="E2E is Linux-only for now; see tests/data/README.md")
    here = Path(__file__).parent
    for item in items:
        if here in Path(str(item.fspath)).parents:
            item.add_marker(pytest.mark.e2e)
            item.add_marker(darwin)


# ── the dataset ──────────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def e2e_source() -> Path:
    """The dataset to stage from: $PALMS_E2E_DATASET, else the committed fixture.

    An env var that does not resolve is an error, not a skip. A typo in a path
    must not read as "no fixture available" and quietly stop testing anything.
    """
    override = os.environ.get("PALMS_E2E_DATASET")
    if override:
        path = Path(override).expanduser()
        if not (path / "experiment.xenium").exists():
            pytest.fail(f"PALMS_E2E_DATASET={override!r} is not a Xenium dataset "
                        f"(no experiment.xenium)")
        return path
    if not (DEFAULT_DATASET / "experiment.xenium").exists():
        pytest.skip(f"no E2E fixture at {DEFAULT_DATASET} and PALMS_E2E_DATASET is unset")
    return DEFAULT_DATASET


@pytest.fixture(scope="session")
def _source_digest(e2e_source) -> dict[str, str]:
    """sha256 of every file in the source, so a stray write cannot go unnoticed.

    Launching writes beside the dataset — ``setup_logging`` drops ``palms.log``,
    ``cache_repair.repair()`` touches the store, every recorded step rewrites
    ``viewer_cache/prov_graph.json``, and the graph restore rewrites
    ``analysis.py``. "The source is read-only" is therefore a property to
    enforce, not a convention to state.
    """
    return {
        str(p.relative_to(e2e_source)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(e2e_source.rglob("*")) if p.is_file()
    }


def _stage(source: Path, into: Path) -> Path:
    """Copy the dataset and repoint the absolute paths it records at the copy.

    Without the repair, ``app.py`` re-emits ``preamble`` for the new
    ``data_path`` on launch and flags every descendant stale — the documented
    rename hazard, which would make every staleness assertion vacuous.
    """
    from palms.scripts.rename_dataset import (
        Report, choose_graph_items, infer_old_path, read_sidecar_items,
        read_session_attrs, repair,
    )

    staged = into / source.name
    shutil.copytree(source, staged, symlinks=True)
    cache = staged / "sdata_cached.zarr"
    items = choose_graph_items(
        read_sidecar_items(staged),
        (read_session_attrs(cache).get("prov_graph") or []) if cache.exists() else [],
    )
    old = infer_old_path(items) or PLACEHOLDER_PATH
    repair(staged, old, str(staged),
           Report(old_path=Path(old), new_path=staged), dry_run=False)
    return staged


@pytest.fixture
def e2e_dataset(tmp_path, e2e_source, _source_digest) -> Path:
    """A staged, path-repaired copy of the dataset, in tmp_path."""
    yield _stage(e2e_source, tmp_path)
    _assert_source_untouched(e2e_source, _source_digest)


def _assert_source_untouched(source: Path, digest: dict[str, str]) -> None:
    changed = {
        str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(source.rglob("*")) if p.is_file()
    }
    if changed != digest:
        added = sorted(set(changed) - set(digest))
        removed = sorted(set(digest) - set(changed))
        edited = sorted(k for k in set(changed) & set(digest)
                        if changed[k] != digest[k])
        pytest.fail(f"the source dataset {source} was modified — "
                    f"added={added} removed={removed} edited={edited}")


# ── guards ───────────────────────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def no_dialogs(monkeypatch):
    """Any modal is a test failure, raised immediately.

    ``reporting._headless()`` already suppresses the modals *it* raises, but
    nothing suppresses a ``QFileDialog`` or a ``QMessageBox`` from a tab or from
    ``loader._ask_corrupt_cache``. Unguarded, one of those blocks the runner
    until the job times out — a lost run rather than a named failure.
    """
    from qtpy.QtWidgets import QDialog, QFileDialog, QMessageBox

    def _refuse(name):
        def _raise(*args, **kwargs):
            raise AssertionError(f"the app opened a modal dialog: {name}")
        return _raise

    for cls, names in (
        (QMessageBox, ("exec", "exec_", "question", "warning", "critical", "information")),
        (QFileDialog, ("getExistingDirectory", "getOpenFileName", "getSaveFileName",
                       "getOpenFileNames")),
        (QDialog, ("exec", "exec_")),
    ):
        for name in names:
            if hasattr(cls, name):
                monkeypatch.setattr(cls, name, _refuse(f"{cls.__name__}.{name}"),
                                    raising=False)


@pytest.fixture(autouse=True)
def deadline():
    """Abort with every thread's stack rather than hanging the runner.

    stdlib, so it needs no plugin — and unlike pytest-timeout it prints where
    the process actually was, which is the diagnostic a Qt hang otherwise denies.
    """
    faulthandler.dump_traceback_later(DEADLINE, exit=True)
    yield
    faulthandler.cancel_dump_traceback_later()


# ── the running application ──────────────────────────────────────────────────

@dataclass
class Session:
    """One launch of the application, and the verbs a test needs on it."""

    viewer: Any
    ctx: Any
    app: dict
    rig: Any
    no_cache: bool = False
    _closed: bool = field(default=False, repr=False)

    @property
    def graph(self):
        return self.ctx.state["prov_graph"]

    def obs_columns(self) -> set[str]:
        return set(self.ctx.adata.obs.columns)

    def settle(self, timeout: float = 120.0) -> None:
        """Run the event loop until no background worker is outstanding.

        Both halves of the launch matter here. ``create_viewer`` returns while
        tab restores are still running in ``thread_worker``s — the H&E pyramid
        is one — so a test that asserted immediately would be racing them.

        And at teardown it is not merely tidiness: a worker whose ``returned``
        callback lands *after* Qt has deleted the widgets it writes into raises
        inside a slot, and PyQt6 turns an exception in a slot into ``qFatal()``.
        The symptom is a core dump with no test name, which is precisely what
        this suite must never produce.
        """
        from superqt.utils import WorkerBase

        from palms.testing.rig import process_events
        deadline = time.time() + timeout
        while WorkerBase._worker_set and time.time() < deadline:
            process_events(0.05)
        if WorkerBase._worker_set:
            raise TimeoutError(
                f"{len(WorkerBase._worker_set)} background worker(s) still running "
                f"after {timeout}s: {WorkerBase._worker_set}")
        # The workers are done; their queued `returned` signals are not yet
        # delivered. Two more turns of the loop hand them to the GUI thread.
        process_events(0.05)
        process_events(0.05)

    def shutdown(self) -> None:
        """Exactly what closing the window does, in the production order.

        Both halves are the app's own functions: the layer snapshot has to
        happen while the Qt objects are alive, the store write after. A test
        that inlined either would be asserting against its own copy of the
        sequence rather than the one that ships.
        """
        from palms.app import shutdown_viewer
        self.app["snapshot_layers"]()
        shutdown_viewer(self.app, no_cache=self.no_cache)

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        from palms.testing.rig import process_events
        try:
            self.settle()
        except TimeoutError as exc:  # a hung worker must name itself, not hang
            print(f"  ! {exc}")
        # The linked UMAP scatter is a second top-level window; leaving it open
        # keeps Qt objects alive across the next launch in the same process.
        umap = getattr(self.ctx.umap_viewer, "_viewer", None) if self.ctx else None
        for window in (umap, self.viewer):
            try:
                window.close()
            except Exception:  # noqa: BLE001 - teardown must not mask a failure
                pass
        process_events(0.05)
        _drain()


def _drain():
    """Deliver deferred deletes, so one launch's widgets die before the next."""
    from qtpy.QtCore import QCoreApplication, QEvent
    from qtpy.QtWidgets import QApplication
    QApplication.processEvents()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    QApplication.processEvents()


def _launcher():
    """A ``launch(dataset)`` factory plus the teardown for whatever it started."""
    sessions: list[Session] = []

    def _launch(dataset: Path, no_cache: bool = False) -> Session:
        from palms.app import create_viewer
        from palms.testing.rig import Rig
        viewer, ctx, app_state = create_viewer(dataset, no_cache=no_cache)
        session = Session(viewer, ctx, app_state,
                          Rig.from_app(viewer, ctx, app_state), no_cache=no_cache)
        sessions.append(session)
        session.settle()
        return session

    yield _launch

    for session in reversed(sessions):
        session.close()


@pytest.fixture
def launch(qapp):
    """Factory: ``launch(dataset)`` -> a running :class:`Session`.

    A factory rather than a fixture value because the point of most of these
    tests is the *second* launch — the one that has to find what the first
    persisted.
    """
    yield from _launcher()


@pytest.fixture
def app(e2e_dataset, launch) -> Session:
    """A single running session against the staged dataset."""
    return launch(e2e_dataset)


# ── one shared, read-only session ────────────────────────────────────────────
# A launch is ~30 s, and the tests that only *inspect* a launched application
# have no reason to pay for one each. They must not mutate it — anything that
# writes takes `app` instead, which is function-scoped and gets its own copy.

@pytest.fixture(scope="module")
def _module_dataset(tmp_path_factory, e2e_source, _source_digest) -> Path:
    staged = _stage(e2e_source, tmp_path_factory.mktemp("e2e-module"))
    yield staged
    _assert_source_untouched(e2e_source, _source_digest)


@pytest.fixture(scope="module")
def _module_launch(qapp):
    yield from _launcher()


@pytest.fixture(scope="module")
def launched(_module_dataset, _module_launch) -> Session:
    """One launched application, shared by every read-only test in the module."""
    return _module_launch(_module_dataset)
