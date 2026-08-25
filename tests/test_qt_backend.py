"""Which Qt binding the app runs on, and why it is stated rather than inferred.

A green test suite on the wrong Qt backend is the failure this guards. It is not
hypothetical — it happened, in CI, on the very branch that migrated to PyQt6:

    napari 0.8.0 | qtpy PyQt5 5.15.15
    AssertionError: expected PyQt6, got PyQt5

with ``pyqt6`` sitting in ``environment.yml`` and PyQt6 6.8.1 installed. Two
things combined to produce it, and this module pins both.

**qtpy resolves in the order PyQt5, PySide2, PyQt6, PySide6.** So any environment
that merely *contains* PyQt5 runs on it, whatever else is installed and whatever
the environment file asked for.

**Something quietly installed PyQt5 beside PyQt6.** conda-forge's ``matplotlib``
metapackage bundles a Qt binding, and which one depends on the version the solver
picks: 3.9.1 depends on ``pyqt >=5.10`` (Qt5) while 3.9.3+ depend on ``pyside6``.
CI landed on 3.9.1. Nothing in this project needs that metapackage — scanpy,
squidpy and matplotlib-scalebar all depend on ``matplotlib-base``, and the app
supplies its own binding.

Both tests are pure stdlib: no Qt, no dataset, no network.
"""

from __future__ import annotations

import os
import re
from importlib.util import find_spec
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
ENVIRONMENT_YML = REPO / "environment.yml"


def test_the_qt_backend_is_stated_not_inherited():
    """``import xenium_viewer`` pins QT_API when PyQt6 is available.

    Skipped rather than asserted-around on an env without PyQt6: there the
    correct behaviour is to leave qtpy alone, which the next test covers.
    """
    if find_spec("PyQt6") is None:
        pytest.skip("no PyQt6 in this environment; see the next test")

    import xenium_viewer  # noqa: F401  (the import is the behaviour under test)

    assert os.environ.get("QT_API") == "pyqt6"

    import qtpy
    assert qtpy.API_NAME == "PyQt6", (
        f"qtpy chose {qtpy.API_NAME} despite QT_API=pyqt6 — a second binding is "
        "winning, which is exactly the CI failure this pin exists to prevent"
    )


#: qtpy's env-var spellings and the module each one names.
_BINDINGS = {"pyqt5": "PyQt5", "pyqt6": "PyQt6",
             "pyside2": "PySide2", "pyside6": "PySide6"}


def test_qt_api_never_names_a_binding_that_is_not_installed():
    """The pin must not turn a working install into an ImportError.

    An environment with only PySide6 is legitimate, so the pin is conditional on
    PyQt6 actually being importable. This checks the resulting invariant in *any*
    environment rather than the branch taken in this one.

    Note the assertion is deliberately not "QT_API is unset" — qtpy writes the
    variable back with whatever binding it selected as soon as it is imported, so
    by the time a test runs it is populated regardless of what this package did.
    """
    import xenium_viewer  # noqa: F401

    api = os.environ.get("QT_API")
    if not api:
        return                      # nothing imported qtpy yet; nothing to check

    module = _BINDINGS.get(api.lower())
    assert module is not None, f"QT_API names an unknown binding: {api!r}"
    assert find_spec(module) is not None, (
        f"QT_API={api!r} but {module} is not installed"
    )


def test_environment_asks_for_matplotlib_base_not_the_metapackage():
    """A source guard: the metapackage drags in whichever Qt binding it fancies.

    Reintroducing a bare ``matplotlib`` dependency reintroduces a second Qt
    binding on some solver runs and not others — which is the worst kind of
    dependency bug, because it reproduces only sometimes.
    """
    lines = [
        ln.split("#", 1)[0].strip()
        for ln in ENVIRONMENT_YML.read_text().splitlines()
    ]
    deps = {re.sub(r"[<>=!].*$", "", ln[2:]).strip() for ln in lines if ln.startswith("- ")}

    assert "matplotlib-base" in deps, "environment.yml should depend on matplotlib-base"
    assert "matplotlib" not in deps, (
        "environment.yml asks for the `matplotlib` metapackage. It bundles a Qt "
        "binding whose identity depends on the version resolved (3.9.1 -> pyqt5, "
        "3.9.3+ -> pyside6), so a second binding can appear and qtpy will prefer "
        "PyQt5 over PyQt6. Use matplotlib-base."
    )


def test_environment_ships_the_unversioned_libglx_name():
    """A source guard for the crash that nearly shelved the Qt6 migration.

    conda's ``libglx`` provides only ``libGLX.so.0``. PyOpenGL's loader
    (``OpenGL/platform/ctypesloader.py::_loadLibraryPosix``) tries the
    *unversioned* ``libGLX.so`` first, so without ``libglx-devel`` it misses the
    environment and ``dlopen``s the host's ``/usr/lib/.../libGLX.so`` with
    ``RTLD_GLOBAL``. Two different glvnd builds then share one process, Qt's GLX
    plugin resolves ``glX*`` across both, ``glXGetVisualFromFBConfig()`` returns
    NULL for a config the other library allocated, and Qt6 calls ``qFatal`` —
    ``SIGABRT``, no traceback, nothing to debug.

    Measured here: napari aborted on startup under Qt6 over a remote X display
    until this package was added, and started rendering immediately once it was.
    PyQt5 maps the same two copies and merely tolerates them, so dropping this
    dependency would look harmless right up until the backend changed.
    """
    lines = [
        ln.split("#", 1)[0].strip()
        for ln in ENVIRONMENT_YML.read_text().splitlines()
    ]
    deps = {re.sub(r"[<>=!].*$", "", ln[2:]).strip() for ln in lines if ln.startswith("- ")}

    assert "libglx-devel" in deps, (
        "environment.yml must depend on libglx-devel. Without it PyOpenGL loads "
        "the host's libGLX beside conda's, and Qt6 aborts with "
        '"Could not initialize GLX" and no traceback.'
    )
