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
ENVIRONMENT_LINUX_YML = REPO / "environment-linux.yml"
INSTALL_SH = REPO / "scripts" / "install.sh"


def _env_deps(path: Path) -> set[str]:
    """Top-level package names from a conda env file, comments and pins stripped."""
    lines = [ln.split("#", 1)[0].strip() for ln in path.read_text().splitlines()]
    return {
        re.sub(r"[<>=!].*$", "", ln[2:]).strip()
        for ln in lines
        if ln.startswith("- ")
    }


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
    deps = _env_deps(ENVIRONMENT_YML)

    assert "matplotlib-base" in deps, "environment.yml should depend on matplotlib-base"
    assert "matplotlib" not in deps, (
        "environment.yml asks for the `matplotlib` metapackage. It bundles a Qt "
        "binding whose identity depends on the version resolved (3.9.1 -> pyqt5, "
        "3.9.3+ -> pyside6), so a second binding can appear and qtpy will prefer "
        "PyQt5 over PyQt6. Use matplotlib-base."
    )


def test_the_linux_overlay_ships_the_unversioned_libglx_name():
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

    It lives in the overlay rather than ``environment.yml`` because it is a
    linux-only conda-forge package — see the next test.
    """
    assert "libglx-devel" in _env_deps(ENVIRONMENT_LINUX_YML), (
        "environment-linux.yml must depend on libglx-devel. Without it PyOpenGL "
        "loads the host's libGLX beside conda's, and Qt6 aborts with "
        '"Could not initialize GLX" and no traceback.'
    )


def test_the_main_environment_stays_solvable_off_linux():
    """The other half of the split: the linux-only package must not creep back.

    ``libglx-devel`` has no ``osx-arm64``/``osx-64`` build, so while it sat in
    ``environment.yml`` the solve failed on macOS before it reached the Qt stack —
    ``conda env create -f environment.yml`` was simply impossible there. conda env
    files have no platform selectors (``# [linux]`` is a conda-*build* feature),
    which is why the fix is an overlay file rather than a conditional line.
    """
    assert "libglx-devel" not in _env_deps(ENVIRONMENT_YML), (
        "libglx-devel is back in environment.yml. It is a linux-only package, so "
        "this makes the file unsolvable on macOS. Put it in environment-linux.yml."
    )


class TestLibglxCollisionWarning:
    """``utils.gl_check`` turns a traceback-free SIGABRT into a sentence.

    The overlay split means a Linux user can now create the environment without
    ``libglx-devel`` — ``conda env create -f environment.yml`` alone does it. That
    is only safe because this check exists: the failure it replaces is napari
    calling ``qFatal`` during ``import``, which produces no Python traceback and
    no clue about the cause.

    Both halves of the collision must be present. Warning on a missing package
    alone would fire on every correctly-working machine that simply has no host
    ``libglx-dev`` installed, and a warning that is usually wrong is a warning
    people learn to scroll past.
    """

    @staticmethod
    def _prefix_without_the_name(tmp_path):
        (tmp_path / "lib").mkdir()
        return tmp_path

    def test_warns_when_the_env_lacks_it_and_the_host_has_it(self, tmp_path):
        from xenium_viewer.utils import gl_check

        host = tmp_path / "host-libGLX.so"
        host.write_text("")
        msg = gl_check.libglx_collision_message(
            prefix=self._prefix_without_the_name(tmp_path),
            host_paths=(str(host),),
            platform="linux",
        )
        assert msg is not None
        # The point of the message is the command, not the diagnosis.
        assert "environment-linux.yml" in msg
        assert str(host) in msg

    def test_silent_when_the_env_has_the_unversioned_name(self, tmp_path):
        from xenium_viewer.utils import gl_check

        prefix = self._prefix_without_the_name(tmp_path)
        (prefix / "lib" / "libGLX.so").write_text("")
        host = tmp_path / "host-libGLX.so"
        host.write_text("")
        assert gl_check.libglx_collision_message(
            prefix=prefix, host_paths=(str(host),), platform="linux"
        ) is None

    def test_silent_when_there_is_no_host_copy_to_collide_with(self, tmp_path):
        from xenium_viewer.utils import gl_check

        assert gl_check.libglx_collision_message(
            prefix=self._prefix_without_the_name(tmp_path),
            host_paths=(str(tmp_path / "nope.so"),),
            platform="linux",
        ) is None

    def test_silent_on_macos_even_in_the_defective_shape(self, tmp_path):
        """macOS has no GLX at all — Qt uses cocoa and PyOpenGL uses the framework.

        Without this branch the warning would fire on every Mac, telling users to
        install a package that has no macOS build.
        """
        from xenium_viewer.utils import gl_check

        host = tmp_path / "host-libGLX.so"
        host.write_text("")
        assert gl_check.libglx_collision_message(
            prefix=self._prefix_without_the_name(tmp_path),
            host_paths=(str(host),),
            platform="darwin",
        ) is None

    def test_the_checked_paths_are_the_unversioned_name(self):
        """``libGLX.so.0`` is present on every working box; it proves nothing.

        ``ctypes.util.find_library('GLX')`` returns exactly that versioned soname,
        which is why this module globs known paths instead of using it. A check
        that matched ``.so.0`` would warn universally.
        """
        from xenium_viewer.utils import gl_check

        assert gl_check.HOST_LIBGLX_PATHS
        for path in gl_check.HOST_LIBGLX_PATHS:
            assert path.endswith("/libGLX.so"), path


def test_the_installer_applies_the_linux_overlay():
    """An overlay nothing applies is an overlay nobody gets.

    The split only works because ``scripts/install.sh`` branches on ``uname -s``
    and runs ``env update`` with the overlay on Linux. If that reference is lost,
    every Linux install silently reverts to the pre-``libglx-devel`` state, whose
    symptom is a SIGABRT with no traceback.
    """
    assert INSTALL_SH.exists(), "scripts/install.sh is the OS-aware install path"
    script = INSTALL_SH.read_text()

    assert "environment-linux.yml" in script, (
        "scripts/install.sh no longer references environment-linux.yml, so the "
        "Linux GL overlay would never be applied"
    )
    assert "env update" in script, (
        "scripts/install.sh must apply the overlay with `env update`"
    )
    assert re.search(r"uname\s+-s", script), (
        "scripts/install.sh must branch on the OS; the overlay is Linux-only"
    )
