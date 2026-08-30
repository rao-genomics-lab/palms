"""Detect the libGLX double-load *before* importing napari aborts the process.

conda's ``libglx`` provides only ``libGLX.so.0``. PyOpenGL's loader
(``OpenGL/platform/ctypesloader.py::_loadLibraryPosix``) tries the *unversioned*
``libGLX.so`` first, so without ``libglx-devel`` it misses the environment and
``dlopen``s the host's ``/usr/lib/.../libGLX.so`` (Ubuntu's ``libglx-dev``) with
``RTLD_GLOBAL``. Two different glvnd builds then share one process, Qt's GLX
plugin resolves ``glX*`` across both, ``glXGetVisualFromFBConfig()`` returns NULL
for a config the other library allocated, and Qt6 calls ``qFatal`` — SIGABRT, no
traceback, nothing to debug. See ``docs/pyqt6-migration.md``.

This module lives apart from ``app.py`` for two reasons: the check has to run
before ``import napari`` (that import is what aborts), and keeping it here means
it can be tested without importing napari at all.

An install with no glvnd of its own cannot hit this: a wheel-only install
(``pip install palms``) brings Qt inside the PyQt6 wheel and pulls no conda
``libglx``, so PyOpenGL and Qt load the one host copy and there is no second
build to disagree with. The check requires the environment's own
``libGLX.so.0`` for that reason.

**It reports; it does not repair.** Preloading conda's ``libGLX.so.0`` with
``RTLD_GLOBAL`` does not help — PyOpenGL still ``dlopen``s the host's unversioned
copy as a separate mapping. Only the unversioned name existing *inside* the env
fixes it, which is what ``environment-linux.yml`` installs.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

#: Where a distribution's ``libglx-dev`` puts the unversioned link. Globbed, never
#: ``dlopen``ed — loading it is the bug this detects.
HOST_LIBGLX_PATHS = (
    "/usr/lib/x86_64-linux-gnu/libGLX.so",
    "/usr/lib/aarch64-linux-gnu/libGLX.so",
    "/usr/lib64/libGLX.so",
    "/usr/lib/libGLX.so",
)


def libglx_collision_message(
    prefix: str | os.PathLike[str] | None = None,
    host_paths: "tuple[str, ...] | None" = None,
    platform: str | None = None,
) -> str | None:
    """Return a warning for the defective configuration, or ``None``.

    **Three** things must hold, because the abort needs *two* glvnd builds in one
    process:

    1. the environment ships its own ``libGLX.so.0`` — that is the second build;
    2. it lacks the unversioned ``libGLX.so`` link, so PyOpenGL misses it;
    3. the host has an unversioned copy for PyOpenGL to find instead.

    Condition 1 was missing, and its absence made the warning fire on an install
    that cannot have the problem. Measured 2026-08-30 on one box: the conda dev
    env has both names and was correctly silent; a fresh ``pip install
    palms[cnv]`` env has **neither** — PyQt6 comes from a wheel with Qt bundled
    and nothing pulls conda's ``libglx`` — and the warning fired anyway, at a
    working app. With one glvnd in the process there is nothing to collide.

    Both file checks are for the names they name. Every box with working graphics
    has *a* ``libGLX.so.0``, and ``ctypes.util.find_library('GLX')`` returns
    exactly that, so it cannot answer this question — what matters is whether the
    *environment's* lib directory holds one, and whether the unversioned link
    sits beside it.

    The arguments exist for testing; the defaults describe the running process.
    """
    if not (platform or sys.platform).startswith("linux"):
        return None                     # no GLX on macOS; Qt uses the cocoa plugin

    prefix = Path(prefix if prefix is not None else sys.prefix)
    lib = prefix / "lib"
    if (lib / "libGLX.so").exists():
        return None                     # libglx-devel is installed — nothing to warn about
    if not (lib / "libGLX.so.0").exists():
        return None                     # no glvnd of its own — a pip/venv install, one copy

    for candidate in (host_paths if host_paths is not None else HOST_LIBGLX_PATHS):
        if os.path.exists(candidate):
            host = candidate
            break
    else:
        return None                     # no host copy to collide with

    return (
        "\n"
        "WARNING: this environment has libGLX.so.0 but not the unversioned libGLX.so\n"
        f"         (conda's `libglx-devel`), and a system libGLX was found at\n"
        f"         {host}. napari may abort at startup with\n"
        '         "Could not initialize GLX" / "Aborted (core dumped)" and no traceback,\n'
        "         most likely over a remote display (ThinLinc/VNC/x2go/xrdp).\n"
        "\n"
        "         Fix it with:\n"
        f"{_remedy(prefix)}"
    )


def _remedy(prefix: Path) -> str:
    """The command to run, for the environment actually running.

    The old text named ``-n palms`` and ``./scripts/install.sh`` unconditionally.
    An env is routinely called something else, and a user who installed from PyPI
    has no checkout to run a script from — advice that cannot be followed is
    worse than none, because it reads as though the app is broken.
    """
    env = prefix.name
    if not (prefix / "conda-meta").is_dir():
        return (
            f"             ln -s libGLX.so.0 {prefix / 'lib' / 'libGLX.so'}\n"
            "         (this environment provides libGLX.so.0 without the unversioned\n"
            "          link PyOpenGL looks for first)\n"
        )
    return (
        f"             mamba install -n {env} libglx-devel\n"
        "         or, from a checkout of the repository:\n"
        f"             mamba env update -n {env} -f environment-linux.yml\n"
        "         (./scripts/install.sh applies that automatically)\n"
    )


def warn_if_libglx_will_collide() -> bool:
    """Print the warning if it applies. Returns whether one was printed."""
    message = libglx_collision_message()
    if message is None:
        return False
    print(message, file=sys.stderr)
    return True
