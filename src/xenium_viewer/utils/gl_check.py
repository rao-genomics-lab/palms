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

    Both halves of the collision must be present for it to be reachable, so an
    env merely missing the package on a box with no host ``libGLX.so`` says
    nothing — there is no second copy for PyOpenGL to find.

    Both checks are for the *unversioned* name specifically. Every box with
    working graphics has ``libGLX.so.0``, and ``ctypes.util.find_library('GLX')``
    finds exactly that, so it cannot answer this question. It is the unversioned
    link that decides which copy PyOpenGL loads.

    The arguments exist for testing; the defaults describe the running process.
    """
    if not (platform or sys.platform).startswith("linux"):
        return None                     # no GLX on macOS; Qt uses the cocoa plugin

    prefix = Path(prefix if prefix is not None else sys.prefix)
    if (prefix / "lib" / "libGLX.so").exists():
        return None                     # libglx-devel is installed — nothing to warn about

    for candidate in (host_paths if host_paths is not None else HOST_LIBGLX_PATHS):
        if os.path.exists(candidate):
            host = candidate
            break
    else:
        return None                     # no host copy to collide with

    return (
        "\n"
        "WARNING: this environment is missing `libglx-devel`, and a system libGLX\n"
        f"         was found at {host}. napari may abort at startup with\n"
        '         "Could not initialize GLX" / "Aborted (core dumped)" and no traceback,\n'
        "         most likely over a remote display (ThinLinc/VNC/x2go/xrdp).\n"
        "\n"
        "         Fix it with:\n"
        "             mamba env update -n xenium_viewer -f environment-linux.yml\n"
        "         (or re-run ./scripts/install.sh, which applies it automatically)\n"
    )


def warn_if_libglx_will_collide() -> bool:
    """Print the warning if it applies. Returns whether one was printed."""
    message = libglx_collision_message()
    if message is None:
        return False
    print(message, file=sys.stderr)
    return True
