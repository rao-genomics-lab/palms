"""PALMS — Provenance-Aware Linking of Multimodal Spatial-omics.

A napari-based viewer that registers spatial transcriptomics (Xenium 3.x),
histology and genomic overlays into one coordinate space, and records every
action as replayable code.
"""

import os
from importlib.util import find_spec

__version__ = "0.1.0"

# ── Qt backend selection ─────────────────────────────────────────────────────
# qtpy picks a binding in the order PyQt5, PySide2, PyQt6, PySide6, so *any*
# environment that still contains PyQt5 runs on it — silently, even when PyQt6 is
# installed and intended. That is not hypothetical: conda-forge's `matplotlib`
# metapackage pulls `pyqt` (Qt5) at some versions and `pyside6` at others, which
# is exactly how a CI run with pyqt6 in environment.yml still came up PyQt5.
#
# environment.yml now asks for matplotlib-base so nothing drags a second binding
# in, but a user upgrading an existing env keeps whatever is already there. So
# state the choice rather than inheriting qtpy's default:
#
#   * setdefault, so an explicit QT_API (a developer testing PySide6, say) wins;
#   * only when PyQt6 is actually importable, so an env that legitimately has
#     just PySide6 is left to qtpy rather than pointed at something absent.
#
# find_spec does not import the module, so this costs nothing at startup, and it
# runs here because `import palms` precedes every qtpy import we make.
if find_spec("PyQt6") is not None:
    os.environ.setdefault("QT_API", "pyqt6")
