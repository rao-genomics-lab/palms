# Migrating napari off the deprecated PyQt5 backend → PyQt6

Status: **done** (2026-08-24), issue #15. This file is now the record of what the
migration actually required, because the plan it replaced was wrong in both
directions — it overstated the code edits and understated the packaging.

## Why

napari deprecated the PyQt5 backend for removal in **autumn 2026**, and says so at
every startup:

> napari support for the PyQt5 backend is deprecated and will be removed in fall of 2026.

Under PyQt5 it also prints a second one, which is the more visible symptom:

> System theme detection requires a Qt6 backend. Please switch to PyQt6 or PySide6.

## What the previous plan got wrong

The earlier version of this document framed the work as "8 unscoped enums + 7
`.exec_()` calls, then flip the pin". Three corrections, each found by measurement:

**1. The enum inventory was understated ~12×.** The tree had **98** unscoped enum
sites, not 8: ~29 `Qt.*` plus ~69 on other classes — `QMessageBox.Yes` ×17,
`QMessageBox.Cancel` ×12, `QScrollArea.NoFrame` ×6, `QDialog.Accepted`,
`QHeaderView.Stretch`, `QImage.Format_RGBA8888`, and more. There were 8 `.exec_()`
sites in `src/`, not 7.

**2. None of it was actually a blocker.** qtpy shims both categories away under
PyQt6, so the old code would have run unchanged:

- `qtpy/enums_compat.py::promote_enums()` runs **only** under PyQt6 and is applied
  to `QtCore`, `QtGui`, `QtWidgets` and `QtTest`. It copies every scoped enum member
  back up to its class namespace, restoring `Qt.Horizontal` / `QMessageBox.Yes`.
  Every class this codebase touched unscoped lives in one of those four modules.
- qtpy aliases `exec_` → `exec` on `QDialog`, `QApplication`, `QMenu`,
  `QCoreApplication`, `QEventLoop` and `QThread`. All 8 sites were `QDialog` or
  `QMessageBox` (a `QDialog` subclass).

The edits were made anyway — they are valid under PyQt5 5.11+, they cost one
scripted pass, and relying on a compatibility shim for 98 sites is a liability
rather than a saving. But **the dependency solve was the real work**, and the plan
did not mention it at all.

**3. `QT_API=pyside6 pytest` does not smoke-test strict enums.** PySide6 6.10
runs in forgiveness mode: probed directly, `QMessageBox.Yes` → `16384`,
`QDialog.Accepted` → `1`, `Qt.Horizontal` → `Orientation.Horizontal`, and qtpy
supplies `exec_` on top. Only a real PyQt6 run tests this, which is why CI grew a
PyQt6 leg rather than a PySide6 one.

## Packaging: the part that took the time

- **conda-forge ships PyQt6 as `pyqt6`, not as a 6.x of `pyqt`.** The `pyqt`
  package stops at 5.15.11, so `pyqt=6` fails with *"does not exist (perhaps a
  typo or a missing channel)"*, which reads like PyQt6 is unavailable on
  conda-forge. It is not — the name is just different. (`pyside6` is packaged
  under its own name too, and is already pulled in as a matplotlib dependency.)
- **`napari` had to be pinned.** With `pyqt6` swapped in, the solver silently
  resolved `napari` to **0.7.0** instead of 0.8.0 — not a Qt conflict at all, but
  the solver easing the unrelated `zarr>=3.0,<3.2` pin by walking napari back a
  minor version. `napari=0.8.0` + `pyqt6` solves cleanly on its own. So
  `environment.yml` now says `napari>=0.8`, and a silent downgrade of the thing
  this application *is* can no longer happen.
- Resolved stack, verified in a real env: napari 0.8.0, PyQt6 6.8.1 / Qt 6.8.1,
  qtpy 2.4.3, magicgui 0.10.2, superqt 0.8.2, spatialdata 0.8.0. No PyQt5 present.

## How to reproduce the verification

```bash
mamba env create -f environment.yml     # now resolves the Qt6 stack
conda activate xenium_viewer
python -c "import qtpy; print(qtpy.API_NAME, qtpy.QT_VERSION)"   # PyQt6 6.8.1
pytest
```

CI (`.github/workflows/ci.yml`) does not add a second backend leg — the conda env
*is* the PyQt6 leg now. What it adds instead is an assertion, before the suite
runs, that the solved env really is PyQt6 and that napari did not walk backwards.
Both of those regressions are silent, both have happened once already, and a green
suite on the wrong backend is precisely what this migration was meant to end.
