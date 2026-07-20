# TODO: Migrate napari off the deprecated PyQt5 backend → PyQt6

Status: **deferred / tracked** (as of 2026-07-19). No code changes made yet.

## Why

Napari prints at startup:

> napari support for the PyQt5 backend is deprecated and will be removed in fall of 2026.

This is a deprecation notice, **not** a failure — the viewer runs fine on PyQt5 today
and will keep working until napari actually drops the backend (their stated target:
**fall 2026**). We are keeping PyQt5 for now and migrating deliberately before that
deadline, after exploring dependency consequences.

## Current state (why the migration is small)

The codebase is already well-positioned:

- **All Qt access goes through `qtpy`** — there are no direct `PyQt5` imports anywhere in
  `src/`. The only hard PyQt5 references are the backend *pins*:
  - `environment.yml` (`pyqt=5`)
  - `pyproject.toml` (`"PyQt5"`)
- `QAction` is already imported from the Qt6-correct location (`qtpy.QtGui`), signals use
  qtpy's `Signal`, and there are **no** `QRegExp`, `QFontMetrics.width`, high-DPI
  attribute, or `setMargin` landmines.
- Only two categories of code would break under PyQt6's strict enum mode (below).

The separate CopyKAT env (`environment-copykat.yml` → `xenium_viewer_copykat`) is a headless
R/rpy2 worker with no Qt GUI; it is unaffected by this migration.

## Migration checklist (execute later)

1. **Explore consequences first.** In a scratch env, confirm `napari>=0.5`, `magicgui`,
   `superqt`, and `spatialdata` all resolve cleanly against **PyQt6**. (PySide6 is the
   LGPL alternative if licensing ever matters; PyQt6 keeps the same family as today's
   PyQt5.)

2. **Fix 8 unscoped enums → scoped.** The scoped form also works under PyQt5 5.11+, so
   these edits are backend-agnostic and can be made and tested *before* flipping the pin:
   - `src/xenium_viewer/tabs/tab_patch_overlays.py:183,193,202` — `Qt.Horizontal` → `Qt.Orientation.Horizontal`
   - `src/xenium_viewer/tabs/tab_external_images.py:70,244` — `Qt.Horizontal` → `Qt.Orientation.Horizontal`
   - `src/xenium_viewer/utils/minimap_widget.py:54` — `Qt.WA_TranslucentBackground` → `Qt.WidgetAttribute.WA_TranslucentBackground`
   - `src/xenium_viewer/tabs/tab_notebook.py:205` — `Qt.AlignLeft` → `Qt.AlignmentFlag.AlignLeft`
   - `src/xenium_viewer/app.py:1173` — `_Qt.Vertical` → `_Qt.Orientation.Vertical`

3. **Fix 7 `.exec_()` → `.exec()`** (also valid under PyQt5):
   `loader.py:165`, `app.py:1217,1456`, `tab_crop_dataset.py:212`,
   `tab_patch_overlays.py:64`, `_helpers.py:547`, `adata_persistence.py:62`.
   - Note: `tab_segmentation.py:491` `spec.loader.exec_module(...)` is unrelated (importlib) — leave it.

4. **Flip the backend pin.** `environment.yml` `pyqt=5` → `pyqt` (Qt6 build);
   `pyproject.toml` `"PyQt5"` → `"PyQt6"`.

5. **Verify** the viewer launches and each tab that touches the edited sites works
   (patch overlays, external images, minimap, notebook figure label, console dock split).

6. **Update `CHANGELOG.md`** when the migration lands.

> Line numbers above are from 2026-07-19 and may drift; grep for the enum/`exec_` patterns
> if they no longer match.
