"""Tab: Crop Dataset.

Lets the user draw one or more polygons in the "Crop Regions" napari Shapes
layer and export each as its own standalone, independently-openable
xenium-viewer data directory (cropped morphology image, cell/nucleus labels,
transcripts, and AnnData table). See utils/crop_export.py for the crop logic.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import PushButton
from qtpy.QtWidgets import QLabel, QFileDialog, QInputDialog, QLineEdit, QMessageBox
from qtpy.QtCore import QThread, Signal as QtSignal

from xenium_viewer.tabs._helpers import make_tab, StatusProxy, make_progress_dialog

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


class _CropExportWorker(QThread):
    """Runs crop_and_export for each job in sequence on a background thread."""

    progress = QtSignal(int, str)
    job_done = QtSignal(int, bool, str)   # job index, success, result-path-or-error
    finished_all = QtSignal()

    def __init__(self, ctx, jobs):
        super().__init__()
        self.ctx = ctx
        self.jobs = jobs   # list of (polygon_yx, output_dir, name)

    def run(self):
        from xenium_viewer.utils.crop_export import crop_and_export, CropExportError

        n = len(self.jobs)
        for i, (polygon_yx, output_dir, name) in enumerate(self.jobs):
            base_pct = int(i * 100 / n)

            def _cb(pct, msg, _i=i, _base=base_pct, _n=n):
                overall = _base + int(pct / _n)
                self.progress.emit(overall, f"[{_i + 1}/{_n}] {msg}")

            try:
                result_path = crop_and_export(
                    self.ctx, polygon_yx, output_dir, name, progress_cb=_cb,
                )
                self.job_done.emit(i, True, str(result_path))
            except CropExportError as exc:
                self.job_done.emit(i, False, str(exc))
            except Exception as exc:
                import traceback
                traceback.print_exc()
                self.job_done.emit(i, False, f"Unexpected error: {exc}")
        self.finished_all.emit()


def build_tab(ctx: ViewerContext) -> tuple:
    status = StatusProxy(ctx.viewer)

    instructions = QLabel(
        "Draw one or more polygons in the \"Crop Regions\" layer, then click "
        "\"Crop & Export\". You'll be asked to choose an output folder and a "
        "dataset name for each drawn region, in order. Each region is exported "
        "as its own standalone dataset (image + cell/nucleus labels + "
        "transcripts + table) that can be opened directly with xenium-viewer."
    )
    instructions.setWordWrap(True)

    draw_btn = PushButton(label="Activate Draw Polygon Tool")
    clear_btn = PushButton(label="Clear All Regions")
    export_btn = PushButton(label="Crop && Export")   # "&&" escapes to a literal "&" (Qt mnemonic syntax)

    def _on_draw():
        if ctx.crop_layer is None:
            return
        ctx.viewer.layers.selection.active = ctx.crop_layer
        ctx.crop_layer.mode = "add_polygon"

    def _on_clear():
        if ctx.crop_layer is None:
            return
        ctx.crop_layer.data = []
        status.value = "Crop regions cleared."

    def _prompt_destination(i: int, n: int):
        """Prompt for an output folder + name for region i. Returns (folder, name) or None if skipped/cancelled."""
        folder_str = QFileDialog.getExistingDirectory(
            None, f"Output folder for region {i + 1}/{n}",
        )
        if not folder_str:
            return None
        folder = Path(folder_str)

        default_name = f"crop_{i + 1}"
        while True:
            name_str, ok = QInputDialog.getText(
                None, f"Name for region {i + 1}/{n}",
                "Dataset folder name:", QLineEdit.Normal, default_name,
            )
            if not ok:
                return None
            name_str = name_str.strip()
            if not name_str or any(c in name_str for c in "/\\"):
                QMessageBox.warning(
                    None, "Invalid name",
                    "Enter a non-empty name with no path separators.",
                )
                continue
            if (folder / name_str).exists():
                reply = QMessageBox.question(
                    None, "Overwrite?",
                    f"{folder / name_str} already exists. Overwrite it?",
                    QMessageBox.Yes | QMessageBox.No,
                )
                if reply != QMessageBox.Yes:
                    continue
            return folder, name_str

    def _on_export():
        if ctx.crop_layer is None or len(ctx.crop_layer.data) == 0:
            status.value = "No crop regions drawn."
            return

        polygons = [np.asarray(p, dtype=np.float64) for p in ctx.crop_layer.data]
        n = len(polygons)
        jobs = []
        used_dests = set()

        for i, poly in enumerate(polygons):
            if len(poly) < 3:
                reply = QMessageBox.warning(
                    None, "Invalid region",
                    f"Region {i + 1}/{n} has fewer than 3 points and will be skipped. Continue?",
                    QMessageBox.Ok | QMessageBox.Cancel,
                )
                if reply == QMessageBox.Cancel:
                    return
                continue

            dest = _prompt_destination(i, n)
            if dest is None:
                reply = QMessageBox.question(
                    None, "Skip region?",
                    f"No destination chosen for region {i + 1}/{n}. Skip it and continue with the rest?",
                    QMessageBox.Yes | QMessageBox.Cancel,
                )
                if reply == QMessageBox.Cancel:
                    return
                continue

            folder, name = dest
            if (folder, name) in used_dests:
                QMessageBox.warning(
                    None, "Duplicate destination",
                    f"{folder / name} was already chosen for another region in this batch. "
                    f"Region {i + 1}/{n} will be skipped.",
                )
                continue
            used_dests.add((folder, name))
            jobs.append((poly, folder, name))

        if not jobs:
            status.value = "No regions to export."
            return

        dlg, bar, lbl = make_progress_dialog("Cropping Dataset(s)")
        worker = _CropExportWorker(ctx, jobs)
        results = []

        def _on_progress(pct, msg):
            bar.setValue(pct)
            lbl.setText(msg)

        def _on_job_done(i, success, info):
            results.append((jobs[i][2], success, info))

        def _on_finished_all():
            dlg.accept()
            lines = []
            n_ok = 0
            for name, success, info in results:
                if success:
                    n_ok += 1
                    lines.append(f"✓ {name} -> {info}")
                    ctx.record_code(
                        f"\n# Crop & Export dataset '{name}'\n"
                        f"# from xenium_viewer.utils.crop_export import crop_and_export\n"
                        f"# crop_and_export(ctx, polygon_yx=<region drawn in the \"Crop Regions\" layer>,\n"
                        f"#                 output_dir=Path({str(Path(info).parent)!r}), name={name!r})\n"
                        f"# -> wrote standalone dataset to \"{info}\""
                    )
                else:
                    lines.append(f"✗ {name}: {info}")
            status.value = f"Crop & Export finished: {n_ok}/{len(results)} succeeded."
            QMessageBox.information(None, "Crop & Export Results", "\n".join(lines))

        worker.progress.connect(_on_progress)
        worker.job_done.connect(_on_job_done)
        worker.finished_all.connect(_on_finished_all)
        worker.start()
        dlg.exec_()

    draw_btn.clicked.connect(_on_draw)
    clear_btn.clicked.connect(_on_clear)
    export_btn.clicked.connect(_on_export)

    widget = make_tab(instructions, draw_btn, clear_btn, export_btn)
    return widget, {"restore_session": lambda session: None}
