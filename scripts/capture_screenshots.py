#!/usr/bin/env python
"""
Automated screenshot capture for palms wiki documentation.

Run with:
    conda run -n palms python scripts/capture_screenshots.py /path/to/xenium/output
    PALMS_SCREENSHOT_DATASET=/path/to/xenium/output \
        conda run -n palms python scripts/capture_screenshots.py

The dataset is an argument rather than a constant because these screenshots are
published: two of the panels (Tools > Dataset, Tools > Cache) print the dataset
path into the widget, so whatever is passed here ends up legible in docs/ and on
the wiki. Point it at a dataset whose path you are willing to publish.

Saves PNGs to docs/screenshots/ and exits when done.
"""
import os
import sys
import time
from pathlib import Path

SCREENSHOTS = Path(__file__).parent.parent / "docs/screenshots"


def _dataset_from_args() -> Path:
    """Dataset path from argv[1], else $PALMS_SCREENSHOT_DATASET."""
    raw = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("PALMS_SCREENSHOT_DATASET")
    if not raw:
        sys.exit(
            "error: no dataset given.\n"
            "  usage: python scripts/capture_screenshots.py /path/to/xenium/output\n"
            "     or: PALMS_SCREENSHOT_DATASET=/path/to/xenium/output "
            "python scripts/capture_screenshots.py"
        )
    path = Path(raw).expanduser()
    if not path.is_dir():
        sys.exit(f"error: not a directory: {path}")
    return path

os.environ.setdefault("DISPLAY", ":0")

# ── Tab navigation map ───────────────────────────────────────────────────────
# Outer tabs:  0=Cells  1=Genes  2=Spatial  3=Images  4=Tools
# Inner tabs depend on group — see _build_control_panel() in app.py

# Tab reference screenshots: grab the dock widget only
TAB_SHOTS = [
    # (outer, inner, filename)
    (0, 0, "tab-clustering.png"),
    (0, 1, "tab-cell-coloring.png"),
    (0, 2, "tab-transcripts.png"),
    (0, 3, "tab-umap.png"),
    (1, 0, "tab-rank-genes.png"),
    (1, 1, "tab-markers.png"),
    (1, 2, "tab-gene-correlation.png"),
    (1, 3, "tab-cnv.png"),
    (2, 0, "tab-roi-analysis.png"),
    (2, 1, "tab-ligand-receptor.png"),
    (2, 2, "tab-neighborhood-enrichment.png"),
    (2, 3, "tab-co-occurrence.png"),
    (2, 4, "tab-domains.png"),
    (2, 5, "tab-annot-nhood.png"),
    (2, 6, "tab-annot-distance.png"),
    (3, 0, "tab-he-registration.png"),
    (3, 1, "tab-arms-overlay.png"),
    (3, 2, "tab-external-images.png"),
    (3, 3, "tab-patches.png"),
    (4, 0, "tab-annotations.png"),
    (4, 1, "tab-segmentation.png"),
    # Tools index 2 is Crop Dataset, not Notebook. This entry said
    # "tab-notebook.png", so every capture run overwrote the Notebook screenshot
    # with a picture of the Crop Dataset tab. Indices here are positional into
    # app.py's addTab() order — when a tab is inserted, this list shifts with it.
    (4, 2, "tab-crop-dataset.png"),
    (4, 3, "tab-notebook.png"),
    (4, 4, "tab-dataset.png"),
    (4, 5, "tab-cache.png"),
    (4, 6, "tab-templates.png"),
]

# Tutorial screenshots: full-window grab, navigate to the indicated tab first
TUTORIAL_SHOTS = [
    # Tutorial-Getting-Started
    (0, 0, "tutorial-getting-started-step1.png"),
    (0, 0, "tutorial-getting-started-step3.png"),
    (0, 1, "tutorial-getting-started-step5.png"),
    (0, 2, "tutorial-getting-started-step6.png"),
    (0, 3, "tutorial-getting-started-step7.png"),
    # Tutorial-Clustering
    (0, 0, "tutorial-clustering-step1.png"),
    (0, 1, "tutorial-clustering-step2.png"),
    (0, 3, "tutorial-clustering-step3.png"),
    (1, 0, "tutorial-clustering-step4.png"),
    (0, 3, "tutorial-clustering-step8.png"),
    # Tutorial-HE-Registration
    (3, 0, "tutorial-he-registration-step2.png"),
    (3, 0, "tutorial-he-registration-step4.png"),
    (3, 0, "tutorial-he-registration-step5.png"),
    (3, 0, "tutorial-he-registration-step6.png"),
    (3, 0, "tutorial-he-registration-step8.png"),
    # Tutorial-ROI-Analysis
    (2, 0, "tutorial-roi-analysis-step1.png"),
    (2, 0, "tutorial-roi-analysis-step3.png"),
    (2, 0, "tutorial-roi-analysis-step4.png"),
    (2, 0, "tutorial-roi-analysis-step6.png"),
    (2, 0, "tutorial-roi-analysis-step8.png"),
    # Tutorial-Annotations
    (4, 0, "tutorial-annotations-step2.png"),
    (4, 0, "tutorial-annotations-step4.png"),
    (4, 0, "tutorial-annotations-step6.png"),
    (4, 0, "tutorial-annotations-step7.png"),
    # Tutorial-ARMS-Overlay
    (3, 1, "tutorial-arms-overlay-step2.png"),
    (3, 1, "tutorial-arms-overlay-step4.png"),
    (3, 1, "tutorial-arms-overlay-step6.png"),
    (3, 1, "tutorial-arms-overlay-step7.png"),
    (3, 1, "tutorial-arms-overlay-step11.png"),
    (3, 1, "tutorial-arms-overlay-step12.png"),
]


def _process_events(pause=0.08):
    from qtpy.QtWidgets import QApplication
    QApplication.processEvents()
    time.sleep(pause)
    QApplication.processEvents()


def _navigate(panel, outer, inner):
    """Set outer and inner tab, then flush the event queue."""
    panel.setCurrentIndex(outer)
    _process_events(0.05)
    sub = panel.currentWidget()
    if sub is not None and hasattr(sub, "setCurrentIndex"):
        sub.setCurrentIndex(inner)
    _process_events(0.08)


def _grab_full_window(qt_window):
    from qtpy.QtWidgets import QApplication
    _process_events()
    return QApplication.primaryScreen().grabWindow(qt_window.winId())


def _grab_dock(dock):
    _process_events()
    return dock.widget().grab()


def _save(pixmap, filename):
    path = str(SCREENSHOTS / filename)
    pixmap.save(path)
    print(f"  {filename}")


def capture_all(viewer, dock, panel):
    SCREENSHOTS.mkdir(parents=True, exist_ok=True)
    qt_window = viewer.window._qt_window

    # Raise and activate the window so grabs aren't occluded
    qt_window.raise_()
    qt_window.activateWindow()
    _process_events(0.3)

    # ── Interface overview (full window, Cells > Clustering visible) ──────────
    _navigate(panel, 0, 0)
    _save(_grab_full_window(qt_window), "interface-overview.png")

    # ── Tab reference screenshots (dock panel only) ───────────────────────────
    print("\nTab reference screenshots:")
    for outer, inner, fname in TAB_SHOTS:
        _navigate(panel, outer, inner)
        _save(_grab_dock(dock), fname)

    # ── Tutorial screenshots (full window) ────────────────────────────────────
    print("\nTutorial screenshots:")
    for outer, inner, fname in TUTORIAL_SHOTS:
        _navigate(panel, outer, inner)
        _save(_grab_full_window(qt_window), fname)

    count = 1 + len(TAB_SHOTS) + len(TUTORIAL_SHOTS)
    print(f"\n{count} screenshots saved to {SCREENSHOTS}")
    viewer.close()


def main():
    # Resolve the dataset before the napari import: a missing argument should
    # fail in a second, not after a ten-second Qt/scanpy import.
    dataset = _dataset_from_args()

    sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

    import napari
    from qtpy.QtCore import QTimer
    # Import internal helpers — these are private but stable across the codebase
    from palms.app import _do_full_init  # noqa: PLC2701

    _app = {
        "dock_widget": None,
        "restore_fn": None,
        "snapshot": {},
        "reload_in_progress": False,
    }

    viewer = napari.Viewer(title="PALMS — Screenshot Capture")
    viewer.window.resize(1400, 900)

    print(f"Loading dataset: {dataset}")
    _do_full_init(viewer, dataset, no_cache=False, _app=_app)

    dock = _app["dock_widget"]
    # The dock wraps the outer QTabWidget (5 top-level groups)
    panel = dock.widget()

    # Give the window time to render before starting captures
    QTimer.singleShot(1500, lambda: capture_all(viewer, dock, panel))

    napari.run()
    print("Done.")


if __name__ == "__main__":
    main()
