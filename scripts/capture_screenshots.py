#!/usr/bin/env python
"""
Automated screenshot capture for the PALMS wiki documentation.

Run with:
    conda run -n palms python scripts/capture_screenshots.py /path/to/xenium/output
    PALMS_SCREENSHOT_DATASET=/path/to/xenium/output \
        conda run -n palms python scripts/capture_screenshots.py

The dataset is an argument rather than a constant because these screenshots are
published: two of the panels (Tools > Dataset, Tools > Cache) print the dataset
path into the widget, so whatever is passed here ends up legible in docs/ and on
the wiki. Point it at a dataset whose path you are willing to publish.

Every tutorial shot **drives the app into the state its step describes** before
grabbing. Switching tabs is not enough: consecutive steps of a tutorial name the
same tab, so a capture that only navigated wrote byte-identical files — 30
tutorial files holding 9 distinct pictures, each tutorial illustrating every step
with a photograph of its first one.

The order of SHOTS is therefore load-bearing. Later shots depend on state earlier
ones established (a colouring, a computed registration, a figure in the Plots
dock), and a few actions are deliberately taken late because they overwrite what
an earlier shot needed — the H&E coarse alignment replaces the fine registration
restored from the session, so it is captured after the registration shots.

Steps whose only action is a file dialog (Save Landmarks…, Export GeoJSON…,
Save Volcano Plot…) are not captured: a modal dialog cannot be driven from here,
and a picture of the tab with nothing happening would be another duplicate.
Those tutorial pages carry one fewer image instead.

Saves PNGs to docs/screenshots/ and exits when done.
"""
import os
import sys
import time
from pathlib import Path

import numpy as np
from qtpy.QtCore import QPoint
from qtpy.QtGui import QImage, QPainter
from qtpy.QtWidgets import QApplication, QLineEdit, QPushButton, QTreeWidget, QWidget

SCREENSHOTS = Path(__file__).parent.parent / "docs/screenshots"


def _only_filter() -> str | None:
    """--only <substring>: re-take just the shots whose filename matches.

    A full run drives real analyses and takes half an hour; when one picture
    needs another attempt, redoing the other fifty-one is not the way to get it.
    Note the shots are ordered, and later ones depend on state earlier ones set,
    so a filtered run captures the named shot from a *freshly loaded* session.
    """
    if "--only" in sys.argv:
        i = sys.argv.index("--only")
        return sys.argv[i + 1] if i + 1 < len(sys.argv) else None
    return None


def _dataset_from_args() -> Path:
    """Dataset path from argv[1], else $PALMS_SCREENSHOT_DATASET."""
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--only" in sys.argv:
        i = sys.argv.index("--only")
        if i + 1 < len(sys.argv) and sys.argv[i + 1] in args:
            args.remove(sys.argv[i + 1])
    raw = args[0] if args else os.environ.get("PALMS_SCREENSHOT_DATASET")
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

# Outer tabs:  0=Cells  1=Genes  2=Spatial  3=Images  4=Tools
CELLS, GENES, SPATIAL, IMAGES, TOOLS = 0, 1, 2, 3, 4

# Everything a working session restores that is not the dataset itself is hidden
# by default: registration landmarks, an ARMS scan, patch overlays and prediction
# rasters sit outside the Xenium extent or on top of it, and none of them is what
# a reference image of the viewer should be showing. An allow-list rather than a
# block-list, because a session can restore arbitrarily named overlays and an
# unrecognised one must not end up in a published image.
BASE_LAYERS = ("cell_labels", "morphology_focus")
FRAME_LAYER = "cell_labels"


def _process_events(pause=0.08):
    QApplication.processEvents()
    time.sleep(pause)
    QApplication.processEvents()


class Rig:
    """The handle a shot's setup gets: the app, plus the few verbs it needs."""

    def __init__(self, viewer, ctx, dock, panel):
        self.viewer = viewer
        self.ctx = ctx
        self.dock = dock
        self.panel = panel

    # ── navigation ──────────────────────────────────────────────────────
    def navigate(self, outer, inner):
        self.panel.setCurrentIndex(outer)
        _process_events(0.05)
        sub = self.panel.currentWidget()
        if sub is not None and hasattr(sub, "setCurrentIndex"):
            sub.setCurrentIndex(inner)
        _process_events(0.08)

    @property
    def page(self):
        sub = self.panel.currentWidget()
        return sub.currentWidget() if hasattr(sub, "currentWidget") else sub

    # ── widgets ─────────────────────────────────────────────────────────
    @staticmethod
    def _norm(text):
        return "".join(ch for ch in str(text).lower() if ch.isalnum())

    def mg(self, label):
        """A magicgui widget on the current page, by its label.

        magicgui stores a back-reference on the Qt widget it wraps
        (``native._magic_widget``), so the whole tab can be driven through the
        real widget objects — ``.value = x``, ``.native.click()`` — without the
        app having to export them and without matching on Qt layout structure.
        """
        want = self._norm(label)
        for w in self.page.findChildren(QWidget):
            mw = getattr(w, "_magic_widget", None)
            if mw is not None and self._norm(getattr(mw, "label", "")) == want:
                return mw
        raise LookupError(f"no magicgui widget labelled {label!r} on this page")

    def qbtn(self, text):
        """A plain QPushButton on the current page, by its text."""
        want = self._norm(text)
        for b in self.page.findChildren(QPushButton):
            if self._norm(b.text()) == want:
                return b
        raise LookupError(f"no QPushButton {text!r} on this page")

    def click(self, label):
        try:
            self.mg(label).native.click()
        except LookupError:
            self.qbtn(label).click()
        _process_events(0.2)

    def wait_idle(self, label, timeout=600):
        """Wait for a button that disables itself while its worker runs."""
        widget = None
        try:
            widget = self.mg(label).native
        except LookupError:
            widget = self.qbtn(label)
        deadline = time.time() + timeout
        # give the worker a moment to actually disable it
        for _ in range(10):
            _process_events(0.1)
            if not widget.isEnabled():
                break
        while time.time() < deadline:
            _process_events(0.25)
            if widget.isEnabled():
                return True
        print(f"    ! timed out after {timeout}s waiting for {label!r}")
        return False

    def wait_for(self, cond, timeout=600, what=""):
        deadline = time.time() + timeout
        while time.time() < deadline:
            _process_events(0.25)
            if cond():
                return True
        print(f"    ! timed out after {timeout}s waiting for {what}")
        return False

    # ── layers and camera ───────────────────────────────────────────────
    def show_only(self, *prefixes):
        prefixes = prefixes or BASE_LAYERS
        for layer in self.viewer.layers:
            layer.visible = layer.name.startswith(tuple(prefixes))

    def show_also(self, *prefixes):
        for layer in self.viewer.layers:
            if layer.name.startswith(tuple(prefixes)):
                layer.visible = True

    def select(self, name):
        if name in self.viewer.layers:
            self.viewer.layers.selection = {self.viewer.layers[name]}

    def frame(self, name=FRAME_LAYER, margin=0.95, zoom_factor=1.0):
        """Put the camera on one layer's own extent.

        Not viewer.reset_view(): napari's fit_to_view measures
        layers._extent_world_augmented, which ignores `visible`, so a hidden
        ARMS scan far outside the tissue would still set the frame.
        """
        if name not in self.viewer.layers:
            return
        lo, hi = self.viewer.layers[name].extent.world
        size = np.maximum(hi - lo, 1.0)
        self.viewer.camera.center = (0.0, *((lo + hi) / 2.0))
        base = float(np.min(np.array(self.viewer._canvas_size) / size))
        self.viewer.camera.zoom = margin * base * zoom_factor

    def zoom_in(self, factor=4.0, name=FRAME_LAYER, offset=(0.0, 0.0)):
        self.frame(name)
        cy, cx = self.viewer.camera.center[1:]
        lo, hi = self.viewer.layers[name].extent.world
        span = hi - lo
        self.viewer.camera.center = (0.0, cy + offset[0] * span[0], cx + offset[1] * span[1])
        self.viewer.camera.zoom *= factor


# ── the shots ────────────────────────────────────────────────────────────────
# (outer, inner, filename, setup) — setup may return a viewer to grab instead of
# the main window ("the UMAP window is the picture the step is about").


def _launch(rig):
    rig.show_only()
    rig.select(FRAME_LAYER)
    rig.frame()


def _navigate_canvas(rig):
    rig.zoom_in(factor=6.0, offset=(-0.1, -0.15))


def _colour_by_gene(rig):
    rig.show_only()
    rig.frame()
    rig.mg("Color cells by").value = "Gene Expression"
    rig.mg("Gene").value = rig.ctx._shot_gene
    rig.mg("Colormap").value = "viridis"
    rig.click("Apply Cell Coloring")
    rig.wait_idle("Apply Cell Coloring", timeout=180)


def _colour_by_cluster(rig):
    rig.show_only()
    rig.frame()
    rig.mg("Color cells by").value = "Cluster"
    _process_events(0.2)
    rig.mg("Clustering").value = rig.ctx._shot_clustering
    rig.click("Apply Cell Coloring")
    rig.wait_idle("Apply Cell Coloring", timeout=180)


def _transcripts(rig):
    rig.mg("Transcript gene").value = rig.ctx._shot_gene
    rig.click("Add Gene")
    rig.mg("Show transcripts").value = True
    rig.click("Apply Transcripts")
    rig.wait_idle("Apply Transcripts", timeout=300)
    # Only the cells and the points: at this zoom the morphology channels fill
    # the canvas and the transcript spots are lost in them.
    # Only the cells and the points: at tissue zoom the morphology channels fill
    # the canvas and the spots are lost in them. The layer's default size is 4
    # data pixels — 0.85 µm, one screen pixel at anything short of cell-level
    # zoom — so the shot goes in close and widens the points, which is what
    # looking at transcripts in the app actually involves.
    rig.show_only("cell_labels", "transcripts")
    points = rig.viewer.layers["transcripts"]
    points.size = 10
    print(f"    transcripts: {len(points.data)} points, visible={points.visible}")
    # Centre on the median spot rather than a fraction of the tissue box: at this
    # zoom a fixed offset can land in a gland lumen, and a picture of transcripts
    # with no cells under them is not the step either.
    rig.frame()
    med = np.median(np.asarray(points.data), axis=0) * np.asarray(points.scale[-2:], dtype=float)
    rig.viewer.camera.center = (0.0, float(med[0]), float(med[1]))
    rig.viewer.camera.zoom *= 20.0


def _umap_window(rig):
    rig.click("Show UMAP Window")
    rig.wait_for(lambda: rig.ctx.umap_viewer._viewer is not None, 120, "the UMAP window")
    _process_events(1.0)
    return rig.ctx.umap_viewer._viewer


def _umap_by_cluster(rig):
    # The linked window mirrors whatever the cells are coloured by, and the shot
    # before this one coloured them by cluster — which is the whole point of the
    # step, and why this one is captured after it.
    rig.click("Show UMAP Window")
    _process_events(0.8)
    return rig.ctx.umap_viewer._viewer


def _roi_gene_colour(rig):
    rig.navigate(CELLS, 1)
    _colour_by_gene(rig)
    rig.navigate(SPATIAL, 0)
    rig.show_only(*BASE_LAYERS, "ROIs")
    rig.frame()


def _run_leiden(rig):
    rig.show_only()
    rig.frame()
    rig.mg("Resolution").value = 1.0
    rig.click("Run Leiden Clustering")
    # A step whose picture is "this is running" is still the picture of the step,
    # so a timeout here is not a failure — it is captured mid-run.
    rig.wait_idle("Run Leiden Clustering", timeout=900)


def _rank_genes(rig):
    rig.show_only()
    rig.frame()
    rig.mg("Clustering").value = rig.ctx._shot_clustering
    rig.mg("Method").value = "wilcoxon"
    rig.click("Run Rank Genes")
    rig.wait_idle("Run Rank Genes", timeout=900)


def _plot_umap_figure(rig):
    rig.click("Plot UMAP by cluster")
    rig.wait_idle("Plot UMAP by cluster", timeout=300)
    _process_events(1.0)


def _rois_drawn(rig):
    rig.show_only(*BASE_LAYERS, "ROIs")
    # The ROI layer's default 2-unit edge is 0.4 µm — sub-pixel at tissue scale,
    # so the regions the step is about would not be in the picture of them. 20
    # is what the ARMS tile layer already uses for the same reason.
    rois = rig.viewer.layers["ROIs"]
    rois.edge_width = [20.0] * len(rois.data)   # per shape: the scalar form is
    rois.current_edge_width = 20.0              # only the *next* shape's width
    rig.select("ROIs")
    rig.frame("ROIs", zoom_factor=0.75)


def _roi_expression(rig):
    rig.click("Calculate Expression")
    rig.wait_idle("Calculate Expression", timeout=600)


def _roi_deg(rig):
    rig.click("Run ROI DEG")
    rig.wait_idle("Run ROI DEG", timeout=900)


def _he_loaded(rig):
    rig.show_only(*BASE_LAYERS, "H&E (")
    rig.frame()
    rig.mg("Opacity").value = 70


def _he_landmarks(rig):
    rig.show_only(*BASE_LAYERS, "H&E (", "Xenium Landmarks", "H&E Landmarks")
    rig.select("Xenium Landmarks")
    rig.frame()


def _he_registration(rig):
    rig.show_only(*BASE_LAYERS, "H&E (", "Xenium Landmarks", "H&E Landmarks")
    try:
        button = rig.mg("Compute Registration").native
    except LookupError:
        button = None
    if button is not None and button.isEnabled():
        rig.click("Compute Registration")
        rig.wait_idle("Compute Registration", timeout=300)
    else:
        print("    note: Compute Registration is disabled (landmarks came from the "
              "session, not from this run) — showing the registered overlay instead")
    rig.mg("Opacity").value = 50
    rig.zoom_in(factor=2.5)


def _he_coarse_align(rig):
    # Last of the H&E shots on purpose: this replaces the fine registration the
    # session restored, which the two shots above are about.
    rig.show_only(*BASE_LAYERS, "H&E (")
    rig.frame()
    try:
        button = rig.mg("Coarse Align").native
    except LookupError:
        button = None
    if button is not None and button.isEnabled():
        rig.click("Coarse Align")
        rig.wait_idle("Coarse Align", timeout=300)
    else:
        print("    note: Coarse Align is disabled — showing the H&E over the cells")
    rig.mg("Opacity").value = 70


def _arms_loaded(rig):
    rig.show_only(*BASE_LAYERS, "ARMS H&E (")
    arms = next((lyr.name for lyr in rig.viewer.layers
                 if lyr.name.startswith("ARMS H&E (")), None)
    rig.frame(arms or FRAME_LAYER, zoom_factor=0.9)
    if arms:
        rig.select(arms)


def _arms_landmarks(rig):
    rig.show_only(*BASE_LAYERS, "ARMS H&E (", "ARMS Xenium Landmarks", "ARMS H&E Landmarks")
    rig.select("ARMS Xenium Landmarks")
    rig.frame()


def _arms_tiles(rig):
    rig.show_only(*BASE_LAYERS, "ARMS Tiles")
    rig.select("ARMS Tiles")
    rig.frame("ARMS Tiles", zoom_factor=0.8)


def _arms_deg(rig):
    rig.click("Run ARMS Tile DEG")
    rig.wait_idle("Run ARMS Tile DEG", timeout=900)


def _annotations_drawn(rig):
    layer = rig.viewer.layers["Annotations"]
    lo, hi = rig.viewer.layers[FRAME_LAYER].extent.world
    span = hi - lo
    # A shapes layer holds DATA coordinates and carries the pixel size in
    # layer.scale, so world µm handed to add_polygons come out 1/0.2125 too
    # small and pinned to the origin.
    scale = np.asarray(layer.scale[-2:], dtype=float)

    def box(fy, fx, fh, fw):
        y0, x0 = lo[0] + fy * span[0], lo[1] + fx * span[1]
        y1, x1 = y0 + fh * span[0], x0 + fw * span[1]
        return np.array([[y0, x0], [y0, x1], [y1, x1], [y1, x0]]) / scale

    if len(layer.data) == 0:
        layer.add_polygons([box(0.20, 0.08, 0.34, 0.20),
                            box(0.28, 0.42, 0.38, 0.24),
                            box(0.52, 0.72, 0.30, 0.20)])
    layer.edge_width = 20
    rig.show_only(*BASE_LAYERS, "Annotations")
    rig.select("Annotations")
    rig.frame()


def _annotate(rig, indices, name):
    layer = rig.viewer.layers["Annotations"]
    layer.selected_data = set(indices)
    _process_events(0.2)
    edits = rig.page.findChildren(QLineEdit)
    if edits:
        edits[0].setText(name)
    rig.click("Assign to selected shapes")
    _process_events(0.4)


def _annotations_first_type(rig):
    _annotate(rig, [0], "Tumour")


def _annotations_second_type(rig):
    _annotate(rig, [1, 2], "Stroma")


def _overview(rig):
    """The hero shot, taken last: by now a figure exists, so the Plots dock —
    which Interface-Overview.md describes — is part of the interface."""
    rig.show_only()
    rig.select(FRAME_LAYER)
    rig.frame()


WINDOW_SHOTS = [
    # Getting Started
    (CELLS, 0, "tutorial-getting-started-step1.png", _launch),
    (CELLS, 0, "tutorial-getting-started-step3.png", _navigate_canvas),
    (CELLS, 1, "tutorial-getting-started-step5.png", _colour_by_gene),
    (CELLS, 2, "tutorial-getting-started-step6.png", _transcripts),
    (CELLS, 3, "tutorial-getting-started-step7.png", _umap_window),
    # Clustering
    (CELLS, 0, "tutorial-clustering-step1.png", _run_leiden),
    (CELLS, 1, "tutorial-clustering-step2.png", _colour_by_cluster),
    (CELLS, 3, "tutorial-clustering-step3.png", _umap_by_cluster),
    (GENES, 0, "tutorial-clustering-step4.png", _rank_genes),
    (CELLS, 3, "tutorial-clustering-step8.png", _plot_umap_figure),
    # ROI analysis
    (SPATIAL, 0, "tutorial-roi-analysis-step1.png", _roi_gene_colour),
    (SPATIAL, 0, "tutorial-roi-analysis-step3.png", _rois_drawn),
    (SPATIAL, 0, "tutorial-roi-analysis-step4.png", _roi_expression),
    (SPATIAL, 0, "tutorial-roi-analysis-step6.png", _roi_deg),
    # H&E registration
    (IMAGES, 0, "tutorial-he-registration-step2.png", _he_loaded),
    (IMAGES, 0, "tutorial-he-registration-step5.png", _he_landmarks),
    (IMAGES, 0, "tutorial-he-registration-step6.png", _he_registration),
    (IMAGES, 0, "tutorial-he-registration-step4.png", _he_coarse_align),
    # ARMS overlay
    (IMAGES, 1, "tutorial-arms-overlay-step2.png", _arms_loaded),
    (IMAGES, 1, "tutorial-arms-overlay-step4.png", _arms_landmarks),
    (IMAGES, 1, "tutorial-arms-overlay-step7.png", _arms_tiles),
    (IMAGES, 1, "tutorial-arms-overlay-step11.png", _arms_deg),
    # Annotations
    (TOOLS, 0, "tutorial-annotations-step2.png", _annotations_drawn),
    (TOOLS, 0, "tutorial-annotations-step4.png", _annotations_first_type),
    (TOOLS, 0, "tutorial-annotations-step6.png", _annotations_second_type),
    # Hero, last
    (CELLS, 0, "interface-overview.png", _overview),
]


def _scan_dataset(rig):
    rig.click("Scan Dataset")
    rig.wait_for(lambda: any(t.topLevelItemCount() > 0
                             for t in rig.page.findChildren(QTreeWidget)),
                 600, "the dataset scan")
    _process_events(0.5)


def _pick_template(rig):
    trees = rig.page.findChildren(QTreeWidget)
    if not trees:
        return
    tree = trees[0]
    for i in range(tree.topLevelItemCount()):
        group = tree.topLevelItem(i)
        if group.childCount():
            tree.setCurrentItem(group.child(0))
            tree.expandItem(group)
            break
    _process_events(0.6)


def _check_celldega(rig):
    """Ask the Publish tab whether the optional dependency is there.

    Its readout is an empty black box until something has been asked of it, and
    what it says either way is the useful half of the picture: the version when
    celldega is installed, and the exact ``--no-deps`` install line when it is
    not — which is the state most readers of this page will be in.
    """
    rig.click("Check Celldega")
    _process_events(0.6)


# Tab reference screenshots: grab the dock widget only. Taken after the tutorial
# shots so the panels show results rather than empty boxes.
#
# Kept as plain (outer, inner, filename) literals, with the setups in the map
# below: the indices are positional into app.py's addTab order, and
# tests/test_docs_links.py reads this list with ast.literal_eval to check them
# against that order — a guard that exists because tab-notebook.png was a
# picture of the Crop Dataset tab for as long as the list was one entry short.
TAB_SHOTS = [   # literal ints, not the names above: literal_eval reads this

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
    # QC was inserted at Tools index 1 in 2026-09, beside Segmentation, and
    # every entry below shifted by one.
    (4, 1, "tab-qc.png"),
    # Preprocess was inserted at Tools index 2 in 2026-09, after QC so the tab
    # order reads filter -> normalise; everything below shifted by one again.
    (4, 2, "tab-preprocess.png"),
    (4, 3, "tab-segmentation.png"),
    # Tools index 4 is Crop Dataset, not Notebook. This entry said
    # "tab-notebook.png", so every capture run overwrote the Notebook screenshot
    # with a picture of the Crop Dataset tab. Indices here are positional into
    # app.py's addTab() order — when a tab is inserted, this list shifts with it.
    (4, 4, "tab-crop-dataset.png"),
    (4, 5, "tab-publish.png"),
    (4, 6, "tab-notebook.png"),
    (4, 7, "tab-dataset.png"),
    (4, 8, "tab-cache.png"),
    (4, 9, "tab-templates.png"),
]


# Three tabs photograph as empty boxes unless something has been asked of them.
TAB_SETUPS = {
    "tab-dataset.png": _scan_dataset,
    "tab-templates.png": _pick_template,
    "tab-publish.png": _check_celldega,
}


def _grab_dock(dock):
    _process_events()
    return dock.widget().grab()


def _grab_window(viewer):
    """Render a whole viewer window into a pixmap, with a freshly drawn canvas.

    Two things had to change here, and both were producing wrong images.

    Deliberately NOT QScreen.grabWindow(winId): under Qt6 that argument form is
    unsupported on several platforms and silently returns a fragment of the root
    window instead of the app. It produced 31 published screenshots of desktop
    wallpaper. QWidget.grab() renders the widget tree itself, so it is
    independent of the compositor and works over remote X. It also excludes the
    WM title bar, which is why viewer.title can never leak the dataset folder
    name into a published image.

    But a widget grab reads the vispy canvas's *last painted* framebuffer, and
    part of that framebuffer is not repainted when only the camera and layer
    visibility change: shots came back with a block of an earlier draw frozen in
    the corner, identical in every frame while the rest of the canvas tracked the
    layers correctly. Calling repaint() on the canvas widget does not clear it.
    So the canvas is rendered separately through vispy — screenshot(canvas_only=
    True), which draws on demand rather than reading what is on screen — and
    painted over the canvas widget's rectangle in the grab.
    """
    _process_events()
    qt_window = viewer.window._qt_window
    pixmap = qt_window.grab()
    native = viewer.window._qt_viewer.canvas.native
    arr = np.ascontiguousarray(viewer.screenshot(canvas_only=True, flash=False))
    height, width = arr.shape[:2]
    image = QImage(arr.data, width, height, 4 * width,
                   QImage.Format.Format_RGBA8888).copy()
    painter = QPainter(pixmap)
    painter.drawImage(native.mapTo(qt_window, QPoint(0, 0)),
                      image.scaled(native.width(), native.height()))
    painter.end()
    return pixmap


def _save(pixmap, filename):
    pixmap.save(str(SCREENSHOTS / filename))
    print(f"  {filename}")


def _pick_subjects(ctx):
    """One gene and one clustering, chosen from the data rather than hardcoded."""
    counts = np.asarray(ctx.adata.X.sum(axis=0)).ravel()
    order = np.argsort(counts)[::-1]
    gene = next((str(ctx.adata.var_names[i]) for i in order
                 if str(ctx.adata.var_names[i]) in set(ctx.gene_names)), ctx.gene_names[0])
    names = list(ctx.clustering_names)
    clustering = next((n for n in names if "leiden" in n), names[0] if names else None)
    ctx._shot_gene = gene
    ctx._shot_clustering = clustering
    print(f"  subjects: gene={gene!r} clustering={clustering!r}")


def capture_all(viewer, ctx, dock, panel):
    SCREENSHOTS.mkdir(parents=True, exist_ok=True)
    only = _only_filter()
    window_shots = [s for s in WINDOW_SHOTS if only is None or only in s[2]]
    tab_shots = [s for s in TAB_SHOTS if only is None or only in s[2]]
    if only:
        print(f"  --only {only!r}: {len(window_shots) + len(tab_shots)} shot(s)")
    rig = Rig(viewer, ctx, dock, panel)
    qt_window = viewer.window._qt_window
    qt_window.raise_()
    qt_window.activateWindow()
    _process_events(0.3)

    _pick_subjects(ctx)
    print("  layers:", ", ".join(lyr.name for lyr in viewer.layers))

    print("\nTutorial screenshots (each drives the app into its step's state):")
    for outer, inner, fname, setup in window_shots:
        rig.navigate(outer, inner)
        target = None
        if setup is not None:
            t0 = time.time()
            try:
                target = setup(rig)
            except Exception as exc:                       # noqa: BLE001
                print(f"    ! {fname}: setup failed: {type(exc).__name__}: {exc}")
            print(f"    ({time.time() - t0:.0f}s) {fname}")
        _process_events(0.4)
        _save(_grab_window(target or viewer), fname)

    print("\nTab reference screenshots:")
    for outer, inner, fname in tab_shots:
        rig.navigate(outer, inner)
        setup = TAB_SETUPS.get(fname)
        if setup is not None:
            try:
                setup(rig)
            except Exception as exc:                       # noqa: BLE001
                print(f"    ! {fname}: setup failed: {type(exc).__name__}: {exc}")
        _save(_grab_dock(dock), fname)

    count = len(window_shots) + len(tab_shots)
    print(f"\n{count} screenshots saved to {SCREENSHOTS}")

    # napari.run() returns when the last top-level window closes, and the linked
    # UMAP window is a second one — closing only the main viewer left the process
    # alive after the final grab, with every file already written.
    umap = getattr(ctx.umap_viewer, "_viewer", None)
    if umap is not None:
        try:
            umap.close()
        except Exception:                                  # noqa: BLE001,S110
            pass
    viewer.close()
    QApplication.quit()


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
    # Wider than the old 1400x900: at that size the Controls dock and the layer
    # panel left the canvas ~400px across, so the tissue was a thumbnail in the
    # middle of a reference image about a spatial viewer.
    viewer.window.resize(1800, 1000)

    print(f"Loading dataset: {dataset}")
    ctx = _do_full_init(viewer, dataset, no_cache=False, _app=_app)

    dock = _app["dock_widget"]
    # The dock wraps the outer QTabWidget (5 top-level groups)
    panel = dock.widget()

    # Give the window time to render before starting captures
    QTimer.singleShot(2500, lambda: capture_all(viewer, ctx, dock, panel))

    napari.run()
    print("Done.")


if __name__ == "__main__":
    main()
