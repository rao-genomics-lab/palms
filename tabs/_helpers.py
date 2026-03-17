"""Shared helpers used across tab modules."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
from qtpy.QtWidgets import QWidget, QVBoxLayout

if TYPE_CHECKING:
    from utils.viewer_context import ViewerContext


# ── Tab layout helper ────────────────────────────────────────────────────────

def make_tab(*widgets_and_natives) -> QWidget:
    """Pack magicgui widgets and raw QWidgets into a single QWidget."""
    tab = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(4, 4, 4, 4)
    for w in widgets_and_natives:
        if hasattr(w, "native"):
            layout.addWidget(w.native)
        else:
            layout.addWidget(w)
    layout.addStretch()
    tab.setLayout(layout)
    return tab


# ── Status proxy ─────────────────────────────────────────────────────────────

class StatusProxy:
    """Redirect `.value = msg` to napari viewer.status."""

    def __init__(self, viewer):
        self._viewer = viewer

    @property
    def value(self):
        return self._viewer.status

    @value.setter
    def value(self, msg):
        self._viewer.status = msg


# ── Progress helpers ──────────────────────────────────────────────────────────

class ProgressMailbox:
    """Thread-safe single-slot message passing (background → main thread)."""
    def __init__(self):
        self._msg = None

    def post(self, msg: str):       # called from background thread
        self._msg = msg

    def read(self) -> str | None:   # called from main thread
        msg = self._msg
        self._msg = None
        return msg


def attach_spinner(worker, set_status_fn, initial_msg: str):
    """Animate a spinner in the status bar while *worker* runs.

    Returns ``(timer, update_msg_fn)`` where ``update_msg_fn(msg)`` changes the
    animated text (connect to ``worker.yielded`` for stage messages).
    """
    from qtpy.QtCore import QTimer
    _FRAMES = ["|", "/", "-", "\\"]
    _idx = [0]
    _msg = [initial_msg]

    timer = QTimer()

    def _tick():
        set_status_fn(f"{_FRAMES[_idx[0] % 4]} {_msg[0]}")
        _idx[0] += 1

    def update_msg(msg: str):
        _msg[0] = msg

    timer.timeout.connect(_tick)
    timer.start(150)
    worker.finished.connect(timer.stop)
    return timer, update_msg


def attach_tqdm_progress(worker, set_status_fn, base_msg: str = ""):
    """Wire a ProgressMailbox + QTimer to relay tqdm updates to the status bar.

    Returns a ``post_fn`` callable safe to call from the background thread;
    pass it into :func:`qt_tqdm_context` inside the worker.
    """
    from qtpy.QtCore import QTimer
    mailbox = ProgressMailbox()

    timer = QTimer()

    def _poll():
        msg = mailbox.read()
        if msg:
            set_status_fn(msg)

    timer.timeout.connect(_poll)
    timer.start(100)
    worker.finished.connect(timer.stop)
    return mailbox.post


from contextlib import contextmanager


@contextmanager
def qt_tqdm_context(post_fn, base_msg: str = ""):
    """Context manager (run inside a background thread) that routes tqdm
    progress updates to *post_fn* instead of the terminal.
    """
    import tqdm as _tqdm_module
    import tqdm.auto as _tqdm_auto
    import io

    original_tqdm = _tqdm_module.tqdm
    original_auto = _tqdm_auto.tqdm

    class _StatusTqdm(original_tqdm):
        def __init__(self, *args, **kwargs):
            kwargs.setdefault('file', io.StringIO())  # suppress stderr output
            super().__init__(*args, **kwargs)

        def update(self, n=1):
            super().update(n)
            if self.total:
                pct = int(100 * self.n / self.total)
                post_fn(f"{base_msg}{self.n}/{self.total} ({pct}%)")
            else:
                post_fn(f"{base_msg}{self.n} iterations")

    _tqdm_module.tqdm = _StatusTqdm
    _tqdm_auto.tqdm = _StatusTqdm
    try:
        yield
    finally:
        _tqdm_module.tqdm = original_tqdm
        _tqdm_auto.tqdm = original_auto


# ── Shared helper factory ────────────────────────────────────────────────────

def create_shared_helpers(ctx: ViewerContext):
    """Create cross-tab helper functions and attach them to *ctx*.

    Must be called after all tab widgets have been registered on *ctx*
    (clustering_widget, ga_clustering_widget, etc.).
    """
    state = ctx.state

    # ── apply_plot_font_size ──────────────────────────────────────────────
    def _apply_plot_font_size():
        import matplotlib.pyplot as _plt
        _plt.rcParams['font.size'] = state.get("plot_font_size", 10)

    ctx.apply_plot_font_size = _apply_plot_font_size

    # ── set_status ────────────────────────────────────────────────────────
    def _set_status(msg: str):
        ctx.viewer.status = msg

    ctx.set_status = _set_status

    # ── record_code ──────────────────────────────────────────────────────
    def _record_code(code: str, tag: str = None):
        if not state.get("record_code"):
            return
        if tag:
            if tag in state["code_journal_tags"]:
                return
            state["code_journal_tags"].add(tag)
        state["code_journal"].append(code)
        code_path = ctx.data_path / "code.py"
        with open(code_path, 'w') as f:
            f.write("\n".join(state["code_journal"]) + "\n")

    ctx.record_code = _record_code

    # ── record_preamble ──────────────────────────────────────────────────
    def _record_preamble():
        _record_code(
            "import scanpy as sc\n"
            "import squidpy as sq\n"
            "import matplotlib.pyplot as plt\n"
            "import pandas as pd\n"
            "import numpy as np\n"
            f"\nplt.rcParams['font.size'] = {state.get('plot_font_size', 10)}\n"
            f"\n# Load data\n"
            "from spatialdata_io import xenium\n"
            f"sdata = xenium(\"{ctx.data_path}\")\n"
            "adata = sdata[\"table\"].copy()",
            tag="preamble"
        )

    ctx.record_preamble = _record_preamble

    # ── record_normalize ─────────────────────────────────────────────────
    def _record_normalize():
        _record_preamble()
        _record_code(
            "\n# Normalize, log-transform, PCA\n"
            "sc.pp.normalize_total(adata)\n"
            "sc.pp.log1p(adata)\n"
            "sc.pp.pca(adata)",
            tag="normalize"
        )

    ctx.record_normalize = _record_normalize

    # ── record_clustering ────────────────────────────────────────────────
    def _record_clustering(key):
        _record_normalize()
        dir_name = f"gene_expression_{key}"
        csv_path = os.path.join(ctx.data_path, "analysis", "clustering", dir_name, "clusters.csv")
        if not os.path.exists(csv_path):
            csv_path = os.path.join(ctx.data_path, "analysis", "clustering", key, "clusters.csv")
        _record_code(
            f"\n# Add clustering: {key}\n"
            f"clust_df = pd.read_csv(\"{csv_path}\", index_col=0)\n"
            f"adata.obs[\"{key}\"] = pd.Categorical("
            f"clust_df.reindex(adata.obs_names).iloc[:, 0].astype(str).values)",
            tag=f"clustering_{key}"
        )

    ctx.record_clustering = _record_clustering

    # ── record_spatial_neighbors ─────────────────────────────────────────
    def _record_spatial_neighbors(n_neighs):
        _record_code(
            f"\n# Compute spatial neighbors (k={n_neighs})\n"
            "adata.obsm['spatial'] = adata.obsm.get('spatial', "
            "np.column_stack([adata.obs['x_centroid'], adata.obs['y_centroid']]))\n"
            f"sq.gr.spatial_neighbors(adata, n_neighs={n_neighs}, coord_type=\"generic\")",
            tag=f"spatial_neighbors_{n_neighs}"
        )

    ctx.record_spatial_neighbors = _record_spatial_neighbors

    # ── get_plot_save_path ───────────────────────────────────────────────
    def _get_plot_save_path(title: str, default_stem: str) -> str | None:
        from qtpy.QtWidgets import QFileDialog
        fmt = state.get("plot_format", "png")
        filter_str = f"{fmt.upper()} Files (*.{fmt});;All Files (*)"
        path, _ = QFileDialog.getSaveFileName(
            None, title, f"{default_stem}.{fmt}", filter_str,
        )
        return path if path else None

    ctx.get_plot_save_path = _get_plot_save_path

    # ── refresh_clustering_choices ───────────────────────────────────────
    def _refresh_clustering_choices():
        names = list(ctx.clusterings.keys())
        for combo in [ctx.clustering_widget, ctx.ga_clustering_widget,
                      ctx.lr_clustering_widget, ctx.ne_clustering_widget,
                      ctx.co_clustering_widget]:
            if combo is None:
                continue
            old_val = combo.value
            combo.choices = names
            if old_val in names:
                combo.value = old_val

    ctx.refresh_clustering_choices = _refresh_clustering_choices

    # ── repopulate_cluster_checkboxes ────────────────────────────────────
    def _repopulate_cluster_checkboxes():
        from qtpy.QtWidgets import QCheckBox
        grid = ctx.cluster_filter_grid
        while grid.count():
            item = grid.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        state["cluster_checkboxes"].clear()

        key = ctx.clustering_widget.value
        if not key or key not in ctx.clusterings:
            return
        raw_ids = ctx.clusterings[key].dropna().unique().tolist()
        try:
            ids = sorted([int(x) for x in raw_ids])
        except (ValueError, TypeError):
            ids = sorted(raw_ids, key=lambda x: str(x))
        try:
            labels = _get_labels_for(key)
        except Exception:
            labels = {}
        cols = 3
        for i, cid in enumerate(ids):
            display = str(labels.get(cid, labels.get(str(cid), cid)))
            cb = QCheckBox(display)
            cb.setChecked(True)
            cb.setEnabled(ctx.filter_check.value)
            grid.addWidget(cb, i // cols, i % cols)
            state["cluster_checkboxes"][cid] = cb

    ctx.repopulate_cluster_checkboxes = _repopulate_cluster_checkboxes

    # ── get_selected_cluster_ids ─────────────────────────────────────────
    def _get_selected_cluster_ids():
        return {cid for cid, cb in state["cluster_checkboxes"].items() if cb.isChecked()}

    ctx.get_selected_cluster_ids = _get_selected_cluster_ids

    # ── make_cluster_mask ────────────────────────────────────────────────
    def _make_cluster_mask(aligned_values, selected_ids):
        sel = {str(s) for s in selected_ids}
        return np.array([str(v) in sel for v in aligned_values], dtype=bool)

    ctx.make_cluster_mask = _make_cluster_mask

    # ── cluster label helpers ────────────────────────────────────────────
    if "cluster_labels" not in state or not isinstance(state.get("cluster_labels"), dict):
        state["cluster_labels"] = {}

    def _get_active_labels():
        key = state.get("active_clustering_name") or ctx.clustering_widget.value
        if key:
            all_labels = state.get("cluster_labels", {})
            if isinstance(all_labels, dict):
                return all_labels.get(key, {})
        return {}

    ctx.get_active_labels = _get_active_labels

    def _get_labels_for(clustering_key):
        all_labels = state.get("cluster_labels", {})
        if isinstance(all_labels, dict):
            return all_labels.get(clustering_key, {})
        return {}

    ctx.get_labels_for = _get_labels_for

    # ── build_label_editor_dialog ────────────────────────────────────────
    def _build_label_editor_dialog(clustering_key):
        from qtpy.QtWidgets import (
            QDialog, QGridLayout, QLineEdit, QDialogButtonBox,
            QScrollArea, QLabel as QtLabel,
        )
        if not clustering_key or clustering_key not in ctx.clusterings:
            return False
        cluster_series = ctx.clusterings[clustering_key]
        ids = sorted(cluster_series.dropna().unique().tolist(), key=lambda x: (str(x),))
        existing = _get_labels_for(clustering_key)

        dialog = QDialog()
        dialog.setWindowTitle(f"Edit Cluster Labels \u2014 {clustering_key}")

        n_cols = min(3, max(1, (len(ids) + 9) // 10))
        n_per_col = (len(ids) + n_cols - 1) // n_cols

        outer_layout = QVBoxLayout()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_content = QWidget()
        grid = QGridLayout()

        edits = {}
        for i, cid in enumerate(ids):
            col_idx = i // n_per_col
            row_idx = i % n_per_col
            grid.addWidget(QtLabel(f"{cid}:"), row_idx, col_idx * 2)
            edit = QLineEdit(str(existing.get(cid, existing.get(str(cid), cid))))
            grid.addWidget(edit, row_idx, col_idx * 2 + 1)
            edits[cid] = edit

        scroll_content.setLayout(grid)
        scroll.setWidget(scroll_content)
        outer_layout.addWidget(scroll)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        outer_layout.addWidget(buttons)
        dialog.setLayout(outer_layout)
        dialog.resize(min(800, 250 * n_cols), min(600, 30 * n_per_col + 60))

        if dialog.exec_() == QDialog.Accepted:
            new_labels = {cid: e.text() for cid, e in edits.items()}
            if "cluster_labels" not in state or not isinstance(state["cluster_labels"], dict):
                state["cluster_labels"] = {}
            state["cluster_labels"][clustering_key] = new_labels
            return True
        return False

    ctx.build_label_editor_dialog = _build_label_editor_dialog

    # ── get_cluster_ids_per_obs ──────────────────────────────────────────
    def _get_cluster_ids_per_obs(clustering_key):
        import pandas as pd
        cluster_series = ctx.clusterings[clustering_key]
        _adata = ctx.color_manager.adata
        if 'cell_id' in _adata.obs.columns:
            cell_ids = _adata.obs['cell_id'].values
            clusters_aligned = cluster_series.reindex(cell_ids)
        else:
            clusters_aligned = cluster_series.reindex(_adata.obs_names)

        try:
            filled = clusters_aligned.fillna(-1)
            cluster_values = filled.values.astype(np.int32)
            state['_cluster_id_to_raw'] = None
            state['_cluster_raw_to_id'] = None
        except (ValueError, TypeError):
            codes, uniques = pd.factorize(clusters_aligned.values)
            cluster_values = codes.astype(np.int32)
            state['_cluster_id_to_raw'] = {int(i): u for i, u in enumerate(uniques)}
            state['_cluster_raw_to_id'] = {u: int(i) for i, u in enumerate(uniques)}

        label_to_obs = ctx.label_to_obs
        max_label = len(label_to_obs) - 1
        label_to_cluster = np.full(max_label + 1, -1, dtype=np.int32)
        valid_mask = label_to_obs >= 0
        valid_labels = np.where(valid_mask)[0]
        obs_indices = label_to_obs[valid_labels]
        label_to_cluster[valid_labels] = cluster_values[obs_indices]
        return cluster_values, label_to_cluster

    ctx.get_cluster_ids_per_obs = _get_cluster_ids_per_obs

    # ── translate_selected_ids_to_int ────────────────────────────────────
    def _translate_selected_ids_to_int(selected_ids):
        raw_to_id = state.get('_cluster_raw_to_id')
        if raw_to_id is None:
            return list(selected_ids)
        return [raw_to_id[sid] for sid in selected_ids if sid in raw_to_id]

    ctx.translate_selected_ids_to_int = _translate_selected_ids_to_int

    # ── get_cluster_filter ───────────────────────────────────────────────
    def _get_cluster_filter():
        if not state.get("filter_by_cluster"):
            return None
        selected = _get_selected_cluster_ids()
        if not selected:
            return None
        return sorted(str(cid) for cid in selected)

    ctx.get_cluster_filter = _get_cluster_filter


# ── File menu ────────────────────────────────────────────────────────────────

def create_file_menu(ctx: ViewerContext, on_open_dataset, on_preprocess_dataset=None):
    """Build a File menu on the napari menu bar with an Open Dataset action."""
    from qtpy.QtWidgets import QMenu
    from qtpy.QtGui import QAction

    menu_bar = ctx.viewer.window._qt_window.menuBar()
    file_menu = QMenu("File", menu_bar)
    existing = menu_bar.actions()
    if existing:
        menu_bar.insertMenu(existing[0], file_menu)  # insert before Preferences
    else:
        menu_bar.addMenu(file_menu)

    act = QAction("Open Dataset...", file_menu)
    act.setShortcut("Ctrl+O")
    file_menu.addAction(act)
    act.triggered.connect(on_open_dataset)

    if on_preprocess_dataset is not None:
        file_menu.addSeparator()
        act2 = QAction("Preprocess Dataset...", file_menu)
        file_menu.addAction(act2)
        act2.triggered.connect(on_preprocess_dataset)


# ── Preferences menu ─────────────────────────────────────────────────────────

def create_preferences_menu(ctx: ViewerContext):
    """Build the Preferences menu on the napari menu bar."""
    from qtpy.QtWidgets import QActionGroup, QMenu
    from qtpy.QtGui import QAction

    state = ctx.state
    menu_bar = ctx.viewer.window._qt_window.menuBar()
    prefs_menu = QMenu("Preferences", menu_bar)
    menu_bar.addMenu(prefs_menu)

    # Plot format
    format_menu = prefs_menu.addMenu("Plot format")
    format_group = QActionGroup(format_menu)
    format_group.setExclusive(True)
    png_action = QAction("PNG", format_group, checkable=True, checked=True)
    svg_action = QAction("SVG", format_group, checkable=True)
    format_menu.addAction(png_action)
    format_menu.addAction(svg_action)

    def _on_format_changed(action):
        state["plot_format"] = action.text().lower()
    format_group.triggered.connect(_on_format_changed)

    # Font size
    fontsize_menu = prefs_menu.addMenu("Plot font size")
    fontsize_group = QActionGroup(fontsize_menu)
    fontsize_group.setExclusive(True)
    for sz in (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20):
        act = QAction(str(sz), fontsize_group, checkable=True, checked=(sz == 10))
        fontsize_menu.addAction(act)

    def _on_fontsize_changed(action):
        state["plot_font_size"] = int(action.text())
    fontsize_group.triggered.connect(_on_fontsize_changed)

    # Record code checkbox
    record_action = QAction("Record reproducible code", prefs_menu, checkable=True, checked=True)
    prefs_menu.addAction(record_action)

    def _on_record_toggled(checked):
        state["record_code"] = checked
        if checked:
            state["code_journal"].clear()
            state["code_journal_tags"].clear()
    record_action.toggled.connect(_on_record_toggled)

    # Save code action
    save_code_action = QAction("Save recorded code...", prefs_menu)
    prefs_menu.addAction(save_code_action)

    def _on_save_code():
        if not state["code_journal"]:
            return
        from qtpy.QtWidgets import QFileDialog as _QFD
        path, _ = _QFD.getSaveFileName(
            None, "Save Reproducible Code", "analysis.py", "Python Files (*.py)",
        )
        if path:
            with open(path, 'w') as f:
                f.write("\n".join(state["code_journal"]) + "\n")
    save_code_action.triggered.connect(_on_save_code)
