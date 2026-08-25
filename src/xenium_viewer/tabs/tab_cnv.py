"""Tab: CNV inference — inferCNV (in-process) and CopyKAT (detached background)."""

from __future__ import annotations
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from magicgui.widgets import CheckBox, ComboBox, PushButton, SpinBox, FloatSpinBox
from qtpy.QtCore import QTimer
from qtpy.QtWidgets import (
    QTextEdit, QHBoxLayout, QWidget, QLabel, QScrollArea, QGridLayout, QCheckBox,
)
from napari.qt.threading import thread_worker
from xenium_viewer.tabs._helpers import make_tab, StatusProxy, attach_spinner, make_progress_bar, combo_value_kwargs
from xenium_viewer.utils.prov_graph import ARTIFACT, TERMINAL
from xenium_viewer.utils.steps import Step, coerce
from xenium_viewer.utils.step_templates import (
    Preview, builtin_assemble, builtin_spec, step_template as _resolved,
)

_BACKEND_LABELS = {"infercnv": "inferCNV", "copykat": "CopyKAT"}


# inferCNV runs in-process, so it is a real Step: one string, executed and
# recorded. The old recorded cell had drifted in three places — it normalised at
# scanpy's median default rather than the target_sum=1e4 the viewer used, and it
# dropped `lfc_clip` and `dendrogram=False`.
#
# CopyKAT is *not* migrated: it runs detached, in a second conda env, so no
# in-process step can be the code that ran. Its node stays on record_node and is
# labelled as a description of that run rather than as executed source.



# pandas 3 backs strings with PyArrow arrays, which do not support the
# numpy-style fancy indexing infercnvpy does on var_names / var columns:
# `_running_mean` slices the gene list with a *2-D* index array, and
# ArrowStringArray routes that to pyarrow's take(), which raises
# "ArrowInvalid: only handle 1-dimensional arrays".
#
# `future.infer_string` must be off across the *assignment*, not just the
# conversion: AnnData re-infers string dtypes when a frame is assigned to
# .obs/.var, so converting to object and assigning it back under the default
# option leaves the Arrow array exactly where it was. This mirrors
# adata_persistence._convert_adata_arrow_strings, which run_cnv_pipeline uses —
# the two must stay identical (tests/test_cnv_step.py).



CNV_TEMPLATE_ID = "genes.cnv_infercnv"


def _cnv_blocks(subset: bool) -> list[str]:
    # arrow_shim is frozen in the .tmpl: its byte identity with
    # adata_persistence._convert_adata_arrow_strings is pinned by a test, and no
    # validation gate could warn a user that editing it breaks on pandas 3.
    return (["head"] + (["subset"] if subset else [])
            + ["prepare", "arrow_shim", "tail"])


def _cnv_template(subset: bool) -> str:
    return builtin_assemble(CNV_TEMPLATE_ID, _cnv_blocks(subset))


#: The pandas-3 Arrow workaround on its own. Named because its byte identity
#: with ``adata_persistence._convert_adata_arrow_strings`` is an invariant a
#: test executes and pins — two implementations of one conversion that must not
#: drift. Marked non-editable in the .tmpl for the same reason.
_CNV_ARROW_SHIM = builtin_spec(CNV_TEMPLATE_ID).blocks["arrow_shim"].text


def _backend_label(backend: str) -> str:
    return _BACKEND_LABELS.get(backend, backend)


def _copykat_worker_pid(running_marker: "Path") -> "int | None":
    """PID recorded in a CopyKAT RUNNING marker, or ``None`` (legacy/no PID).

    Newer markers hold ``{"pid": ..., "timestamp": ...}``; older ones held a bare
    ISO timestamp string, from which no PID can be recovered.
    """
    try:
        return int(json.loads(running_marker.read_text())["pid"])
    except (OSError, ValueError, TypeError, KeyError):
        return None


def _copykat_run_state(plots_dir: "Path") -> "tuple[str, int | None]":
    """Classify a detached CopyKAT run from its marker files.

    Returns one of:
      * ``("running", pid)``     — worker PID is alive and is our worker;
      * ``("interrupted", pid)`` — marker present but the worker is gone (killed);
      * ``("unknown", None)``    — marker present but no PID to check (legacy);
      * ``("idle", None)``       — no in-progress marker (absent, or already DONE).

    Liveness is confirmed against ``/proc/<pid>/cmdline`` (so a reused PID belonging
    to an unrelated process is treated as gone), falling back to an ``os.kill(pid, 0)``
    probe where ``/proc`` is unavailable.
    """
    running = plots_dir / "copykat_RUNNING.txt"
    done = plots_dir / "copykat_DONE.txt"
    if not running.exists() or done.exists():
        return ("idle", None)
    pid = _copykat_worker_pid(running)
    if pid is None:
        return ("unknown", None)
    # The /proc route is Linux-only, and the platform check has to come *first*:
    # on macOS `open("/proc/...")` raises FileNotFoundError, which is an OSError
    # subclass caught by the "no such process" branch below — so a live worker was
    # reported as "interrupted" there rather than falling through to os.kill().
    alive = None
    if sys.platform.startswith("linux"):
        try:
            with open(f"/proc/{pid}/cmdline", "rb") as f:
                cmdline = f.read().replace(b"\x00", b" ").decode("utf-8", "replace")
            # Confirms identity too, so a reused PID belonging to an unrelated
            # process is treated as gone.
            alive = "cnv_copykat_worker" in cmdline
        except FileNotFoundError:
            alive = False  # no such process
        except OSError:
            alive = None   # /proc unreadable (permissions, container) — probe below
    if alive is None:
        try:
            os.kill(pid, 0)
            alive = True
        except ProcessLookupError:
            alive = False
        except PermissionError:
            alive = True  # exists but owned by another user
    return ("running", pid) if alive else ("interrupted", pid)


def _resolve_copykat_prefix() -> "Path | None":
    """Locate the CopyKAT conda env *prefix* (a separate py3.11 + R env).

    CopyKAT's R stack is incompatible with the viewer's python 3.12, so the
    detached worker runs in a second env (default name ``xenium_viewer_copykat``,
    created from environment-copykat.yml). It is launched via ``conda run -p
    <prefix>`` so conda activates the env's R (rpy2 needs the env's R_HOME, not
    the system R). Resolution order:
      1. ``XENIUM_COPYKAT_PYTHON`` — path to that env's python (prefix inferred);
      2. a sibling conda env named by ``XENIUM_COPYKAT_ENV`` (default
         ``xenium_viewer_copykat``) under the same conda base as this env.
    Returns the env prefix Path, or ``None`` if not found.
    """
    explicit = os.environ.get("XENIUM_COPYKAT_PYTHON")
    if explicit and Path(explicit).exists():
        return Path(explicit).resolve().parents[1]  # <prefix>/bin/python -> <prefix>
    env_name = os.environ.get("XENIUM_COPYKAT_ENV", "xenium_viewer_copykat")
    # This env is <base>/envs/<this>; the copykat env is a sibling <base>/envs/<env_name>.
    cand = Path(sys.prefix).parent / env_name
    if (cand / "bin" / "python").exists():
        return cand
    return None

if TYPE_CHECKING:
    from xenium_viewer.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    cnv_clustering_widget = ComboBox(
        label="Reference clustering", choices=ctx.clustering_names,
        **combo_value_kwargs(ctx.clustering_names),
    )
    ctx.cnv_clustering_widget = cnv_clustering_widget

    cnv_ref_label = QLabel("Select reference (\"normal\") clusters:")
    cnv_select_all_btn = PushButton(label="Select All")
    cnv_deselect_all_btn = PushButton(label="Deselect All")

    cnv_cluster_container = QWidget()
    cnv_cluster_grid = QGridLayout()
    cnv_cluster_grid.setContentsMargins(0, 0, 0, 0)
    cnv_cluster_container.setLayout(cnv_cluster_grid)

    cnv_cluster_scroll = QScrollArea()
    cnv_cluster_scroll.setWidget(cnv_cluster_container)
    cnv_cluster_scroll.setWidgetResizable(True)
    cnv_cluster_scroll.setMaximumHeight(150)

    state["cnv_reference_checkboxes"] = {}

    def _cnv_select_all():
        for cb in state["cnv_reference_checkboxes"].values():
            cb.setChecked(True)

    def _cnv_deselect_all():
        for cb in state["cnv_reference_checkboxes"].values():
            cb.setChecked(False)

    cnv_select_all_btn.clicked.connect(_cnv_select_all)
    cnv_deselect_all_btn.clicked.connect(_cnv_deselect_all)

    def _repopulate_cnv_reference_checkboxes():
        grid = cnv_cluster_grid
        while grid.count():
            item = grid.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        state["cnv_reference_checkboxes"].clear()

        key = cnv_clustering_widget.value
        if not key or key not in ctx.clusterings:
            return
        raw_ids = ctx.clusterings[key].dropna().unique().tolist()
        try:
            ids = sorted([int(x) for x in raw_ids])
        except (ValueError, TypeError):
            ids = sorted(raw_ids, key=lambda x: str(x))
        labels = ctx.get_labels_for(key) if ctx.get_labels_for else {}
        cols = 3
        for i, cid in enumerate(ids):
            display = str(labels.get(cid, labels.get(str(cid), cid)))
            cb = QCheckBox(display)
            cb.setChecked(False)  # reference population must be explicitly chosen
            grid.addWidget(cb, i // cols, i % cols)
            state["cnv_reference_checkboxes"][cid] = cb

    cnv_clustering_widget.changed.connect(lambda _: _repopulate_cnv_reference_checkboxes())
    _repopulate_cnv_reference_checkboxes()

    def _get_cnv_reference_ids():
        """Return list of checked reference cluster IDs (str)."""
        cbs = state["cnv_reference_checkboxes"]
        return [str(cid) for cid, cb in cbs.items() if cb.isChecked()]

    # ── Cell types to analyze (limit CNV to a subset) ───────────────────
    cnv_analyze_label = QLabel("Cell types to analyze (CNV subclones):")
    cnv_analyze_hint = QLabel(
        "Only these cell types plus the reference are included in the analysis;\n"
        "leave all checked to analyze the whole tissue."
    )
    cnv_analyze_select_all_btn = PushButton(label="Select All")
    cnv_analyze_deselect_all_btn = PushButton(label="Deselect All")

    cnv_analyze_container = QWidget()
    cnv_analyze_grid = QGridLayout()
    cnv_analyze_grid.setContentsMargins(0, 0, 0, 0)
    cnv_analyze_container.setLayout(cnv_analyze_grid)

    cnv_analyze_scroll = QScrollArea()
    cnv_analyze_scroll.setWidget(cnv_analyze_container)
    cnv_analyze_scroll.setWidgetResizable(True)
    cnv_analyze_scroll.setMaximumHeight(150)

    state["cnv_analyze_checkboxes"] = {}

    def _cnv_analyze_select_all():
        for cb in state["cnv_analyze_checkboxes"].values():
            cb.setChecked(True)

    def _cnv_analyze_deselect_all():
        for cb in state["cnv_analyze_checkboxes"].values():
            cb.setChecked(False)

    cnv_analyze_select_all_btn.clicked.connect(_cnv_analyze_select_all)
    cnv_analyze_deselect_all_btn.clicked.connect(_cnv_analyze_deselect_all)

    def _repopulate_cnv_analyze_checkboxes():
        grid = cnv_analyze_grid
        while grid.count():
            item = grid.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        state["cnv_analyze_checkboxes"].clear()

        key = cnv_clustering_widget.value
        if not key or key not in ctx.clusterings:
            return
        raw_ids = ctx.clusterings[key].dropna().unique().tolist()
        try:
            ids = sorted([int(x) for x in raw_ids])
        except (ValueError, TypeError):
            ids = sorted(raw_ids, key=lambda x: str(x))
        labels = ctx.get_labels_for(key) if ctx.get_labels_for else {}
        cols = 3
        for i, cid in enumerate(ids):
            display = str(labels.get(cid, labels.get(str(cid), cid)))
            cb = QCheckBox(display)
            cb.setChecked(True)  # analyze all cell types by default
            grid.addWidget(cb, i // cols, i % cols)
            state["cnv_analyze_checkboxes"][cid] = cb

    _repopulate_cnv_analyze_checkboxes()

    def _get_cnv_analyze_ids():
        """Return list of checked analyze cluster IDs (str)."""
        cbs = state["cnv_analyze_checkboxes"]
        return [str(cid) for cid, cb in cbs.items() if cb.isChecked()]

    def _all_cluster_ids(key):
        """All (str) cluster IDs present in clustering ``key``."""
        if not key or key not in ctx.clusterings:
            return set()
        return {str(x) for x in ctx.clusterings[key].dropna().unique().tolist()}

    # The reference grid is already rebuilt on dropdown change (see above);
    # rebuild the analyze grid on the same signal.
    cnv_clustering_widget.changed.connect(lambda _: _repopulate_cnv_analyze_checkboxes())

    # ── Parameters ──────────────────────────────────────────────────────
    # Defaults match InSituCNV's own reference notebook (run_insitucnv.ipynb).
    cnv_n_neighbors = SpinBox(label="Neighbors (expression graph)", min=5, max=100, value=15)
    cnv_smoothing_neighbors = SpinBox(label="Smoothing neighbors", min=5, max=200, value=20)
    cnv_window_size = SpinBox(label="Window size (genes)", min=2, max=200, value=60)
    cnv_step = SpinBox(label="Window step", min=1, max=50, value=10)
    cnv_resolution = FloatSpinBox(label="CNV cluster resolution", min=0.05, max=2.0, step=0.05, value=0.2)
    cnv_resolution.tooltip = (
        "InSituCNV's own notebook evaluates several resolutions (e.g. 0.1, 0.2, 0.3) "
        "and picks one per dataset after reviewing the results — this default may not "
        "be right for your data."
    )
    cnv_resolution_hint = QLabel(
        "Default may need tuning per dataset — check the chromosome heatmap and\n"
        "cluster count after running, and re-run with a different value if clusters\n"
        "look too coarse or too fragmented."
    )

    cnv_backend_widget = ComboBox(
        label="CNV backend",
        choices=[("Both (inferCNV + CopyKAT)", "both"), ("inferCNV only", "infercnv"),
                 ("CopyKAT only", "copykat")],
        value="both",
    )
    cnv_max_cells = SpinBox(label="CopyKAT max cells", min=500, max=500000, value=10000)
    cnv_max_cells.tooltip = (
        "CopyKAT is slow (~2 h/sample) and runs as a detached background job on a random "
        "subsample of at most this many cells (all reference cells are kept). inferCNV uses all cells."
    )
    cnv_extrapolate = CheckBox(label="Extrapolate CopyKAT calls to all cells", value=False)
    cnv_extrapolate.tooltip = (
        "After CopyKAT finishes on the subsample, extend its tumor/normal (cnv_status), "
        "copykat_pred, and CNV-subclone (copykat_leiden_res*) results to every cell. Cells CopyKAT "
        "actually ran keep their real value; each un-run cell is filled with the majority value "
        "among run cells in its reference-clustering group. Filled values are copied cluster-level "
        "values, not per-cell inferred CNV; un-run cells in a group with no run cell are labelled "
        "'unknown'. Adds colorable '<col>_propagated' clusterings for each."
    )

    run_button = PushButton(label="Run CNV Inference", enabled=True)
    heatmap_backend_widget = ComboBox(label="Heatmap backend", choices=[], **combo_value_kwargs([]))
    heatmap_res_widget = ComboBox(label="Heatmap resolution", choices=[], **combo_value_kwargs([]))
    heatmap_button = PushButton(label="Save Chromosome Heatmap (PDF/PNG)", enabled=False)
    score_color_button = PushButton(label="Color Cells by CNV Score", enabled=False)

    state.setdefault("cnv_results", {})
    state.setdefault("_cnv_bg_jobs", [])

    def _key_to_res_label(key: str) -> str:
        """'cnv_leiden_res0.2' -> 'res 0.2' (falls back to the raw key)."""
        m = re.search(r"res([0-9.]+)$", str(key))
        return f"res {m.group(1)}" if m else str(key)

    def _res_from_key(key):
        m = re.search(r"res([0-9.]+)$", str(key))
        try:
            return float(m.group(1)) if m else None
        except (TypeError, ValueError):
            return None

    def _refresh_heatmap_res_choices(select=None):
        """Populate the resolution combo from the selected backend's result."""
        backend = heatmap_backend_widget.value
        result = state.get("cnv_results", {}).get(backend)
        keys = list(result.get("cluster_keys", [])) if result else []
        adata_cnv = result.get("adata_cnv") if result else None
        if adata_cnv is not None:
            cols = set(adata_cnv.obs.columns)
            keys = [k for k in keys if k in cols]
        heatmap_res_widget.choices = [(_key_to_res_label(k), k) for k in keys]
        if select in keys:
            heatmap_res_widget.value = select
        elif keys:
            heatmap_res_widget.value = keys[-1]
        heatmap_button.enabled = bool(keys) and adata_cnv is not None
        score_color_button.enabled = bool(result and result.get("cnv_score") is not None)

    def _refresh_heatmap_backend_choices(select=None):
        results = state.get("cnv_results", {})
        backends = [b for b in ("infercnv", "copykat") if b in results]
        heatmap_backend_widget.choices = [(_backend_label(b), b) for b in backends]
        if select in backends:
            heatmap_backend_widget.value = select
        elif backends:
            heatmap_backend_widget.value = backends[-1]
        _refresh_heatmap_res_choices()

    heatmap_backend_widget.changed.connect(lambda _: _refresh_heatmap_res_choices())

    def _cnv_signature(result):
        """Core params a CNV profile depends on (resolution excluded, backend included).

        Must stay term-for-term identical to _cnv_signature_from_info in
        adata_persistence.py: leading backend, then core params, trailing analyzed set.
        """
        p = result["params"]
        return (
            result.get("backend", "infercnv"),
            result["reference_clustering_name"],
            tuple(result["reference_categories"]),
            p.get("n_neighbors"),
            p.get("smoothing_neighbors"),
            p.get("window_size"),
            p.get("step"),
            p.get("lfc_clip"),
            tuple(sorted(result.get("analyze_categories") or [])),
        )

    results_text = QTextEdit()
    results_text.setReadOnly(True)
    results_text.setFontFamily("monospace")
    results_text.setMaximumHeight(220)

    cnv_status = StatusProxy(ctx.viewer)
    cnv_progress = make_progress_bar()

    def _align_score_to_obs_order(score_series):
        adata = ctx.adata
        if 'cell_id' in adata.obs.columns:
            cell_ids = adata.obs['cell_id'].values
            aligned = score_series.reindex(cell_ids)
        else:
            aligned = score_series.reindex(adata.obs_names)
        return aligned.to_numpy(dtype=float)

    # ── Provenance + status helpers ─────────────────────────────────────
    def _record_cnv_node(result):
        """Record the CopyKAT run.

        inferCNV is not recorded here — it runs through ``ctx.run_step``, which
        records the source it executed. CopyKAT cannot: it runs detached, in a
        second conda env (``xenium_viewer_copykat``), because its R stack needs
        python 3.11. No in-process step can be the code that ran, so this node is
        a *reconstruction* of that run rather than executed source, and it is
        labelled as such in the notebook.
        """
        backend = result.get("backend", "infercnv")
        if backend != "copykat":
            return
        ref_name = result["reference_clustering_name"]
        ctx.record_clustering(ref_name)
        p = result["params"]
        ref_obs = result["reference_obs_key"]
        ref_repr = repr(result["reference_categories"])
        var = "adata_copykat"
        analyze_cats = result.get("analyze_categories") or []
        L = [
            "\n# CNV inference (CopyKAT)",
            "# NOTE: the viewer ran this in a detached process in the",
            "# 'xenium_viewer_copykat' conda env (R + rpy2 need python 3.11), so",
            "# unlike every other cell this is a reconstruction of that run.",
            "from insitucnv.tl import prepare_cnv_input, compute_cnv_neighbors, cluster_cnv_resolutions, run_copykat",
            "import scanpy as sc",
            f"{var} = adata.copy()",
        ]
        if analyze_cats:
            include_repr = repr(sorted(set(analyze_cats) | set(result["reference_categories"])))
            L.append(f"{var} = {var}[{var}.obs['{ref_name}'].astype(str).isin({include_repr})].copy()"
                     f"  # limit to selected cell types + reference")
        max_cells = result.get("max_cells") or p.get("max_cells") or result.get("n_cells")
        L.append(f"if {var}.n_obs > {max_cells}: {var} = sc.pp.subsample({var}, n_obs={max_cells}, "
                 f"random_state=0, copy=True)  # CopyKAT is slow — subsample")
        L += [
            f"{var}.obs['{ref_obs}'] = {var}.obs['{ref_name}']  # reference clustering",
            f"{var}.layers['raw_counts'] = {var}.X.copy()",
            f"sc.pp.normalize_total({var}, target_sum=1e4); sc.pp.log1p({var}); sc.pp.pca({var})",
            f"sc.pp.neighbors({var}, n_neighbors={p['n_neighbors']})",
            f"{var} = prepare_cnv_input({var}, raw_layer='raw_counts', smoothing_neighbors="
            f"{p['smoothing_neighbors']}, add_gene_positions=True, drop_unmapped_genes=True, copy=False)",
        ]
        n_cores = p.get("n_cores", 1)
        L.append(f"run_copykat({var}, reference_key='{ref_obs}', reference_categories={ref_repr}, "
                 f"input_layer='M', n_cores={n_cores}, copy=False)  # fork defaults (genome hg20, win_size 25, ...)")
        prefix = "copykat_leiden_res"
        L += [
            f"compute_cnv_neighbors({var}, copy=False)",
            f"cluster_cnv_resolutions({var}, {result['resolutions']!r}, key_prefix='{prefix}', "
            f"dendrogram=False, copy=False)",
            f"# result: {var}.obs[{result['cluster_keys']!r}]",
        ]
        ctx.record_node("cnv:copykat", "\n".join(L), deps=[f"clustering:{ref_name}"],
                        label="CNV inference (CopyKAT)")

    def _write_results_text(result, note=""):
        backend = result.get("backend", "infercnv")
        ref_name = result["reference_clustering_name"]
        series = result["cluster_series"]
        analyze_cats = result.get("analyze_categories") or []
        if analyze_cats:
            labels = ctx.get_labels_for(ref_name) if ctx.get_labels_for else {}
            analyze_line = ", ".join(str(labels.get(c, labels.get(str(c), c))) for c in analyze_cats)
            cells_line = f"  Cells analyzed: {result.get('n_cells')} (cell types: {analyze_line})"
        else:
            cells_line = f"  Cells analyzed: {result.get('n_cells')} (all cell types)"
        if backend == "copykat":
            cells_line += f" [CopyKAT subsample ≤ {result.get('max_cells') or result.get('n_cells')}]"
        res_line = ", ".join(str(r) for r in result.get("resolutions", []))
        try:
            score_rng = f"{result['cnv_score'].min():.3f} - {result['cnv_score'].max():.3f}"
        except Exception:
            score_rng = "n/a"
        results_text.setPlainText(
            f"CNV inference complete — {_backend_label(backend)}{note}\n"
            f"  Reference clustering: {ref_name}\n"
            f"  Reference clusters: {', '.join(result['reference_categories'])}\n"
            f"{cells_line}\n"
            f"  Genes: {result.get('n_genes_mapped')} / {result.get('n_genes_total')} mapped to genome\n"
            f"  CNV windows: {result.get('n_windows')}\n"
            f"  CNV clusters found (res {result['params'].get('resolution')}): {series.nunique()}\n"
            f"  CNV score range: {score_rng}\n"
            f"\nResult stored as: {result['cluster_key']}\n"
            f"  Resolutions available for heatmap: {res_line}\n"
            f"Pick a backend + resolution below to save its chromosome heatmap, or color by CNV score."
        )

    # ── Shared ingest: fold one backend's result into the registry ──────
    def _ingest_cnv_result(result, note=""):
        backend = result.get("backend", "infercnv")
        key = result["cluster_key"]
        series = result["cluster_series"]

        # Accumulate resolutions across same-signature runs of THIS backend.
        sig = _cnv_signature(result)
        prev = state["cnv_results"].get(backend)
        param_change_note = note
        if (prev is not None and prev.get("signature") == sig
                and prev.get("adata_cnv") is not None):
            shared = prev["adata_cnv"]
            shared.obs[key] = (
                pd.Series(result["adata_cnv"].obs[key].values,
                          index=result["adata_cnv"].obs_names)
                .reindex(shared.obs_names).values
            )
            result["adata_cnv"] = shared
            result["cluster_keys"] = list(dict.fromkeys(prev.get("cluster_keys", []) + [key]))
        else:
            result["cluster_keys"] = list(dict.fromkeys(result.get("cluster_keys") or [key]))
            if prev is not None and prev.get("signature") not in (None, sig):
                param_change_note = (note or "") + " (parameters changed — previous resolutions cleared)"
        result["signature"] = sig
        result["resolutions"] = sorted(
            {r for r in (_res_from_key(k) for k in result["cluster_keys"]) if r is not None}
        )
        state["cnv_results"][backend] = result
        state["cnv_result"] = result  # alias to the most-recently-updated backend

        # Register the new clustering(s) so Cell Coloring lists them.
        ctx.color_manager.invalidate_cluster_cache(key)
        ctx.clusterings[key] = series
        state.setdefault("custom_clusterings", {})[key] = series
        ctx.refresh_clustering_choices()
        _refresh_heatmap_backend_choices(select=backend)
        _refresh_heatmap_res_choices(select=key)

        # Auto-apply cluster coloring in a background thread.
        @thread_worker
        def _apply_colors():
            color_arr, cluster_to_color = ctx.color_manager.get_cluster_colors(series)
            cluster_ids_per_obs, label_to_cluster = ctx.get_cluster_ids_per_obs(key)
            colormap = ctx.color_manager.build_direct_label_colormap(color_arr)
            return colormap, color_arr, cluster_to_color, label_to_cluster, cluster_ids_per_obs

        def _on_colors_ready(color_result):
            colormap, color_arr, cluster_to_color, label_to_cluster, cluster_ids_per_obs = color_result
            state["cluster_to_color"] = cluster_to_color
            state["label_to_cluster"] = label_to_cluster
            state["active_clustering_name"] = key
            if ctx.cell_labels_layer is not None:
                ctx.cell_labels_layer.colormap = colormap
                ctx.cell_labels_layer.refresh()
            ctx.umap_viewer.color_by_cluster(
                key, color_arr, ctx.label_to_obs, cluster_ids_per_obs=cluster_ids_per_obs,
            )

        color_worker = _apply_colors()
        color_worker.returned.connect(_on_colors_ready)
        color_worker.start()

        # Persist + record provenance + status.
        from xenium_viewer.utils.adata_persistence import save_clustering_to_adata, save_cnv_results_to_adata
        save_clustering_to_adata(ctx, key, series)
        save_cnv_results_to_adata(ctx, result)
        _record_cnv_node(result)
        _record_cnv_clustering_node(result, key)
        _write_results_text(result, param_change_note)

    def _record_cnv_clustering_node(result, key):
        """Record ``clustering:<key>`` for the CNV-derived cluster column.

        The CNV node computes these labels on ``adata_cnv`` and binds
        ``cnv_clusters``; the viewer then publishes them onto the main table
        with ``save_clustering_to_adata``. Without a node of that id, any tab
        that later analysed this clustering fell through to the generic
        ``record_clustering`` fallback, which used to emit a ``read_csv`` of an
        ``analysis/clustering/`` file that was never written — so the exported
        notebook failed outright. This is the code the viewer actually ran,
        expressed against what the CNV cell already bound.
        """
        backend = result.get("backend", "infercnv")
        head = ("\n# Publish the CNV clustering onto the main table\n"
                "_ids = (adata.obs['cell_id'].values\n"
                "        if 'cell_id' in adata.obs.columns else adata.obs_names)\n")
        if backend == "copykat":
            # CopyKAT leaves its labels on the (subsampled) adata_copykat, so the
            # column covers only the cells it ran on; the rest stay unlabelled.
            var = "adata_copykat"
            body = (
                f"_ck_ids = ({var}.obs['cell_id'].values\n"
                f"           if 'cell_id' in {var}.obs.columns else {var}.obs_names)\n"
                f"_ck = pd.Series({var}.obs[{key!r}].values, index=_ck_ids)\n"
                f"adata.obs[{key!r}] = pd.Categorical(_ck.reindex(_ids).values)"
            )
        else:
            # The inferCNV cell binds ``cnv_clusters``, indexed by cell_id.
            body = f"adata.obs[{key!r}] = pd.Categorical(cnv_clusters.reindex(_ids).values)"
        ctx.record_node(
            f"clustering:{key}", head + body,
            deps=[f"cnv:{backend}"],
            label=f"Clustering: {key}",
        )

    def _result_from_copykat_output(sidecar, adata_cnv):
        cluster_key = sidecar["cluster_key"]
        series = adata_cnv.obs[cluster_key].copy()
        if "cell_id" in adata_cnv.obs.columns:
            series.index = adata_cnv.obs["cell_id"].values
        series.name = cluster_key
        idx = adata_cnv.obs["cell_id"].values if "cell_id" in adata_cnv.obs.columns else adata_cnv.obs_names
        score = (pd.Series(adata_cnv.obs["cnv_score"].values, index=idx, name="cnv_score")
                 if "cnv_score" in adata_cnv.obs.columns else None)
        return {
            "backend": "copykat", "adata_cnv": adata_cnv, "cluster_key": cluster_key,
            "cluster_series": series, "cnv_score": score,
            "cluster_keys": list(sidecar.get("cluster_keys", [cluster_key])),
            "reference_obs_key": sidecar.get("reference_obs_key"),
            "reference_clustering_name": sidecar.get("reference_clustering_name", ""),
            "reference_categories": list(sidecar.get("reference_categories", [])),
            "analyze_categories": list(sidecar.get("analyze_categories", [])),
            "n_cells": sidecar.get("n_cells"), "n_genes_total": sidecar.get("n_genes_total"),
            "n_genes_mapped": sidecar.get("n_genes_mapped"), "n_windows": sidecar.get("n_windows"),
            "max_cells": sidecar.get("max_cells"), "params": dict(sidecar.get("params", {})),
            "extrapolate": bool(sidecar.get("extrapolate", False)),
        }

    def _record_extrapolation_node(result, label_cols, ref_key):
        backend = result.get("backend", "copykat")
        var = f"adata_{backend}"
        keys_repr = repr(tuple(label_cols))
        prop_cols = [f"{c}_propagated" for c in label_cols]
        L = [
            "\n# Extrapolate CopyKAT results to all cells (fill un-run cells; keep run cells' real values)",
            "from insitucnv.tl import propagate_cnv_labels",
            f"# adata = full dataset; {var} = CopyKAT-run subsample (from the CNV step above)",
            f"adata.obs['{ref_key}'] = adata.obs['{ref_key}'].astype(str)",
            f"propagate_cnv_labels(adata, {var}, label_keys={keys_repr},",
            f"    method='cluster', cluster_key='{ref_key}', suffix='_propagated', copy=False)",
            f"# result: adata.obs[{prop_cols!r}] (all cells; un-run empty groups -> 'unknown')",
            "# Each propagated column is finished off by its own cell below, which",
            "# overlays the run cells' real values.",
        ]
        ctx.record_node("cnv:copykat_propagated", "\n".join(L), deps=[f"cnv:{backend}"],
                        label="Extrapolate CopyKAT calls (all cells)")
        # One clustering node per propagated column, rather than a loop hidden
        # inside the node above: each is a clustering the user can select, so
        # each needs an id that dependents (rank genes, nhood, ...) can name.
        for col in label_cols:
            pkey = f"{col}_propagated"
            ctx.record_node(
                f"clustering:{pkey}",
                f"\n# Overlay each CopyKAT-run cell's real {col!r} so its true call/subclone\n"
                f"# isn't collapsed to the reference-cluster majority (extrapolation only\n"
                f"# fills the cells CopyKAT never ran):\n"
                f"_real = {var}.obs[{col!r}].astype(str)\n"
                f"adata.obs[{pkey!r}] = adata.obs[{pkey!r}].astype(str)\n"
                f"adata.obs.loc[_real.index, {pkey!r}] = _real.values",
                deps=["cnv:copykat_propagated"],
                label=f"Clustering: {pkey}",
            )

    def _extrapolate_copykat(result):
        """Propagate CopyKAT results from the run subsample to ALL cells.

        Cells CopyKAT actually ran keep their real value; only the un-run cells are
        filled, by majority CopyKAT value among run cells in their reference-clustering
        group (``propagate_cnv_labels(method="cluster")``, then the run cells' real
        values overlaid back on top). Propagates both the per-cell calls (``cnv_status``,
        ``copykat_pred``) and the CNV-subclone clustering(s) (``copykat_leiden_res*``),
        registering each as a full-coverage ``<col>_propagated`` colorable clustering.
        Filled values are copied cluster-level values, not per-cell inferred CNV;
        un-run cells in a reference group with no run cell are labelled ``unknown``.
        """
        ref_key = result.get("reference_clustering_name")
        adata_cnv = result.get("adata_cnv")
        if not ref_key or ref_key not in ctx.clusterings or adata_cnv is None:
            ctx.set_status("Extrapolation skipped: reference clustering unavailable.")
            return
        call_cols = [c for c in ("cnv_status", "copykat_pred") if c in adata_cnv.obs.columns]
        cluster_cols = [k for k in (result.get("cluster_keys") or []) if k in adata_cnv.obs.columns]
        label_cols = call_cols + cluster_cols
        if not label_cols:
            ctx.set_status("Extrapolation skipped: no CopyKAT result columns found.")
            return
        try:
            from insitucnv.tl import propagate_cnv_labels
        except Exception as e:  # noqa: BLE001
            ctx.set_status(f"Extrapolation unavailable (insitucnv import failed): {e}")
            return
        import anndata as _ad

        ref_full = ctx.clusterings[ref_key].astype(str)  # indexed by cell_id, all cells
        sub_idx = [str(c) for c in (adata_cnv.obs["cell_id"].values
                   if "cell_id" in adata_cnv.obs.columns else adata_cnv.obs_names)]
        labeled_obs = pd.DataFrame({ref_key: ref_full.reindex(sub_idx).astype(str).values},
                                   index=sub_idx)
        for c in label_cols:
            labeled_obs[c] = [str(v) for v in adata_cnv.obs[c].values]
        adata_labeled = _ad.AnnData(obs=labeled_obs)
        full_ids = [str(c) for c in ref_full.index]
        adata_full = _ad.AnnData(obs=pd.DataFrame({ref_key: ref_full.values}, index=full_ids))

        propagate_cnv_labels(
            adata_full, adata_labeled, label_keys=tuple(label_cols),
            method="cluster", cluster_key=ref_key, suffix="_propagated", copy=False,
        )

        from xenium_viewer.utils.adata_persistence import save_clustering_to_adata
        n_filled = 0
        for c in label_cols:
            pkey = f"{c}_propagated"
            series = pd.Series(adata_full.obs[pkey].astype(str).values, index=full_ids, name=pkey)
            # Extrapolation should only *fill in* the cells CopyKAT never ran; keep each
            # SAMPLED cell's real value instead of overwriting it with its reference-cluster
            # majority (which would collapse the true 0/1/2/3 subclones to the dominant one).
            real = pd.Series([str(v) for v in adata_cnv.obs[c].values], index=sub_idx)
            common = real.index.intersection(series.index)
            series.loc[common] = real.loc[common]
            n_filled = int(series.index.difference(common).size)
            ctx.color_manager.invalidate_cluster_cache(pkey)
            ctx.clusterings[pkey] = series
            state.setdefault("custom_clusterings", {})[pkey] = series
            save_clustering_to_adata(ctx, pkey, series)
        ctx.refresh_clustering_choices()

        n_labeled = len(set(sub_idx) & set(full_ids))
        ctx.set_status(
            f"Extrapolated CopyKAT results: kept {n_labeled} run cells' real calls, "
            f"filled {n_filled} un-run cells by reference-cluster majority "
            f"({', '.join(label_cols)})."
        )
        _record_extrapolation_node(result, label_cols, ref_key)

    # ── Background-job polling (detached CopyKAT worker) ────────────────
    def _stop_poll_timer():
        timer = state.pop("_cnv_poll_timer", None)
        if timer is not None:
            timer.stop()

    def _ensure_poll_timer():
        if state.get("_cnv_poll_timer") is None:
            timer = QTimer()
            timer.setInterval(5000)
            timer.timeout.connect(_poll_bg_jobs)
            timer.start()
            state["_cnv_poll_timer"] = timer

    def _poll_bg_jobs():
        jobs = state.get("_cnv_bg_jobs", [])
        for job in list(jobs):
            done_file = job["done_file"]
            if not done_file.exists():
                proc = job.get("proc")
                if proc is not None and proc.poll() is not None:
                    jobs.remove(job)
                    # The worker died without finishing (crash/kill); its finally-block
                    # cleanup didn't run, so clear the stale RUNNING marker ourselves.
                    marker = job.get("running_marker")
                    if marker is not None:
                        try:
                            Path(marker).unlink(missing_ok=True)
                        except OSError:
                            pass
                    cnv_status.value = "CopyKAT background job exited without a result (see log)."
                continue
            jobs.remove(job)
            try:
                status = json.loads(done_file.read_text())
            except Exception:
                status = {"status": "unknown"}
            if status.get("status") != "ok":
                cnv_status.value = f"CopyKAT background analysis failed: {status.get('error', 'see log')}"
                results_text.setPlainText(f"CopyKAT failed:\n{status.get('error', '')}")
                continue
            if ctx.dataset_generation != job.get("gen"):
                continue  # dataset switched while it ran
            try:
                import scanpy as sc
                sidecar = json.loads(Path(job["result_json"]).read_text())
                adata_cnv = sc.read_h5ad(job["out_h5ad"])
                result = _result_from_copykat_output(sidecar, adata_cnv)
                _ingest_cnv_result(result, note=" (background run)")
                if result.get("extrapolate"):
                    try:
                        _extrapolate_copykat(result)
                    except Exception as e:  # noqa: BLE001 — extrapolation is best-effort
                        cnv_status.value = f"CopyKAT done; extrapolation failed: {e}"
                cnv_status.value = f"CopyKAT background analysis complete ({result.get('n_cells')} cells)."
            except Exception as e:
                cnv_status.value = f"CopyKAT result load error: {e}"
        if not jobs:
            _stop_poll_timer()

    # ── Run dispatch ────────────────────────────────────────────────────
    def _infercnv_preview(reference_key=None, reference_ids=None,
                          analyze_categories=None) -> Preview:
        """What an inferCNV run would execute with the widgets as they stand.

        One expression of the current settings, called by ``_run_infercnv`` with
        the selections it has validated, and by the Templates tab's preview pane
        with none. Restricting the analysis to a subset of cell types selects a
        block as well as filling a param, so both halves travel together.

        inferCNV only. CopyKAT runs detached in another conda env and stays on
        ``record_node`` (its cell says in line that it is a reconstruction), so
        there is no in-process template for this pane to preview.
        """
        if reference_key is None:
            reference_key = cnv_clustering_widget.value
            reference_ids = _get_cnv_reference_ids()
            analyze_ids = _get_cnv_analyze_ids()
            analyze_categories = (
                None if not analyze_ids or set(analyze_ids) >= _all_cluster_ids(reference_key)
                else list(analyze_ids)
            )
        from xenium_viewer.utils.cnv_analysis import CNV_REFERENCE_OBS_KEY

        reference_categories = [str(c) for c in (reference_ids or [])]
        analyze_cats = [str(c) for c in analyze_categories] if analyze_categories else []
        params = {
            "reference_clustering": reference_key,
            "reference_obs_key": CNV_REFERENCE_OBS_KEY,
            "reference_categories": reference_categories,
            "n_neighbors": coerce(cnv_n_neighbors.value),
            "smoothing_neighbors": coerce(cnv_smoothing_neighbors.value),
            "window_size": coerce(cnv_window_size.value),
            "step": coerce(cnv_step.value),
            "lfc_clip": 4.0,
            "resolution": coerce(cnv_resolution.value),
        }
        if analyze_cats:
            params["include"] = sorted(set(analyze_cats) | set(reference_categories))
        return Preview(_cnv_blocks(bool(analyze_cats)), params)

    ctx.state.setdefault(
        "template_preview", {})[CNV_TEMPLATE_ID] = _infercnv_preview

    def _run_infercnv(reference_key, reference_ids, analyze_categories):
        from xenium_viewer.utils.cnv_analysis import (
            CNV_REFERENCE_OBS_KEY, _patch_matplotlib_cm_compat,
        )
        from xenium_viewer.utils.gene_analysis import add_clustering_to_obs

        _adata = ctx.adata
        n_genes_total = _adata.n_vars
        reference_categories = [str(c) for c in reference_ids]
        analyze_cats = [str(c) for c in analyze_categories] if analyze_categories else []

        # The reference clustering must exist as a node, and in adata.obs, before
        # the step can declare it as a dependency and read it. Both on the GUI thread.
        ctx.record_clustering(reference_key)
        add_clustering_to_obs(_adata, _adata, ctx.clusterings[reference_key], reference_key)

        blocks, step_params, _ = _infercnv_preview(
            reference_key, reference_ids, analyze_categories)
        # The result summary reports the settings back to the user; n_cores is a
        # runtime detail of this machine, not part of the recorded analysis.
        params = {k: step_params[k] for k in (
            "n_neighbors", "smoothing_neighbors", "window_size", "step",
            "lfc_clip", "resolution")}
        params["n_cores"] = 1

        step = Step(
            id="cnv:infercnv",
            **_resolved(CNV_TEMPLATE_ID, blocks),
            params=step_params,
            deps=[f"clustering:{reference_key}"],
            kind=ARTIFACT,
            label="CNV inference (inferCNV)",
            outputs=["adata_cnv", "cnv_cluster_keys", "cnv_clusters", "cnv_score"],
        )

        gen = ctx.dataset_generation
        run_button.enabled = False

        @thread_worker
        def _run():
            _patch_matplotlib_cm_compat()
            out = ctx.run_step(step)
            adata_cnv = out["adata_cnv"]
            cluster_key = out["cnv_cluster_keys"][0]
            return {
                "adata_cnv": adata_cnv,
                "cluster_key": cluster_key,
                "cluster_series": out["cnv_clusters"],
                "cnv_score": out["cnv_score"],
                "n_genes_total": int(n_genes_total),
                "n_genes_mapped": int(adata_cnv.n_vars),
                "n_windows": int(adata_cnv.obsm["X_cnv"].shape[1]),
                "n_cells": int(adata_cnv.n_obs),
                "backend": "infercnv",
                "reference_obs_key": CNV_REFERENCE_OBS_KEY,
                "reference_clustering_name": reference_key,
                "reference_categories": reference_categories,
                "analyze_categories": analyze_cats,
                "params": params,
            }

        worker = _run()
        _timer, _ = attach_spinner(
            worker, lambda m: setattr(cnv_status, "value", m),
            "Running inferCNV...", progress_bar=cnv_progress,
        )
        state["_cnv_spinner_timer"] = _timer

        def _done(result):
            run_button.enabled = True
            if ctx.dataset_generation != gen:
                return
            _ingest_cnv_result(result)

        def _err(exc):
            run_button.enabled = True
            cnv_status.value = f"inferCNV error: {exc}"
            results_text.setPlainText(f"Error running inferCNV:\n{exc}")

        worker.returned.connect(_done)
        worker.errored.connect(_err)
        worker.start()

    def _launch_copykat(reference_key, reference_ids, analyze_ids, analyze_categories):
        if ctx.no_cache or ctx.sdata is None or ctx.sdata.path is None:
            return False, "CopyKAT needs the zarr cache — not available with --no-cache."
        copykat_prefix = _resolve_copykat_prefix()
        if copykat_prefix is None:
            return False, (
                "CopyKAT backend needs the separate 'xenium_viewer_copykat' env. Create it with: "
                "conda env create -f environment-copykat.yml  (or set XENIUM_COPYKAT_PYTHON)."
            )
        conda_exe = os.environ.get("CONDA_EXE") or "conda"
        from xenium_viewer.utils.cnv_analysis import subsample_indices

        import anndata as _ad

        adata = ctx.adata
        reference_series = ctx.clusterings[reference_key]
        obs_index = adata.obs["cell_id"].values if "cell_id" in adata.obs.columns else adata.obs_names
        keep = subsample_indices(reference_series, reference_ids, analyze_categories,
                                 obs_index, int(cnv_max_cells.value))
        sub = adata[keep]
        # Write a MINIMAL, clean input h5ad (raw counts + cell_id + reference column
        # + gene names only). The copykat env pins an older anndata (python 3.11), so
        # extra obs/obsm/uns from the main env's newer anndata can use encodings it
        # can't read (e.g. all-null columns -> encoding_type='null').
        cell_ids = [str(c) for c in (sub.obs["cell_id"].values if "cell_id" in sub.obs.columns else sub.obs_names)]
        ref_src = [str(v) for v in reference_series.reindex(cell_ids).to_numpy()]
        names = [str(x) for x in sub.obs_names]
        obs = pd.DataFrame({"cell_id": cell_ids, "_cnv_ref_src": ref_src}, index=names)
        var = pd.DataFrame(index=[str(g) for g in sub.var_names])
        X = sub.X
        subset = _ad.AnnData(X=(X.copy() if hasattr(X, "copy") else X), obs=obs, var=var)

        # Beside the zarr store, not inside it. The detached worker writes here
        # while the GUI is live; files in the store root make zarr's hierarchy
        # walk warn on every consolidation, and a cache rebuild would delete
        # hours of CopyKAT compute along with the cache.
        from xenium_viewer.utils.adata_persistence import sidecar_dir
        cache_dir = sidecar_dir(ctx.data_path, create=True)
        plots_dir = Path(ctx.data_path) / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        input_h5ad = cache_dir / "cnv_copykat_input.h5ad"
        params_json = cache_dir / "cnv_copykat_params.json"
        done_file = plots_dir / "copykat_DONE.txt"
        running_marker = plots_dir / "copykat_RUNNING.txt"
        result_json = cache_dir / "cnv_copykat_result.json"
        out_h5ad = cache_dir / "adata_cnv_cache_copykat.h5ad"
        for stale in (done_file, running_marker, result_json):
            try:
                stale.unlink(missing_ok=True)
            except Exception:
                pass
        try:
            from xenium_viewer.utils.adata_persistence import _convert_adata_arrow_strings
            _convert_adata_arrow_strings(subset)
        except Exception:
            pass
        subset.write_h5ad(input_h5ad)

        params = {
            "reference_categories": list(reference_ids),
            "reference_clustering_name": reference_key,
            "reference_obs_col": "_cnv_ref_src",
            "analyze_categories": list(analyze_categories) if analyze_categories else [],
            "resolution": cnv_resolution.value,
            "n_neighbors": cnv_n_neighbors.value,
            "smoothing_neighbors": cnv_smoothing_neighbors.value,
            "window_size": cnv_window_size.value,
            "step": cnv_step.value,
            "lfc_clip": 4.0,
            "n_cores": int(ctx.state.get("n_cores", 1)),
            "extrapolate": bool(cnv_extrapolate.value),
            "max_cells": int(cnv_max_cells.value),
            "n_input_cells": int(subset.n_obs),
            "output_h5ad": str(out_h5ad),
            "result_json": str(result_json),
            "plots_dir": str(plots_dir),
            "done_file": str(done_file),
            "running_marker": str(running_marker),
            "copykat_workdir": str(cache_dir / "copykat"),
        }
        params_json.write_text(json.dumps(params, indent=2))

        # Run the worker in the CopyKAT env via `conda run -p <prefix>` so conda
        # activates that env's R (rpy2 needs the env's R_HOME, not the system R).
        # Make THIS repo's xenium_viewer importable there via PYTHONPATH, so the
        # copykat env needn't install the viewer package. <pkg>/../ is the src root.
        import xenium_viewer as _xv
        src_root = str(Path(_xv.__file__).resolve().parents[1])
        child_env = dict(os.environ)
        child_env["PYTHONPATH"] = os.pathsep.join(
            [src_root] + ([child_env["PYTHONPATH"]] if child_env.get("PYTHONPATH") else [])
        )
        log_file = open(plots_dir / "copykat_worker.log", "w")
        proc = subprocess.Popen(
            [conda_exe, "run", "-p", str(copykat_prefix), "--no-capture-output",
             "python", "-m", "xenium_viewer.cnv_copykat_worker",
             str(input_h5ad), str(params_json), str(cache_dir)],
            stdout=log_file, stderr=subprocess.STDOUT, start_new_session=True, env=child_env,
        )
        try:
            pgid = os.getpgid(proc.pid)
        except Exception:
            pgid = proc.pid
        state["_cnv_bg_jobs"].append({
            "backend": "copykat", "proc": proc, "pgid": pgid,
            "done_file": done_file, "running_marker": running_marker,
            "result_json": result_json, "out_h5ad": out_h5ad,
            "gen": ctx.dataset_generation, "started_at": time.time(),
        })
        _ensure_poll_timer()
        return True, f"CopyKAT started in the background on {subset.n_obs} cells (~2 h)."

    def on_run_cnv():
        if ctx.adata is None:
            cnv_status.value = "No AnnData loaded"
            return
        reference_key = cnv_clustering_widget.value
        reference_ids = _get_cnv_reference_ids()
        if not reference_key or reference_key not in ctx.clusterings:
            cnv_status.value = "Select a reference clustering first"
            return
        if not reference_ids:
            cnv_status.value = "Select at least one reference (\"normal\") cluster"
            return
        analyze_ids = _get_cnv_analyze_ids()
        if not analyze_ids:
            cnv_status.value = "Select at least one cell type to analyze"
            return
        analyze_categories = (
            None if set(analyze_ids) >= _all_cluster_ids(reference_key) else list(analyze_ids)
        )

        backends = {"both": ["infercnv", "copykat"], "infercnv": ["infercnv"],
                    "copykat": ["copykat"]}[cnv_backend_widget.value]
        started = []
        if "infercnv" in backends:
            results_text.setPlainText("Running inferCNV...")
            _run_infercnv(reference_key, reference_ids, analyze_categories)
            started.append("inferCNV")
        if "copykat" in backends:
            ok, msg = _launch_copykat(reference_key, reference_ids, analyze_ids, analyze_categories)
            cnv_status.value = msg
            if ok:
                started.append("CopyKAT (background)")
        if started:
            cnv_status.value = "Started: " + ", ".join(started)

    def on_show_heatmap():
        backend = heatmap_backend_widget.value
        result = state.get("cnv_results", {}).get(backend)
        if result is None:
            cnv_status.value = "No CNV result to plot"
            return
        adata_cnv = result.get("adata_cnv")
        if adata_cnv is None:
            cnv_status.value = "No cached CNV profile available — rerun CNV inference"
            return
        cluster_key = heatmap_res_widget.value or result["cluster_key"]
        if cluster_key not in adata_cnv.obs.columns:
            cnv_status.value = (
                f"Resolution '{_key_to_res_label(cluster_key)}' is not in the "
                f"{_backend_label(backend)} profile — re-run at that resolution"
            )
            return

        cnv_status.value = f"Building {_backend_label(backend)} heatmap ({_key_to_res_label(cluster_key)})..."
        heatmap_button.enabled = False
        safe = re.sub(r"[^0-9A-Za-z._-]", "_", cluster_key)
        # inferCNV binds adata_cnv (the step's own name); the CopyKAT
        # reconstruction node still uses adata_copykat.
        var = "adata_cnv" if backend == "infercnv" else "adata_copykat"

        @thread_worker
        def _build():
            from xenium_viewer.utils.cnv_analysis import make_cnv_heatmap
            ctx.apply_plot_font_size()
            fig = make_cnv_heatmap(adata_cnv, cluster_key)
            plots_dir = os.path.join(ctx.data_path, "plots")
            os.makedirs(plots_dir, exist_ok=True)
            png_path = os.path.join(plots_dir, f"cnv_heatmap_{backend}_{safe}.png")
            pdf_path = os.path.join(plots_dir, f"cnv_heatmap_{backend}_{safe}.pdf")
            fig.savefig(png_path, dpi=200, bbox_inches="tight")
            fig.savefig(pdf_path, bbox_inches="tight")
            import matplotlib.pyplot as _plt
            _plt.close(fig)
            return png_path, pdf_path

        def _on_ready(paths):
            heatmap_button.enabled = True
            png_path, pdf_path = paths
            cnv_status.value = f"{_backend_label(backend)} heatmap saved to {png_path} and {pdf_path}"
            ctx.record_node(
                f"plot:cnv_heatmap:{backend}:{cluster_key}",
                f"\n# CNV chromosome heatmap ({_backend_label(backend)}, {_key_to_res_label(cluster_key)})\n"
                f"import infercnvpy as cnv\n"
                f"cnv.pl.chromosome_heatmap({var}, groupby='{cluster_key}', dendrogram=True, "
                f"vmin=-0.4, vmax=0.4, show=False)\n"
                f"plt.savefig('cnv_heatmap_{backend}_{safe}.png', dpi=200, bbox_inches='tight')\n"
                f"plt.savefig('cnv_heatmap_{backend}_{safe}.pdf', bbox_inches='tight')",
                deps=[f"cnv:{backend}"],
                kind=TERMINAL,
                label=f"CNV heatmap ({_backend_label(backend)}, {_key_to_res_label(cluster_key)})",
            )

        def _on_error(exc):
            heatmap_button.enabled = True
            cnv_status.value = f"Heatmap error: {exc}"

        worker = _build()
        worker.returned.connect(_on_ready)
        worker.errored.connect(_on_error)
        worker.start()

    def on_color_by_score():
        backend = heatmap_backend_widget.value
        result = state.get("cnv_results", {}).get(backend)
        if result is None:
            return
        score_series = result.get("cnv_score")
        if score_series is None:
            cnv_status.value = "No CNV score available — rerun CNV inference"
            return
        cnv_status.value = f"Coloring cells by {_backend_label(backend)} CNV score..."
        score_color_button.enabled = False

        @thread_worker
        def _build():
            values = _align_score_to_obs_order(score_series)
            color_arr = ctx.color_manager.get_continuous_colors(
                values, colormap="viridis", cache_key=f"cnv_score_{backend}",
            )
            return ctx.color_manager.build_direct_label_colormap(color_arr)

        def _on_ready(colormap):
            score_color_button.enabled = True
            if ctx.cell_labels_layer is not None:
                ctx.cell_labels_layer.colormap = colormap
                ctx.cell_labels_layer.refresh()
            cnv_status.value = f"Cells colored by {_backend_label(backend)} CNV score (viridis)"

        def _on_error(exc):
            score_color_button.enabled = True
            cnv_status.value = f"CNV score coloring error: {exc}"

        worker = _build()
        worker.returned.connect(_on_ready)
        worker.errored.connect(_on_error)
        worker.start()

    run_button.clicked.connect(on_run_cnv)
    heatmap_button.clicked.connect(on_show_heatmap)
    score_color_button.clicked.connect(on_color_by_score)

    cnv_sel_btn_row = QWidget()
    cnv_sel_btn_layout = QHBoxLayout()
    cnv_sel_btn_layout.setContentsMargins(0, 0, 0, 0)
    cnv_sel_btn_layout.addWidget(cnv_select_all_btn.native)
    cnv_sel_btn_layout.addWidget(cnv_deselect_all_btn.native)
    cnv_sel_btn_row.setLayout(cnv_sel_btn_layout)

    cnv_analyze_btn_row = QWidget()
    cnv_analyze_btn_layout = QHBoxLayout()
    cnv_analyze_btn_layout.setContentsMargins(0, 0, 0, 0)
    cnv_analyze_btn_layout.addWidget(cnv_analyze_select_all_btn.native)
    cnv_analyze_btn_layout.addWidget(cnv_analyze_deselect_all_btn.native)
    cnv_analyze_btn_row.setLayout(cnv_analyze_btn_layout)

    widget = make_tab(
        cnv_clustering_widget,
        cnv_ref_label,
        cnv_sel_btn_row,
        cnv_cluster_scroll,
        cnv_analyze_label,
        cnv_analyze_hint,
        cnv_analyze_btn_row,
        cnv_analyze_scroll,
        cnv_n_neighbors,
        cnv_smoothing_neighbors,
        cnv_window_size,
        cnv_step,
        cnv_resolution,
        cnv_resolution_hint,
        cnv_backend_widget,
        cnv_max_cells,
        cnv_extrapolate,
        run_button,
        cnv_progress,
        results_text,
        heatmap_backend_widget,
        heatmap_res_widget,
        heatmap_button,
        score_color_button,
    )

    def _restore_session(session):
        results = state.get("cnv_results") or {}
        if not results:
            return
        # Ensure every restored backend's cluster columns are colorable, even a
        # background CopyKAT run the GUI never folded into the main table.
        for backend, result in results.items():
            adata_cnv = result.get("adata_cnv")
            if adata_cnv is None:
                continue
            from xenium_viewer.utils.adata_persistence import save_clustering_to_adata
            for key in result.get("cluster_keys", []):
                if key in ctx.clusterings or key not in adata_cnv.obs.columns:
                    continue
                series = adata_cnv.obs[key].copy()
                if "cell_id" in adata_cnv.obs.columns:
                    series.index = adata_cnv.obs["cell_id"].values
                series.name = key
                ctx.clusterings[key] = series
                state.setdefault("custom_clusterings", {})[key] = series
                try:
                    save_clustering_to_adata(ctx, key, series)
                except Exception:
                    pass
            result["resolutions"] = sorted(
                {r for r in (_res_from_key(k) for k in result.get("cluster_keys", [])) if r is not None}
            )
        ctx.refresh_clustering_choices()

        default_b = "infercnv" if "infercnv" in results else next(iter(results))
        _refresh_heatmap_backend_choices(select=default_b)
        state["cnv_result"] = results.get(default_b)

        lines = ["CNV inference (restored from previous session)"]
        for backend, result in results.items():
            res_line = ", ".join(str(r) for r in result.get("resolutions", []))
            analyze_cats = result.get("analyze_categories") or []
            scope = (", ".join(analyze_cats) if analyze_cats else "all cell types")
            lines.append(f"  [{_backend_label(backend)}] ref={result.get('reference_clustering_name')} "
                         f"cells={result.get('n_cells')} ({scope}); resolutions: {res_line or '(none)'}")
        # Surface a possibly-still-running detached CopyKAT job, and clear the
        # stale marker a killed worker leaves behind (its finally-block cleanup
        # never runs on SIGTERM/SIGKILL).
        try:
            plots_dir = Path(ctx.data_path) / "plots"
            st, pid = _copykat_run_state(plots_dir)
            if st == "running":
                lines.append(f"  ⚠ a background CopyKAT job is still running (pid {pid}).")
            elif st == "interrupted":
                try:
                    (plots_dir / "copykat_RUNNING.txt").unlink(missing_ok=True)
                except OSError:
                    pass
                lines.append("  ⚠ a previous CopyKAT job was interrupted; cleared its stale marker.")
            elif st == "unknown":
                lines.append("  ⚠ a background CopyKAT job appears to be in progress or was interrupted.")
        except Exception:
            pass
        lines.append("\nPick a backend + resolution below to save a chromosome heatmap, "
                     "or color cells by CNV score.")
        results_text.setPlainText("\n".join(lines))
        print(f"  Restored CNV results (backends: {list(results.keys())})")

    return widget, {"cnv_clustering_widget": cnv_clustering_widget, "restore_session": _restore_session}
