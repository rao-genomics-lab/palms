"""Tools -> QC: the standard Xenium quality-control panel, and cell/gene filtering.

Two things, in the order a user needs them. First a look at the data --- count
distributions and the negative-control rates --- so a cutoff is chosen from a
picture rather than from habit. Then the filter itself, which drops near-empty
cells and rarely-detected genes and **rebinds the cell set every later analysis
is about**.

That second part is why this is its own tab rather than a group box in
Clustering, where issue #77 suggested it. A filter applied here changes what
Rank Genes, CNV, ROI DEG, the neighbourhood statistics and the exported notebook
are all about; a control with that reach should not sit inside one analysis.
Tools -> Segmentation is its neighbour for the same reason: they are the two
places that answer "which cells?".

Nothing is filtered until Apply. The cutoffs are pre-filled with the
conventional Xenium starting points, and a dataset that has never been through
this tab is untouched --- no node in the provenance graph, nothing stale.
"""
from __future__ import annotations

import numpy as np
from magicgui.widgets import CheckBox, PushButton, SpinBox
from napari.qt.threading import thread_worker
from qtpy.QtCore import QTimer
from qtpy.QtWidgets import QLabel, QTextEdit
from superqt.utils import ensure_main_thread

from palms.tabs._helpers import (
    StatusProxy, attach_spinner, make_progress_bar, make_tab, qc_filter_preview,
    toolbar_row,
)
from palms.utils.prov_graph import TERMINAL
from palms.utils.step_templates import Preview, builtin_assemble, step_template as _resolved
from palms.utils.steps import Step, coerce
from palms.utils.viewer_context import ViewerContext

FILTER_TEMPLATE_ID = "qc.filter"
METRICS_TEMPLATE_ID = "qc.metrics"

#: 10x's own Xenium starting points, and what scanpy's tutorials use.
DEFAULT_MIN_COUNTS = 10
DEFAULT_MIN_CELLS = 3

#: The positions ``sc.pp.calculate_qc_metrics`` is asked to summarise. Clipped
#: to the panel before use --- scanpy raises "Positions outside range of
#: features" for any position past ``n_vars``, and a custom panel, or one after
#: gene filtering, need not reach 150.
PERCENT_TOP = (10, 20, 50, 150)

#: The obs columns the four-panel figure needs beyond the count metrics. Xenium
#: writes both; a custom segmentation table need not.
AREA_COLUMNS = ("cell_area", "nucleus_area")


def percent_top_for(n_vars: int) -> tuple:
    """The subset of :data:`PERCENT_TOP` a panel of *n_vars* genes supports."""
    return tuple(p for p in PERCENT_TOP if p <= int(n_vars))


def _metrics_blocks(has_areas: bool) -> list[str]:
    """Which blocks the available obs columns call for."""
    return (["head", "controls", "plot4", "areas", "save"] if has_areas
            else ["head", "controls", "plot2", "save"])


def _qc_metrics_template(has_areas: bool) -> str:
    """The *shipped* metrics text, for tests to pin without reading overrides."""
    return builtin_assemble(METRICS_TEMPLATE_ID, _metrics_blocks(has_areas))


def _qc_filter_template(filter_cells: bool, filter_genes: bool) -> str:
    """The *shipped* filter text, for tests to pin without reading overrides."""
    blocks = qc_filter_preview(1 if filter_cells else None,
                               1 if filter_genes else None).blocks
    return builtin_assemble(FILTER_TEMPLATE_ID, blocks)


# ── Keep/drop arithmetic ─────────────────────────────────────────────────────

def counts_per_cell(adata) -> np.ndarray:
    """Transcripts per cell, as a float array."""
    X = adata.X
    total = X.sum(axis=1)
    return np.asarray(total, dtype=np.float64).ravel()


def cells_per_gene(adata, cell_mask=None) -> np.ndarray:
    """How many of the *kept* cells detect each gene.

    Counted without slicing the matrix. A boolean row-slice of a 500k-row CSR is
    a full copy of it; repeating the mask over ``indptr`` and binning ``indices``
    reads the same answer out of the arrays that are already there.
    """
    import scipy.sparse as sp

    X = adata.X
    if cell_mask is None:
        cell_mask = np.ones(adata.n_obs, dtype=bool)
    if sp.issparse(X):
        X = X.tocsr()
        rows = np.repeat(cell_mask, np.diff(X.indptr))
        nz = X.data != 0
        return np.bincount(X.indices[rows & nz], minlength=X.shape[1])
    dense = np.asarray(X)
    return (dense[cell_mask] > 0).sum(axis=0).astype(np.int64)


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    status = StatusProxy(ctx.viewer)

    # ── Widgets ──────────────────────────────────────────────────────────
    # Labels read as English; the tooltip names the template parameter, the
    # convention the Clustering tab established.
    cells_check = CheckBox(
        label="Filter cells", value=True,
        tooltip="Drop cells with too few transcripts to cluster reliably.",
    )
    counts_spin = SpinBox(
        label="Min transcripts per cell", min=1, max=100000,
        value=DEFAULT_MIN_COUNTS,
        tooltip="Cells with fewer total transcripts than this are dropped.\n"
                "10 is the usual Xenium starting point.\n\n"
                "Template parameter: min_counts",
    )
    genes_check = CheckBox(
        label="Filter genes", value=True,
        tooltip="Drop genes detected in too few cells to be informative.",
    )
    cells_spin = SpinBox(
        label="Min cells per gene", min=1, max=1000000,
        value=DEFAULT_MIN_CELLS,
        tooltip="Genes detected in fewer cells than this are dropped,\n"
                "counted over the cells that survive the cell filter.\n\n"
                "Template parameter: min_cells",
    )

    readout = QLabel("Keep: —")
    readout.setWordWrap(True)

    metrics_button = PushButton(label="QC metrics && plots")
    apply_button = PushButton(label="Apply filter")
    revert_button = PushButton(label="Revert to all cells")
    progress = make_progress_bar()

    summary = QTextEdit()
    summary.setReadOnly(True)
    summary.setFontFamily("monospace")
    summary.setMaximumHeight(140)

    def _full():
        """The unfiltered table --- what a cutoff is judged against.

        ``is not None`` rather than ``or``: AnnData defines ``__len__``, so a
        table with no cells would test falsy and silently fall through.
        """
        full = getattr(ctx, "full_adata", None)
        return full if full is not None else ctx.adata

    # ── The keep/drop readout ────────────────────────────────────────────
    # Two sorted vectors, computed once per bound table and then answered with
    # searchsorted: exact, and fast enough to run on every spin-box tick without
    # a debounce. Only the per-gene half depends on the cell cutoff, so only it
    # is recomputed when that moves --- hence the two cache slots.
    cache = {"src": None, "counts_sorted": None, "min_counts": None,
             "per_gene_sorted": None}

    def _ensure_counts():
        adata = _full()
        if adata is None:
            return False
        if cache["src"] != id(adata):
            cache.update(src=id(adata), min_counts=None, per_gene_sorted=None,
                         counts_sorted=np.sort(counts_per_cell(adata)))
        return True

    def _ensure_per_gene(min_counts):
        adata = _full()
        if cache["min_counts"] == min_counts and cache["per_gene_sorted"] is not None:
            return
        mask = (counts_per_cell(adata) >= min_counts if min_counts is not None
                else np.ones(adata.n_obs, dtype=bool))
        cache["per_gene_sorted"] = np.sort(cells_per_gene(adata, mask))
        cache["min_counts"] = min_counts

    def _would_keep():
        """``(cells_kept, cells_total, genes_kept, genes_total)`` or ``None``."""
        if not _ensure_counts():
            return None
        adata = _full()
        counts = cache["counts_sorted"]
        min_counts = int(counts_spin.value) if cells_check.value else None
        n_cells = (len(counts) - int(np.searchsorted(counts, min_counts))
                   if min_counts is not None else len(counts))
        _ensure_per_gene(min_counts)
        per_gene = cache["per_gene_sorted"]
        min_cells = int(cells_spin.value) if genes_check.value else None
        n_genes = (len(per_gene) - int(np.searchsorted(per_gene, min_cells))
                   if min_cells is not None else adata.n_vars)
        return n_cells, len(counts), n_genes, adata.n_vars

    def _refresh_readout():
        counts_spin.enabled = cells_check.value
        cells_spin.enabled = genes_check.value
        apply_button.enabled = bool(cells_check.value or genes_check.value)
        revert_button.enabled = bool(state.get("qc_filter"))
        try:
            kept = _would_keep()
        except Exception:            # noqa: BLE001 - a readout must not raise
            kept = None
        if kept is None:
            readout.setText("Keep: —")
            return
        n_cells, all_cells, n_genes, all_genes = kept
        pct = 100.0 * n_cells / all_cells if all_cells else 0.0
        readout.setText(
            f"Keep {n_cells:,} of {all_cells:,} cells ({pct:.1f}%) "
            f"· {n_genes:,} of {all_genes:,} genes"
        )

    for widget in (cells_check, counts_spin, genes_check, cells_spin):
        widget.changed.connect(lambda *_: _refresh_readout())

    # ── Providers ────────────────────────────────────────────────────────
    def _qc_preview() -> Preview:
        """What "Apply filter" would run with the boxes as they stand.

        Block selection is delegated to ``_helpers.qc_filter_preview`` so this
        and the launch-time restore in ``app.py`` turn the same two numbers into
        the same blocks --- two expressions of that would drift.

        With neither box ticked there is nothing to render: no assembly filters
        nothing, because that is not a step. The pane shows what ticking both
        would run, and the note says so.
        """
        min_counts = coerce(counts_spin.value) if cells_check.value else None
        min_cells = coerce(cells_spin.value) if genes_check.value else None
        if min_counts is None and min_cells is None:
            return qc_filter_preview(
                coerce(counts_spin.value), coerce(cells_spin.value),
            )._replace(
                note="nothing is filtered until you tick a filter and press Apply"
            )
        return qc_filter_preview(min_counts, min_cells)

    def _has_areas():
        adata = _full()
        return adata is not None and all(c in adata.obs.columns for c in AREA_COLUMNS)

    def _qc_metrics_preview() -> Preview:
        """What "QC metrics & plots" would draw for the table now bound."""
        adata = _full()
        n_vars = int(adata.n_vars) if adata is not None else 0
        return Preview(
            _metrics_blocks(_has_areas()),
            {
                "title": f"QC — {getattr(ctx.data_path, 'name', ctx.data_path)}",
                "paths": list(ctx.plot_paths("qc_metrics")),
                "percent_top": percent_top_for(n_vars),
            },
        )

    ctx.state.setdefault("template_preview", {})[FILTER_TEMPLATE_ID] = _qc_preview
    ctx.state.setdefault("template_preview", {})[METRICS_TEMPLATE_ID] = _qc_metrics_preview

    # ── QC metrics & plots ───────────────────────────────────────────────
    def _on_metrics():
        blocks, params, _ = _qc_metrics_preview()
        missing = [c for c in ("control_probe_counts", "control_codeword_counts")
                   if c not in _full().obs.columns]
        if missing:
            status.value = f"QC metrics need obs columns {missing}"
            summary.setPlainText(
                "This table has no negative-control columns "
                f"({', '.join(missing)}), so the control rates cannot be "
                "computed. They come with the 10x output; a table built from a "
                "custom segmentation may not carry them."
            )
            return
        step = Step(
            id="plot:qc_metrics",
            **_resolved(METRICS_TEMPLATE_ID, blocks),
            params=params,
            deps=[ctx.cell_root()],
            kind=TERMINAL,
            label="QC metrics",
            outputs=["fig", "cprobes", "cwords"],
        )
        gen = ctx.dataset_generation
        metrics_button.enabled = False

        @thread_worker
        def _run():
            return ctx.run_step(step)

        def _done(out):
            metrics_button.enabled = True
            if ctx.dataset_generation != gen:
                return
            # save=False: the template's own `save` block already wrote them.
            ctx.show_plot(out["fig"], "qc_metrics", title="QC metrics",
                          save=False, paths=params["paths"])
            summary.setPlainText(
                f"Negative DNA probe count % : {out['cprobes']:.4f}\n"
                f"Negative decoding count %  : {out['cwords']:.4f}"
            )
            status.value = "QC metrics drawn"

        worker = _run()
        worker.returned.connect(_done)
        worker.errored.connect(lambda e: (setattr(metrics_button, "enabled", True),
                                          setattr(status, "value", f"QC metrics failed: {e}")))
        timer, _ = attach_spinner(worker, lambda m: setattr(status, "value", m),
                                  "Computing QC metrics...", progress_bar=progress)
        state["_qc_spinner"] = timer  # prevent GC
        worker.start()

    metrics_button.clicked.connect(_on_metrics)

    # ── Apply / Revert ───────────────────────────────────────────────────
    def _describe():
        adata, full = ctx.adata, _full()
        if not state.get("qc_filter"):
            return f"QC: off — all {full.n_obs:,} cells, {full.n_vars:,} genes."
        qc = state["qc_filter"]
        return (
            f"QC applied: {adata.n_obs:,} of {full.n_obs:,} cells, "
            f"{adata.n_vars:,} of {full.n_vars:,} genes.\n"
            f"  min_counts (per cell) : {qc.get('min_counts')}\n"
            f"  min_cells  (per gene) : {qc.get('min_cells')}\n"
            "Cells the filter dropped are no longer coloured in the image."
        )

    def _run_in_worker(fn, busy, done_msg):
        gen = ctx.dataset_generation
        apply_button.enabled = revert_button.enabled = False

        @thread_worker
        def _work():
            fn()

        def _done(_=None):
            if ctx.dataset_generation != gen:
                return
            cache["src"] = None          # the bound table changed under us
            _refresh_readout()
            summary.setPlainText(_describe())
            status.value = done_msg

        def _failed(exc):
            apply_button.enabled = revert_button.enabled = True
            status.value = f"{busy.rstrip('.')} failed: {exc}"

        worker = _work()
        worker.returned.connect(_done)
        worker.errored.connect(_failed)
        timer, _ = attach_spinner(worker, lambda m: setattr(status, "value", m),
                                  busy, progress_bar=progress)
        state["_qc_spinner"] = timer
        worker.start()

    def _on_apply():
        preview = _qc_preview()
        _run_in_worker(lambda: ctx.apply_qc_filter(preview),
                       "Applying QC filter...", "QC filter applied")

    def _on_revert():
        _run_in_worker(ctx.clear_qc_filter,
                       "Reverting to all cells...", "QC filter removed")

    apply_button.clicked.connect(_on_apply)
    revert_button.clicked.connect(_on_revert)

    # ── restore_session ──────────────────────────────────────────────────
    def _restore_session(session):
        """Display only. ``app.ensure_qc_filter()`` has already applied it.

        Deliberately not the thing that applies the filter: this handler runs
        after every other tab's, and one of them reaching ``ensure_normalized``
        first would normalise the unfiltered table.
        """
        qc = session.get("qc_filter") or state.get("qc_filter")
        if qc:
            if qc.get("min_counts") is not None:
                counts_spin.value = int(qc["min_counts"])
            if qc.get("min_cells") is not None:
                cells_spin.value = int(qc["min_cells"])
            cells_check.value = qc.get("min_counts") is not None
            genes_check.value = qc.get("min_cells") is not None
        cache["src"] = None
        _refresh_readout()
        summary.setPlainText(_describe())

    # The Clustering tab (and anyone else) wants to know when the cell set moves.
    # Bounced to the GUI thread: rebind_cells fires these, and Apply runs in a
    # napari worker -- so the callback that reaches a widget must not be the one
    # the worker calls directly.
    @ensure_main_thread
    def _on_cells_rebound():
        cache["src"] = None
        _refresh_readout()

    state.setdefault("qc_listeners", []).append(_on_cells_rebound)

    widget = make_tab(
        cells_check, counts_spin,
        genes_check, cells_spin,
        readout,
        toolbar_row(metrics_button.native, apply_button.native, revert_button.native),
        progress,
        summary,
    )

    # Deferred past the constructor: the first readout reads the whole count
    # matrix, and a launch should not pay for a tab nobody has opened yet.
    QTimer.singleShot(0, _refresh_readout)

    return widget, {"restore_session": _restore_session}
