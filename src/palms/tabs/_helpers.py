"""Shared helpers used across tab modules."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
from qtpy.QtWidgets import QWidget, QVBoxLayout
from superqt.utils import ensure_main_thread

from palms.utils.prov_graph import (
    ProvGraph, CycleError, SETUP, ARTIFACT, TERMINAL,
)
from palms.utils.environment import environment_code, same_environment
from palms.utils.fig_render import to_figure
from palms.utils.plot_output import (
    plot_formats, recorded_paths, save_figure, save_paths,
)
from palms.utils.reporting import get_logger, report_recording_failure
from palms.utils.step_templates import (
    Preview, builtin_assemble, builtin_spec, builtin_text, check_base_namespace,
    step_template as _resolved,
)
from palms.utils.steps import Step, StepExecutor, coerce

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext

log = get_logger(__name__)

# The provenance graph, written beside the store on every recorded step. The
# zarr session attr is only written by ``save_session`` (dataset switch / exit),
# so this is what a mid-session reader — the verification script, the next
# launch after a crash — should believe.
PROV_GRAPH_SIDECAR = "prov_graph.json"


# Template text now lives in ``utils/step_templates/builtin/*.tmpl``. These
# bindings go through ``builtin_text``, which reads only the shipped files and
# never consults an override path — so the tests that pin template text stay
# immune to a developer's own customisations by construction rather than by
# remembering to isolate them.
#
# Binds ``adata_norm`` rather than mutating ``adata``, so a step that makes its
# own copy cannot end up normalising twice. Mirrors what the viewer has always
# actually run (``utils.gene_analysis.get_normalized_adata``) — including the
# ``target_sum`` the old recorded cell omitted.
#
# ``builtin_assemble``, not ``builtin_text``: the template gained a second
# scaling block when the target became settable, and ``builtin_text`` returns
# *every* block in file order — which for this one means two
# ``normalize_total`` calls. This constant is the default assembly, the fixed
# target the viewer has always used.
_NORMALIZE_BLOCKS = ("copy", "scale.fixed", "tail")
_NORMALIZE_TEMPLATE = builtin_assemble("normalize", list(_NORMALIZE_BLOCKS))


# Built on ``adata_norm`` rather than ``adata``: every consumer of the spatial
# graph (nhood enrichment, co-occurrence, ligrec) works on the normalised copy,
# and squidpy stores the graph in ``.obsp`` of whichever object it was given.
# The old recorded cell built it on ``adata`` — so the notebook's consumers read
# a graph that did not exist on the object they were passed.
_SPATIAL_NEIGHBORS_TEMPLATE = builtin_text("spatial_neighbors")


_ANNOT_POLYGONS_TEMPLATE = builtin_text("annot.polygons")


def annotation_polygons_preview(ctx) -> Preview:
    """The ``annot.polygons`` step as the annotation layer stands right now.

    Shared rather than per-tab: the annotation-neighbourhood and
    annotation-distance tabs both depend on the same drawn shapes, and two
    expressions of "what has been drawn" would drift the moment one of them
    learned something the other did not.

    Shapes with no type assigned are left out. They are not annotations of type
    ``""`` — the Annotations tab has simply not been told what they are yet, and
    inlining them would put an unnamed group in the notebook's results.
    """
    layer = getattr(ctx, "annotation_layer", None)
    shapes = list(layer.data) if layer is not None else []
    raw = list(layer.properties.get("annotation_type", [])) if layer is not None else []
    raw += [""] * (len(shapes) - len(raw))
    typed = [(s, str(t)) for s, t in zip(shapes, raw) if str(t).strip()]
    return Preview(
        list(builtin_spec("annot.polygons").blocks),
        {
            "polygons": [np.round(np.asarray(s), 2).tolist() for s, _ in typed],
            "types": [t for _, t in typed],
            "pixel_size": float(ctx.pixel_size),
        },
    )


# ── Normalisation ────────────────────────────────────────────────────────────

#: Counts per cell after normalisation, when a fixed target is used. scanpy's
#: own default is the median count across cells; this is the other convention,
#: and the one the viewer has always used.
DEFAULT_TARGET_SUM = 1e4


def normalize_preview(target_sum) -> Preview:
    """The ``normalize`` step for a scaling target; ``None`` means the median.

    Module level, beside :func:`qc_filter_preview`, and for the same reason: the
    Preprocess tab's provider and :func:`ViewerContext.ensure_normalized` both
    have to turn one setting into the same blocks, and two expressions of that
    would drift — which is precisely the defect that made this step a template
    in the first place (the GUI scaled to 1e4 while the recorded cell used
    scanpy's median default).
    """
    if target_sum is None:
        return Preview(["copy", "scale.median", "tail"], {})
    return Preview(["copy", "scale.fixed", "tail"],
                   {"target_sum": coerce(float(target_sum))})


# ── QC filtering ─────────────────────────────────────────────────────────────

#: The node that filters the cell set, when one is in force.
QC_NODE_ID = "qc_filter"

def qc_filter_preview(min_counts, min_cells) -> Preview:
    """The ``qc.filter`` step for a pair of cutoffs; ``None`` switches one off.

    Module level, beside :func:`annotation_polygons_preview`, and for the same
    reason: the QC tab's provider and the launch-time restore in ``app.py`` both
    have to turn the same two numbers into the same blocks, and two expressions
    of that would drift.

    There is no assembly that filters nothing, because that is not a step. A
    caller with both cutoffs off wants :func:`ViewerContext.clear_qc_filter`.
    """
    blocks = [name for name, on in (("cells", min_counts is not None),
                                    ("genes", min_cells is not None)) if on]
    blocks.append("bind")
    params = {}
    if min_counts is not None:
        params["min_counts"] = int(min_counts)
    if min_cells is not None:
        params["min_cells"] = int(min_cells)
    return Preview(blocks, params)


def qc_label(min_counts, min_cells) -> str:
    """The node label for a pair of cutoffs, as the Notebook tab shows it."""
    halves = []
    if min_counts is not None:
        halves.append(f"≥{int(min_counts)} counts/cell")
    if min_cells is not None:
        halves.append(f"≥{int(min_cells)} cells/gene")
    return "QC filter: " + ", ".join(halves)


# ── magicgui ComboBox default helper ─────────────────────────────────────────

def combo_value_kwargs(choices, index: int = 0) -> dict:
    """Return kwargs for a magicgui ComboBox 'value' that won't crash on empty/short choices.

    magicgui raises ``ValueError: None is not a valid choice`` if you pass
    ``value=None`` (or any non-member) against its ``choices`` list, which
    happens for datasets that legitimately have no clusterings (e.g. a Crop
    Dataset export, which has no ``analysis/`` folder). Returns
    ``{"value": choices[index]}`` only when that element exists, else ``{}`` —
    with no ``value`` kwarg, magicgui defaults to the first choice, or ``None``
    when choices is empty, neither of which raises.
    """
    seq = list(choices)
    if len(seq) > index:
        return {"value": seq[index]}
    return {}


# ── Tab layout helper ────────────────────────────────────────────────────────

def scrollable(widget: QWidget) -> QWidget:
    """Wrap a widget so it cannot impose its minimum size on the dock.

    The control panel is a ``QTabWidget`` of ``QTabWidget``s, and a stacked
    widget's minimum is the maximum over *all* its pages, hidden ones included.
    So a single page that reports a wide or tall minimum becomes the floor for
    the whole Controls dock, and the separator between the dock and the
    canvas stops moving. A ``QScrollArea`` reports a fixed ~68px minimum
    whatever it contains, which is what keeps that floor low: every page must
    go through here (or through `make_tab`, which does).
    """
    from qtpy.QtWidgets import QScrollArea
    scroll = QScrollArea()
    scroll.setWidget(widget)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QScrollArea.Shape.NoFrame)
    return scroll


def toolbar_row(*buttons) -> QWidget:
    """Pack buttons into a fixed-height row that scrolls sideways when narrow.

    A plain ``QHBoxLayout`` of buttons cannot shrink below the sum of their
    label widths, so a button bar outside a scroll area sets the dock's minimum
    width all by itself (the six-button Notebook toolbar cost 528px). Keeping
    the row pinned but horizontally scrollable costs ~68px instead, and unlike
    wrapping the whole tab it leaves the bar visible while the content below
    scrolls.
    """
    from qtpy.QtCore import Qt
    from qtpy.QtWidgets import QHBoxLayout

    row = QWidget()
    layout = QHBoxLayout()
    layout.setContentsMargins(0, 0, 0, 0)
    for button in buttons:
        layout.addWidget(button)
    layout.addStretch()
    row.setLayout(layout)

    scroll = scrollable(row)
    scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    scroll.setFixedHeight(row.sizeHint().height() + 2)
    return scroll


def labelled(w) -> QWidget:
    """Return the Qt widget for a magicgui widget, *with* its label attached.

    ``w.native`` is the bare control. magicgui keeps a widget's caption in a
    separate ``_LabeledWidget`` wrapper that only a ``Container`` ever creates,
    so dropping ``.native`` into a plain layout silently discards it — which is
    how every ``Slider(label="n_neighbors", ...)`` in this package rendered as
    an anonymous slider reading ``15``. Found by screenshotting the running
    viewer over ``--mcp``; no test caught it, because none of them render a tab.

    A one-widget ``Container`` restores the label through public API only. The
    private ``_LabeledWidget`` renders identically (measured: same 142px
    minimum, same 22px height) but would break on a magicgui rename.

    Two widgets are returned bare, and both cases are load-bearing:

    - **No label.** A widget built without ``label=`` reports ``''``, never a
      name derived from anything, so this cannot invent a caption for a control
      that never had one.
    - **A ``ButtonWidget``** (``CheckBox``, ``PushButton``, ``RadioButton``).
      Qt paints their text *on* the control, so a caption would show it twice.
      magicgui's own ``Container`` skips them for the same reason — this mirrors
      upstream rather than inventing a rule.
    """
    from magicgui.widgets import Container
    from magicgui.widgets.bases import ButtonWidget

    if isinstance(w, (Container, ButtonWidget)) or not w.label:
        return w.native
    box = Container(widgets=[w], labels=True)
    box.margins = (0, 0, 0, 0)   # the row must sit flush like a bare .native did
    return box.native


def make_tab(*widgets_and_natives) -> QWidget:
    """Pack magicgui widgets and raw QWidgets into a scrollable container."""
    inner = QWidget()
    layout = QVBoxLayout()
    layout.setContentsMargins(4, 4, 4, 4)
    for w in widgets_and_natives:
        if hasattr(w, "native"):
            layout.addWidget(labelled(w))
        else:
            layout.addWidget(w)
    layout.addStretch()
    inner.setLayout(layout)
    return scrollable(inner)


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

def make_progress_dialog(title: str, parent=None):
    """Return (dialog, progress_bar, label) — caller shows and manages them."""
    from qtpy.QtWidgets import QDialog, QVBoxLayout as _QVBoxLayout, QLabel, QProgressBar
    from qtpy.QtCore import Qt
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.setWindowFlags(dlg.windowFlags() & ~Qt.WindowType.WindowCloseButtonHint)
    dlg.setMinimumWidth(500)
    layout = _QVBoxLayout(dlg)
    lbl = QLabel("Starting…")
    lbl.setWordWrap(True)
    bar = QProgressBar()
    bar.setRange(0, 100)
    layout.addWidget(lbl)
    layout.addWidget(bar)
    return dlg, bar, lbl


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


def attach_spinner(worker, set_status_fn, initial_msg: str, progress_bar=None):
    """Animate a spinner in the status bar while *worker* runs.

    Returns ``(timer, update_msg_fn)`` where ``update_msg_fn(msg)`` changes the
    animated text (connect to ``worker.yielded`` for stage messages).
    If *progress_bar* is given it is shown immediately and hidden on finish.
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
    if progress_bar is not None:
        progress_bar.setVisible(True)
        worker.finished.connect(lambda: progress_bar.setVisible(False))
    return timer, update_msg


def attach_tqdm_progress(worker, set_status_fn, base_msg: str = "", progress_bar=None):
    """Wire a ProgressMailbox + QTimer to relay tqdm updates to the status bar.

    Returns a ``post_fn`` callable safe to call from the background thread;
    pass it into :func:`qt_tqdm_context` inside the worker.
    If *progress_bar* is given it is shown immediately and hidden on finish.
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
    if progress_bar is not None:
        progress_bar.setVisible(True)
        worker.finished.connect(lambda: progress_bar.setVisible(False))
    return mailbox.post, timer   # caller must hold timer reference to prevent GC


def make_progress_bar():
    """Return a hidden indeterminate QProgressBar for embedding in tab layouts."""
    from qtpy.QtWidgets import QProgressBar
    bar = QProgressBar()
    bar.setRange(0, 0)       # indeterminate / marquee animation
    bar.setMaximumHeight(16)
    bar.setVisible(False)
    return bar


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

    import tqdm.std as _tqdm_std
    original_std = _tqdm_std.tqdm

    _tqdm_module.tqdm = _StatusTqdm
    _tqdm_auto.tqdm = _StatusTqdm
    _tqdm_std.tqdm = _StatusTqdm
    try:
        yield
    finally:
        _tqdm_module.tqdm = original_tqdm
        _tqdm_auto.tqdm = original_auto
        _tqdm_std.tqdm = original_std


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

    # ── code recording (provenance DAG) ──────────────────────────────────
    # The graph in state["prov_graph"] is the source of truth. record_node
    # upserts a step keyed by a stable artifact id; re-recording revises it in
    # place and flags dependents stale. The flat state["code_journal"] and the
    # on-disk .py file are kept as a derived, append-style view so the existing
    # Notebook tab and Save-code action keep working during the tab-by-tab
    # migration (task 5 switches those consumers to the derived, topo-ordered
    # graph output; task 4 makes the file persistent across sessions).
    state.setdefault("prov_graph", ProvGraph())

    def _write_code_file():
        code_path = ctx.data_path / state.get("code_file", "code.py")
        try:
            with open(code_path, 'w') as f:
                f.write("\n".join(state["code_journal"]) + "\n")
        except (PermissionError, OSError) as e:
            from palms.utils.adata_persistence import _maybe_show_permission_dialog
            _maybe_show_permission_dialog(e, "code journal")

    def _emit_flat(code: str):
        """Append to the derived flat journal + file + Notebook tab."""
        state["code_journal"].append(code)
        sync_fn = state.get("_notebook_sync_fn")
        if sync_fn:
            sync_fn()
        _write_code_file()
        _save_prov_graph()

    def _save_prov_graph():
        """Persist the provenance graph as soon as it changes.

        The graph used to reach disk only inside ``save_session``, which runs on
        a dataset switch or at viewer exit — while the *artifacts* it explains
        (``clustering_*`` columns, ``uns['rank_genes_*']``) are persisted
        the moment they are produced. A session verified, inspected, or killed
        in between therefore had results on disk whose code was nowhere: measured
        on a real session, the store held a 16-minute-old three-node graph while
        the table already carried two Leiden clusterings and a rank-genes result.

        Written as a sidecar rather than into the store: it costs one small
        atomic file write per recorded step, where updating the zarr group means
        copying every parquet under ``viewer_session/``. ``save_session`` still
        writes the attr, and the sidecar takes precedence on load.

        Gated on ``prov_graph_restored``, which ``app.py`` sets once the session
        has been restored. Tabs seed a preamble node while the viewer is still
        being built — writing at that point replaced a 13-node graph on disk
        with a one-node stub, which the next launch then preferred over the
        session attr, and the DAG came up empty.
        """
        graph = state.get("prov_graph")
        if not state.get("prov_graph_restored"):
            return
        if graph is None or not len(graph) or ctx.data_path is None:
            return
        try:
            from palms.utils.adata_persistence import sidecar_write_path
            from palms.utils.zarr_safe import atomic_json
            atomic_json(sidecar_write_path(ctx, PROV_GRAPH_SIDECAR), graph.to_list())
        except (OSError, TypeError, ValueError) as e:
            # Never let a persistence hiccup abort the user's analysis action;
            # save_session remains the backstop.
            log.warning("could not persist the provenance graph: %s", e)

    ctx.save_prov_graph = _save_prov_graph

    def _record_node(node_id: str, code: str, deps=(), kind: str = ARTIFACT,
                     label: str = None, params: dict = None):
        """Record one analysis step as a node in the provenance graph.

        ``node_id`` is the stable identity of the artifact produced (the upsert
        key); all parameters belong in ``code``/``params``, never in the id, so
        re-running the same artifact revises its node rather than duplicating.
        Dependencies must already be recorded — a missing dep is surfaced here,
        at record time, instead of as a NameError at replay.
        """
        if not state.get("record_code"):
            return
        graph = state.setdefault("prov_graph", ProvGraph())
        prev = graph.get(node_id)
        prev_code = prev.code if prev is not None else None
        try:
            graph.upsert(node_id, code, deps=deps, kind=kind,
                         label=label, params=params)
        except (KeyError, CycleError) as e:
            # Never let a recorder bug abort the user's analysis action; degrade
            # to appending the snippet so the code isn't silently lost — but say
            # so, because from here on the notebook is incomplete.
            report_recording_failure(node_id, e)
            _emit_flat(code)
            return
        if prev_code != code:  # new node, or revised → show the (new) code
            _emit_flat(code)

    ctx.record_node = _record_node

    # ── run_step (the executed-code-is-recorded-code path) ───────────────
    @ensure_main_thread
    def _emit_flat_main_thread(code: str):
        _emit_flat(code)

    def _get_executor():
        """Lazily build the StepExecutor, sharing the session's ProvGraph.

        The namespace mirrors the exported notebook's globals, seeded with the
        objects the viewer already loaded. Note the one documented divergence:
        the ``preamble`` node records ``xenium(data_path)``, but the viewer
        reaches the same objects via the zarr cache rather than re-reading the
        raw output. Every *other* step executes exactly what it records.
        """
        if ctx.executor is None:
            import matplotlib.pyplot as plt
            import numpy as _np
            import pandas as pd
            import scanpy as sc
            import spatialdata as sd
            import squidpy as sq
            from pathlib import Path as _Path

            # A placeholder graph, deliberately not state["prov_graph"]. The
            # executor is also built for ``preview_step``, and a preview must
            # not have side effects — seeding a session graph (or creating one)
            # from a call that only draws a pane is exactly the kind of side
            # effect that rule exists to forbid. ``_run_step`` points it at the
            # real graph on every run, so nothing here needs to be right yet.
            graph = ProvGraph()
            base_namespace = {
                "sc": sc, "sq": sq, "sd": sd, "pd": pd, "np": _np, "plt": plt,
                "Path": _Path,
                "data_path": ctx.data_path,
                "sdata": ctx.sdata,
                "adata": ctx.adata,
            }
            # The declared set is what template validation checks against, so a
            # name added here and not there (or vice versa) is caught now rather
            # than as a NameError when someone replays the notebook.
            check_base_namespace(base_namespace)
            ctx.executor = StepExecutor(namespace=base_namespace, graph=graph)
        return ctx.executor

    def _step_graph():
        """The graph ``run_step`` records into — and its root, guaranteed.

        Two things this settles, each of which was a real failure.

        **Which graph.** ``StepExecutor.run`` always records; recording is not a
        mode of ``run``, it is what ``run`` is (see :meth:`StepExecutor.preview`
        for why that is deliberate). So Preferences -> "Record reproducible code"
        off cannot mean *don't record*, only *don't record here*: the executor is
        pointed at a throwaway graph that nothing exports and nothing persists,
        which leaves the session's own graph exactly as the user left it.

        Letting the session graph take those nodes instead would be worse than
        the crash it replaces. ``_record_node`` still returns early with the
        toggle off, so ``environment``, every ``clustering:<key>`` and every
        ``plot:*`` terminal would be missing while the migrated steps were
        present — and ``save_session`` writes that graph at exit. A persisted
        graph describing an analysis nobody ran that way is the shape of
        provenance defect this codebase is most careful about.

        **The root.** Every migrated step declares ``deps=["preamble"]`` and
        ``ProvGraph.upsert`` refuses an unknown dep, so a graph without a
        preamble makes ``run_step`` raise *after* the step's code has already
        executed — the user sees a failure for a computation that succeeded, and
        its output binding is thrown away. ``app.py`` seeds the preamble at
        launch only when recording is on, which is precisely the case that
        cannot be relied on here, so ``run_step`` establishes its own
        precondition rather than inheriting a conditional launch-time side
        effect. That also self-heals the recording-on path if the seed is ever
        missed again.
        """
        if state.get("record_code"):
            graph = state.setdefault("prov_graph", ProvGraph())
            if graph.get("preamble") is None:
                # Through the real recorder, so this seed is indistinguishable
                # from app.py's: it brings the environment node with it and
                # reaches the flat journal.
                _record_preamble()
            return graph

        graph = state.get("_unrecorded_prov_graph")
        if graph is None:
            graph = state["_unrecorded_prov_graph"] = ProvGraph()
        if graph.get("preamble") is None:
            # Not through _record_preamble: that goes via _record_node, which
            # returns early with recording off — which is the whole reason this
            # graph has no root. Upsert the same text directly.
            graph.upsert("preamble", _preamble_code(), kind=SETUP,
                         label="Setup & data loading")
        return graph

    def _run_step(step: Step, progress=None) -> dict:
        """Execute *step* and record the same source, then sync viewer state.

        Returns the step's declared outputs. Raises ``StepError`` on failure —
        callers surface it rather than swallowing it, so a broken step is
        visible instead of producing a node for an artifact that never existed.
        """
        executor = _get_executor()
        # Keep the executor pointed at the live graph as well as the live
        # objects. Both can be swapped out from under it: a dataset reload
        # rebinds ctx.adata/ctx.sdata, and toggling recording replaces
        # state["prov_graph"] wholesale.
        executor.graph = _step_graph()
        # Keep the namespace pointed at the live objects: a dataset reload
        # rebinds ctx.adata/ctx.sdata without going through the executor.
        if executor.ns.get("adata") is not ctx.adata:
            executor.ns["adata"] = ctx.adata
        if executor.ns.get("sdata") is not ctx.sdata:
            executor.ns["sdata"] = ctx.sdata

        outputs = executor.run(step, progress=progress)

        # A step may rebind rather than mutate (``adata = adata[:, mask]``);
        # follow it so the GUI and the notebook stay on the same object.
        if executor.ns.get("adata") is not ctx.adata:
            ctx.adata = executor.ns["adata"]

        if state.get("record_code"):
            # Steps run in napari worker threads; the flat journal writes a file
            # and refreshes the Notebook tab widget, so bounce it to the GUI thread.
            _emit_flat_main_thread(step.render())
        return outputs

    ctx.run_step = _run_step

    # ── preview_step (display only, records nothing) ─────────────────────
    def _preview_step(step: Step, bindings: dict = None) -> dict:
        """Execute *step* for display only. Nothing is recorded, nothing persists.

        Note what is deliberately missing relative to :func:`_run_step`: no
        ``_emit_flat_main_thread``, no graph write, and no ``ctx.adata``
        follow-up — a preview executes in a scratch namespace, so there is
        nothing to follow. See :meth:`StepExecutor.preview` for why the scratch
        namespace is a correctness requirement rather than tidiness.
        """
        executor = _get_executor()
        if executor.ns.get("adata") is not ctx.adata:
            executor.ns["adata"] = ctx.adata
        if executor.ns.get("sdata") is not ctx.sdata:
            executor.ns["sdata"] = ctx.sdata
        return executor.preview(step, bindings=bindings)

    ctx.preview_step = _preview_step

    # ── record_code (backward-compat shim) ───────────────────────────────
    def _record_code(code: str, tag: str = None):
        """Legacy string-append recorder, mapped onto the provenance graph.

        A ``tag`` becomes the node id (so the old tag-dedup becomes upsert);
        untagged calls get a fresh opaque id so they still append. Superseded by
        :func:`record_node` — call sites are migrated to declare real deps.
        """
        if not state.get("record_code"):
            return
        if tag:
            node_id = tag
        else:
            state["_legacy_counter"] = state.get("_legacy_counter", 0) + 1
            node_id = f"_legacy:{state['_legacy_counter']}"
        _record_node(node_id, code, deps=(), kind=ARTIFACT)

    ctx.record_code = _record_code

    # ── record_environment ───────────────────────────────────────────────
    def _record_environment():
        """Record what the analysis is being run with.

        A separate node from ``preamble``, and deliberately one with **no
        dependents**: an environment that differs between recording and replay
        is something to read, not a reason to flag every downstream result
        stale. It has no deps either, so it sorts first among the setup nodes
        (ties break on id, and ``environment`` < ``preamble``).
        """
        code = state.get("_environment_code")
        if code is None:
            code = state["_environment_code"] = environment_code()
        graph = state.get("prov_graph")
        previous = graph.get("environment") if graph is not None else None
        if previous is not None and same_environment(previous.code, code):
            return                       # same versions — keep the original stamp
        _record_node("environment", code, kind=SETUP,
                     label="Environment & seeds")

    ctx.record_environment = _record_environment

    # ── record_preamble ──────────────────────────────────────────────────
    def _preamble_code() -> str:
        """Build the preamble cell's source.

        Split out from :func:`_record_preamble` because ``_step_graph`` needs the
        *text* without the recording: a step's graph must have its root before
        ``upsert`` will accept ``deps=["preamble"]``, and with recording off
        ``_record_node`` is a no-op.
        """
        # A Crop Dataset export has no raw 10x output — the zarr store *is* the
        # data — so ``spatialdata_io.xenium()`` cannot read it, and a notebook
        # recorded with that call fails on its very first cell. Branch on
        # ``is_cache_only``, which prefers what the export *declared* in its
        # cache manifest over what the directory listing suggests: an export
        # that also wrote raw-shaped files would otherwise flip this branch back
        # to ``xenium(data_path)`` and silently drop every derived layer the
        # crop carried. Imported lazily: this module is on the tab import path
        # and ``loader`` pulls in spatialdata.
        from palms.loader import is_cache_only

        cache_only = is_cache_only(ctx.data_path)
        if not cache_only:
            load = (
                "from spatialdata_io import xenium\n"
                f"data_path = Path(r\"{ctx.data_path}\")\n"
                "sdata = xenium(data_path)\n"
            )
        else:
            # Derived from data_path, not recorded as an absolute cache path, so
            # palms-rename-dataset's single ``data_path = Path(r"…")`` rewrite
            # keeps moving this notebook with its dataset.
            load = (
                "# This dataset is a Crop Dataset export: it has no raw 10x output,\n"
                "# so the zarr store written by the crop is the thing to read.\n"
                f"data_path = Path(r\"{ctx.data_path}\")\n"
                "sdata = sd.read_zarr(data_path / \"sdata_cached.zarr\")\n"
            )

        # Which *cells* the analysis is about. Tools → Segmentation swaps the
        # Xenium segmentation for a custom one, rebinding ``ctx.adata`` to a
        # different set of cells — and every node recorded after that swap is
        # about those cells while the preamble still said ``sdata["table"]``.
        # A replay then reproduced the whole notebook against the Xenium cells
        # and reported no error, which is the worst shape a provenance defect
        # can take: the numbers come out, and they are about something else.
        #
        # The custom table is not in the 10x output — it is built by
        # ``scripts/extract_seurat_segmentation.R`` + ``build_custom_segmentation.py``
        # and cached into the store by the tab, which is the only copy a replay
        # can rely on (the h5ad the user picked may be anywhere, or gone).
        if state.get("segmentation_source") == "custom":
            source = "sdata" if cache_only else (
                "sd.read_zarr(data_path / \"sdata_cached.zarr\")"
            )
            table = (
                "\n# Custom cell segmentation (Tools -> Segmentation), not part of the\n"
                "# 10x output: PALMS caches it in the viewer's store, which is where a\n"
                "# replay has to read it from.\n"
                f"adata = {source}.tables[\"custom_table\"].copy()"
            )
        else:
            table = "adata = sdata[\"table\"].copy()"

        return (
            "import scanpy as sc\n"
            "import squidpy as sq\n"
            "import spatialdata as sd\n"
            "import matplotlib.pyplot as plt\n"
            "import pandas as pd\n"
            "import numpy as np\n"
            "from pathlib import Path\n"
            f"\nplt.rcParams['font.size'] = {state.get('plot_font_size', 10)}\n"
            f"\n# Load data\n"
            + load
            + table
        )

    def _record_preamble():
        _record_environment()
        _record_node("preamble", _preamble_code(), kind=SETUP,
                     label="Setup & data loading")

    ctx.record_preamble = _record_preamble

    # ── QC filter: cell_root, apply, ensure, clear ───────────────────────
    def _cell_root() -> str:
        """The node that says which *cells* a step is about.

        ``"preamble"`` normally, ``"qc_filter"`` once a filter is in force.
        Every step that reads ``obs``, ``var`` or ``X`` roots here rather than
        naming the preamble directly, so the cell set a step was computed on is
        recorded once, at the step -- and stays recorded: a filter applied
        later does not move it. That is what makes the two lineages (the
        unfiltered work, then the filtered work) both true in the graph.

        Steps that read only images (H&E, ARMS, external images, patch
        overlays) or only the dataset path keep depending on ``"preamble"``:
        they are about the dataset, not about the cells.

        Membership is checked rather than assumed, because naming a node that
        is not in the graph raises out of ``upsert``.
        """
        graph = state.get("prov_graph")
        if state.get("qc_filter") and graph is not None and QC_NODE_ID in graph:
            return QC_NODE_ID
        return "preamble"

    ctx.cell_root = _cell_root

    def _cell_scoped_id(base: str) -> str:
        """The id of an *unkeyed* artifact under the current cell set.

        ``normalize``, ``spatial_neighbors`` and ``roi_deg`` have no key of
        their own, so without this a normalisation run under the filter would
        upsert the node the unfiltered clusterings depend on and flag them
        stale -- when nothing about them changed. Suffixing ``:qc`` gives the
        filtered lineage its own node, and a dependent asks for the same
        scoped id rather than writing the base name, so the two can never
        cross. Keyed ids (``clustering:<key>``) do not need it: re-running one
        under the filter is a genuine revision of that result.
        """
        return base if _cell_root() == "preamble" else f"{base}:qc"

    ctx.cell_scoped_id = _cell_scoped_id

    def _apply_qc_filter(preview):
        """Run and record a QC filter, then re-point the viewer at its cells.

        Takes a live ``Preview`` the way ``ensure_annotations`` does, so the tab
        and this function cannot disagree about what the widgets meant.

        The step is a **barrier**: it rebinds ``adata``, which every earlier
        cell-rooted step read, so the notebook has to run all of them before
        it. ``topo_sort`` orders that from the flag alone. Nothing recorded
        before the filter is touched -- those results were computed on the
        full table and still were; only work recorded from here on roots at
        ``qc_filter``.
        """
        from palms.utils.adata_persistence import full_table
        from palms.utils.rebind_cells import (
            kept_mask, rebind_cells, repoint_label_to_obs,
        )

        blocks, params, _ = preview
        if not params:
            raise ValueError("a QC filter needs at least one cutoff")

        full = full_table(ctx)
        full_l2o = ctx.full_label_to_obs
        if full_l2o is None:
            full_l2o = ctx.label_to_obs

        # Always filter *the full table*, never the current one: re-running with
        # a looser cutoff must be able to bring cells back, and a filter of a
        # filter could not. It is also what makes the recorded step the one that
        # replays -- the notebook's adata is the whole table at this point.
        #
        # Restored on failure, because the intermediate state is the dangerous
        # one: ``adata`` would be the full table while ``label_to_obs`` still
        # described the previous filter, and that mismatch paints cells with
        # other cells' values rather than raising.
        previous = ctx.adata
        ctx.adata = full
        try:
            _record_preamble()
            _run_step(Step(
                id=QC_NODE_ID,
                **_resolved("qc.filter", blocks),
                params=params,
                deps=["preamble"],
                kind=SETUP,
                label=qc_label(params.get("min_counts"), params.get("min_cells")),
                outputs=["adata"],
                barrier=True,
            ))
        except Exception:
            ctx.adata = previous
            if ctx.executor is not None:
                ctx.executor.ns["adata"] = previous
            raise
        state["qc_filter"] = {
            "min_counts": params.get("min_counts"),
            "min_cells": params.get("min_cells"),
        }

        rebind_cells(ctx, ctx.adata, repoint_label_to_obs(full_l2o, kept_mask(full, ctx.adata)))
        state["_qc_applied_key"] = (tuple(sorted(params.items())), id(full))
        _save_prov_graph()

    ctx.apply_qc_filter = _apply_qc_filter

    def _ensure_qc_filter():
        """Re-apply the stored cutoffs if they are not already in force.

        Idempotent, and what ``app.py`` calls at launch so a session comes back
        filtered without any tab having to run first.
        """
        from palms.utils.adata_persistence import full_table

        qc = state.get("qc_filter")
        if not qc:
            return
        params = {k: v for k, v in qc.items() if v is not None}
        key = (tuple(sorted(params.items())), id(full_table(ctx)))
        if state.get("_qc_applied_key") == key:
            return
        _apply_qc_filter(qc_filter_preview(qc.get("min_counts"), qc.get("min_cells")))

    ctx.ensure_qc_filter = _ensure_qc_filter

    def _clear_qc_filter():
        """Revert to every cell, recording nothing.

        There is no code for "un-filter", and none is needed: work recorded
        from here on roots at ``preamble`` again, and the barrier sorts it
        *before* ``qc_filter`` in the notebook, where ``adata`` is still the
        full table. The node is removed only if nothing depends on it -- a
        filtered lineage that exists keeps its step, and ``ProvGraph.remove``
        refusing is how that is decided, not a list of ids that might rot.
        """
        from palms.utils.adata_persistence import full_table
        from palms.utils.rebind_cells import rebind_cells

        if not state.get("qc_filter"):
            return
        state["qc_filter"] = None
        state.pop("_qc_applied_key", None)
        graph = state.get("prov_graph")
        if graph is not None and QC_NODE_ID in graph:
            try:
                graph.remove(QC_NODE_ID)
            except ValueError:
                pass  # something was computed under it; its step stays

        full = full_table(ctx)
        full_l2o = ctx.full_label_to_obs
        rebind_cells(ctx, full, full_l2o if full_l2o is not None else ctx.label_to_obs)
        _save_prov_graph()

    ctx.clear_qc_filter = _clear_qc_filter

    # ── ensure_normalized ────────────────────────────────────────────────
    def _ensure_normalized():
        """Run the ``normalize`` step if needed and return ``adata_norm``.

        Replaces the old ``record_normalize`` + ``get_normalized_adata`` pair,
        which had drifted: the GUI normalised with ``target_sum=1e4`` while the
        recorded cell used scanpy's median default.

        The step binds a *copy* rather than mutating ``adata`` in place. The old
        node was ``kind=SETUP``, so it sorted ahead of every artifact node and
        silently log-normalised the object other cells then copied — a notebook
        containing both it and a self-contained step could double-normalise.
        Consumers now name ``adata_norm`` explicitly.

        Idempotent: re-running is skipped while the source ``adata`` *and* the
        scaling target are unchanged. The target is part of the memo key for the
        same reason ``_ensure_spatial_neighbors`` keys on ``n_neighs`` — without
        it, changing the setting in the Preprocess tab would be a silent no-op
        for the rest of the session, since every later caller would be handed
        the ``adata_norm`` the old setting produced.
        """
        target_sum = state.get("normalize_target_sum", DEFAULT_TARGET_SUM)
        blocks, params, _ = normalize_preview(target_sum)
        key = (id(ctx.adata), target_sum)
        executor = _get_executor()
        if state.get("_norm_src_id") == key and "adata_norm" in executor.ns:
            return executor.ns["adata_norm"]
        _record_preamble()
        outputs = _run_step(Step(
            id=_cell_scoped_id("normalize"),
            **_resolved("normalize", blocks),
            params=params,
            deps=[_cell_root()],
            kind=SETUP,
            label="Normalize, log-transform, PCA",
            outputs=["adata_norm"],
        ))
        state["_norm_src_id"] = key
        return outputs["adata_norm"]

    ctx.ensure_normalized = _ensure_normalized

    # ── record_clustering ────────────────────────────────────────────────
    def _clustering_csv(key):
        """The 10x ``analysis/clustering`` CSV backing *key*, or None.

        Only a clustering that came with the dataset has one. Anything the
        viewer derived (Leiden, CNV, Novae, an import) does not, and the
        producer is responsible for recording the code that made it.
        """
        root = os.path.join(ctx.data_path, "analysis", "clustering")
        for dir_name in (f"gene_expression_{key}", key):
            candidate = os.path.join(root, dir_name, "clusters.csv")
            if os.path.exists(candidate):
                return candidate
        return None

    def _record_clustering(key):
        """Ensure a ``clustering:<key>`` node exists so dependents can name it.

        If the true producer already recorded one (the Leiden tab, the CNV tab,
        Novae, a file import, or a prior session), leave it alone — don't
        overwrite real code with a loader, and don't flag dependents stale.
        """
        graph = state.get("prov_graph")
        if graph is not None and f"clustering:{key}" in graph:
            return
        _record_preamble()
        csv_path = _clustering_csv(key)
        if csv_path is not None:
            code = (
                f"\n# Add clustering: {key}\n"
                f"clust_df = pd.read_csv(r\"{csv_path}\", index_col=0)\n"
                f"adata.obs[\"{key}\"] = pd.Categorical("
                f"clust_df.reindex(adata.obs_names).iloc[:, 0].astype(str).values)"
            )
        else:
            # No CSV and no producer node: the column exists only in the viewer's
            # cache, from a session recorded before its producer recorded code.
            # Reload it, and say so in the cell — the previous version emitted a
            # read_csv of this path regardless, so the exported notebook died with
            # FileNotFoundError on every viewer-derived clustering.
            from palms.utils.clustering_code import reload_clustering_code

            code = reload_clustering_code(
                key,
                os.path.join(ctx.data_path, "sdata_cached.zarr"),
                reason=(f"Clustering '{key}' was computed in an earlier session, "
                        "before its producer recorded code."),
            )
        _record_node(
            f"clustering:{key}",
            code,
            deps=[_cell_root()],  # puts labels into obs; needs no normalisation
            kind=ARTIFACT,
            label=f"Clustering: {key}",
        )

    ctx.record_clustering = _record_clustering

    # ── ensure_spatial_neighbors ─────────────────────────────────────────
    def _ensure_spatial_neighbors(n_neighs):
        """Build the spatial graph on ``adata_norm`` if needed, and record it.

        Replaces ``record_spatial_neighbors``, which recorded a graph built on
        ``adata`` while the viewer built one on the normalised copy the analyses
        were actually handed — so a replayed notebook ran nhood/co-occurrence
        against an object with no ``.obsp`` graph on it.

        Idempotent per ``(adata_norm, n_neighs)``. Re-running with a different
        ``n_neighs`` upserts the node and flags its dependents stale, which is
        the intended semantics.
        """
        n_neighs = coerce(n_neighs)
        _ensure_normalized()
        executor = _get_executor()
        cache_key = (id(executor.ns.get("adata_norm")), n_neighs)
        if state.get("_spatial_neighbors_key") == cache_key:
            return
        _run_step(Step(
            id=_cell_scoped_id("spatial_neighbors"),
            **_resolved("spatial_neighbors", list(builtin_spec("spatial_neighbors").blocks)),
            params={"n_neighs": n_neighs},
            deps=[_cell_scoped_id("normalize")],
            kind=ARTIFACT,
            label="Spatial neighbors",
        ))
        state["_spatial_neighbors_key"] = cache_key

    ctx.ensure_spatial_neighbors = _ensure_spatial_neighbors

    # ── ensure_annotations ───────────────────────────────────────────────
    def _ensure_annotations(preview):
        """Bind and record ``annotations`` from the shapes drawn in the viewer.

        The counterpart of ``tab_roi._record_rois``, and the answer to what was
        recorded here for months as a blocker: the drawn geometry does not have
        to be *reachable* from the notebook, because it is inlined into the
        recorded cell as literals. Both annotation analyses declare this node as
        a dependency, so one drawing produces one node they share.

        Not memoised. Re-running upserts the node with whatever is on the canvas
        now and flags its dependents stale, which is exactly right: an analysis
        run against a different set of shapes is a different analysis.

        Returns ``None`` when nothing typed has been drawn, so a caller can
        report that rather than run an analysis over an empty region.
        """
        blocks, params, _ = preview
        if not params["polygons"]:
            return None
        _record_preamble()
        return _run_step(Step(
            id="annotations",
            **_resolved("annot.polygons", blocks),
            params=params,
            deps=["preamble"],
            kind=SETUP,
            label=f"Annotation polygons ({len(params['polygons'])})",
            outputs=["annotations"],
        ))["annotations"]

    ctx.ensure_annotations = _ensure_annotations

    # ── show_plot ────────────────────────────────────────────────────────
    @ensure_main_thread
    def _add_to_plots_panel(fig, title: str, paths: list):
        """Append to the Plots dock and reveal it. GUI thread only.

        Bounced through ``ensure_main_thread`` because building the thumbnail
        touches ``QPixmap``: figures are routinely *built* in a napari worker,
        and a caller that forgets to hop back would crash Qt rather than show a
        plot.
        """
        panel = ctx.plots_panel
        if panel is None:
            return
        panel.add_figure(fig, title, paths)
        # Through the same route as the View menu: the dock may have been closed
        # (and therefore destroyed) since the last plot, in which case this
        # re-creates it around the panel rather than writing to a dead pointer.
        if ctx.reveal_plots_dock is not None:
            ctx.reveal_plots_dock()

    def _show_plot(fig, stem: str, title: str = None, save: bool = True,
                   paths=None) -> list[str]:
        """Save *fig* under ``<data_path>/plots/`` and show it in the Plots dock.

        The single route every figure takes. Replaces the old
        ``auto_save_plot`` + ``plt.show(block=False)`` pair, which six tabs used,
        four tabs did differently, and four tabs skipped entirely.

        ``save=False`` is for a figure its own Step template already wrote — the
        template must remain the thing that writes the file, or the recorded
        code would stop being the code that ran. Pass the ``paths`` it wrote so
        the card can still say where the figure went.

        Returns the paths written, for the status line and for the recorded
        ``savefig`` argument.

        *fig* may be a scanpy plot object rather than a Figure — the rank-genes
        dotplot is one — so it is resolved here, once, instead of at each of the
        three things that go on to use it.
        """
        figure = to_figure(fig)
        written = [str(p) for p in (paths or [])]
        if save:
            written = save_figure(
                figure, save_paths(ctx.data_path, stem, state=state))
        _add_to_plots_panel(figure, title or stem, written)
        return written

    ctx.show_plot = _show_plot

    def _plot_paths(stem: str) -> list:
        """Where ``stem`` will be written — for a template's ``paths`` param."""
        return [str(p) for p in save_paths(ctx.data_path, stem, state=state)]

    ctx.plot_paths = _plot_paths

    def _recorded_plot_paths(paths) -> list[str]:
        return recorded_paths(ctx.data_path, paths)

    ctx.recorded_plot_paths = _recorded_plot_paths

    # ── refresh_gene_choices ─────────────────────────────────────────────
    def _refresh_gene_choices():
        """Re-populate every gene ComboBox from ``ctx.gene_names``.

        The counterpart of ``refresh_clustering_choices``, and needed for the
        same reason one step later: the gene pickers are built once from
        ``var_names``, so a gene filter (or a segmentation swap onto a table
        with a different panel) leaves them offering genes the table no longer
        has, and the first click is a ``KeyError`` out of ``CellColorManager``.

        A combo whose current value survives keeps it; one whose value is gone
        falls back to the first name rather than to ``None``, which magicgui
        rejects as an invalid choice.
        """
        names = list(ctx.gene_names or [])
        for combo in [ctx.gene_widget, ctx.corr_gene_a_widget,
                      ctx.corr_gene_b_widget, ctx.transcript_gene_widget,
                      ctx.transcript_density_gene_widget, ctx.umap_gene_widget]:
            if combo is None:
                continue
            old_val = combo.value
            combo.choices = names
            if old_val in names:
                combo.value = old_val
            elif names:
                combo.value = names[0]

    ctx.refresh_gene_choices = _refresh_gene_choices

    # ── refresh_clustering_choices ───────────────────────────────────────
    def _refresh_clustering_choices():
        names = list(ctx.clusterings.keys())
        for combo in [ctx.clustering_widget, ctx.ga_clustering_widget,
                      ctx.lr_clustering_widget, ctx.ne_clustering_widget,
                      ctx.co_clustering_widget,
                      ctx.annot_nhood_clustering_widget,
                      ctx.annot_dist_clustering_widget,
                      ctx.mg_clustering_widget,
                      ctx.cnv_clustering_widget]:
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
        # Keep the raw category values as checkbox keys — do NOT coerce to int.
        # CNV clusterings carry *string* categories ('0','1','2', or 'tumor'/'unknown'),
        # and get_cluster_ids_per_obs factorizes them into _cluster_raw_to_id keyed by
        # those raw strings. Coercing the checkbox keys to int broke that lookup, so
        # translate_selected_ids_to_int returned [] and the cluster filter blanked every
        # cell. Sort numerically for display order only; preserve the original type.
        try:
            ids = sorted(raw_ids, key=lambda x: int(x))
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
        _ids = cluster_series.dropna().unique().tolist()
        try:
            ids = sorted(_ids, key=lambda x: int(x))
        except (ValueError, TypeError):
            ids = sorted(_ids, key=lambda x: str(x))
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

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        outer_layout.addWidget(buttons)
        dialog.setLayout(outer_layout)
        dialog.resize(min(800, 250 * n_cols), min(600, 30 * n_per_col + 60))

        if dialog.exec() == QDialog.DialogCode.Accepted:
            new_labels = {cid: e.text() for cid, e in edits.items()}
            if "cluster_labels" not in state or not isinstance(state["cluster_labels"], dict):
                state["cluster_labels"] = {}
            state["cluster_labels"][clustering_key] = new_labels
            from palms.utils.adata_persistence import save_cluster_labels_to_sdata
            save_cluster_labels_to_sdata(ctx, clustering_key, new_labels)
            _record_clustering(clustering_key)
            _lbl_map = {str(k): v for k, v in new_labels.items()}
            _record_node(
                f"annotation:{clustering_key}",
                f"\n# Manual cluster labels for '{clustering_key}'\n"
                f"annotation_map = {_lbl_map!r}\n"
                f"adata.obs[\"{clustering_key}_annotated\"] = ("
                f"adata.obs[\"{clustering_key}\"].astype(str).map(annotation_map)"
                f".astype(\"category\"))",
                deps=[f"clustering:{clustering_key}"],
                label=f"Cluster labels: {clustering_key}",
            )
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

def create_file_menu(viewer, on_open_dataset, on_preprocess_dataset=None):
    """Add Xenium actions to napari's existing File menu (call once at startup)."""
    from qtpy.QtGui import QAction

    file_menu = viewer.window.file_menu

    sep = file_menu.addSeparator()
    sep.setObjectName("palms_sep")

    open_act = QAction("Open Dataset...", file_menu)
    open_act.setObjectName("palms_open")
    open_act.setShortcut("Ctrl+O")
    open_act.triggered.connect(on_open_dataset)
    file_menu.addAction(open_act)

    if on_preprocess_dataset is not None:
        pre_act = QAction("Preprocess Dataset...", file_menu)
        pre_act.setObjectName("palms_preprocess")
        pre_act.triggered.connect(on_preprocess_dataset)
        file_menu.addAction(pre_act)


# ── Plots dock lifecycle ─────────────────────────────────────────────────────
#
# napari's dock title bar has a close button, and it does not hide the dock: it
# calls ``destroyOnClose`` -> ``Window.remove_dock_widget``, which reparents the
# inner widget to ``None`` and ``deleteLater()``s the dock. So one click on the
# "x" left the Plots gallery in a state where the dock object was a dangling C++
# pointer and the panel was an orphan — and the View-menu toggle, which called
# ``setVisible`` on that pointer, raised ``RuntimeError`` into a bare
# ``except: pass`` and appeared to do nothing at all.
#
# The panel itself always survives (Python holds it in ``_app["plots_panel"]``),
# so the fix is to treat the *dock* as disposable and re-create it around the
# surviving panel whenever it is needed.

def dock_is_alive(dock) -> bool:
    """Whether *dock* is still a live Qt object rather than a dangling pointer."""
    if dock is None:
        return False
    try:
        dock.isVisible()
        return True
    except RuntimeError:      # wrapped C/C++ object has been deleted
        return False


def _is_on_a_screen(widget) -> bool:
    """Whether any part of *widget* is on a connected display.

    A floating dock can be dragged somewhere unreachable — off the desktop, or
    onto a screen that has since been disconnected. Re-showing it there looks
    exactly like the toggle doing nothing.
    """
    from qtpy.QtGui import QGuiApplication

    frame = widget.frameGeometry()
    return any(screen.availableGeometry().intersects(frame)
               for screen in QGuiApplication.screens())


def ensure_plots_dock(viewer, _app):
    """The live Plots dock, re-created around the surviving panel if needed.

    Returns ``None`` only when there is no panel yet (before the first dataset
    finishes loading).
    """
    dock = _app.get("plots_dock")
    if dock_is_alive(dock):
        return dock

    panel = _app.get("plots_panel")
    if panel is None:
        return None

    dock = viewer.window.add_dock_widget(panel, name="Plots", area="bottom")
    _app["plots_dock"] = dock

    # Make the title bar's "x" *hide* the dock instead of destroying it. Closing
    # a gallery should not throw away the figures in it, and an instance
    # attribute shadows the class method napari's close button calls.
    dock.destroyOnClose = dock.hide

    if hasattr(dock, "visibilityChanged"):
        def _sync(visible, _dock=dock):
            action = _app.get("plots_action")
            if action is not None and action.isChecked() != visible:
                action.setChecked(visible)
        dock.visibilityChanged.connect(_sync)
    return dock


def reveal_plots_dock(viewer, _app):
    """Show the Plots dock and make sure it is somewhere the user can see it."""
    dock = ensure_plots_dock(viewer, _app)
    if dock is None:
        return None
    # A floating dock that was dragged off the desktop cannot be shown where it
    # is; bring it back into the main window rather than leaving the user with a
    # menu item that appears dead.
    if dock.isFloating() and not _is_on_a_screen(dock):
        dock.setFloating(False)
    dock.setVisible(True)
    dock.raise_()
    dock.activateWindow()
    return dock


# ── View menu ─────────────────────────────────────────────────────────────────

def create_view_menu(viewer, _app):
    """Add View menu items to napari's native View menu (call once at startup)."""
    from qtpy.QtGui import QAction

    view_menu = viewer.window.view_menu
    view_menu.addSeparator()

    # Show/hide the Controls dock panel
    controls_action = QAction("Show Controls", view_menu, checkable=True)
    controls_action.setObjectName("palms_controls_toggle")
    controls_action.setChecked(True)
    controls_action.setShortcut("Ctrl+Shift+X")

    def _on_controls_toggled(checked):
        dw = _app.get("dock_widget")
        if dw is None:
            return
        try:
            if checked:
                from qtpy.QtCore import Qt
                viewer.window._qt_window.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dw)
                dw.show()
            else:
                dw.hide()
        except RuntimeError:
            # C++ object deleted (e.g. during dataset reload) — ignore
            pass

    controls_action.toggled.connect(_on_controls_toggled)
    view_menu.addAction(controls_action)
    _app["controls_action"] = controls_action

    minimap_action = QAction("Show Minimap", view_menu, checkable=True)
    minimap_action.setObjectName("palms_minimap_toggle")
    minimap_action.setChecked(False)
    minimap_action.setEnabled(False)

    def _on_toggled(checked):
        minimap = _app.get("minimap")
        if minimap is not None:
            minimap.setVisible(checked)

    minimap_action.toggled.connect(_on_toggled)
    view_menu.addAction(minimap_action)
    _app["minimap_action"] = minimap_action

    # Show/hide the Plots dock. Starts unchecked — the dock reveals itself when
    # the first figure arrives, so an empty gallery never takes up the canvas.
    plots_action = QAction("Show Plots", view_menu, checkable=True)
    plots_action.setObjectName("palms_plots_toggle")
    plots_action.setChecked(False)
    plots_action.setShortcut("Ctrl+Shift+P")

    def _on_plots_toggled(checked):
        if checked:
            # Not ``dw.setVisible(True)``: napari may have destroyed the dock
            # (its close button does), and calling into the dangling pointer is
            # what made this menu item look broken.
            reveal_plots_dock(viewer, _app)
        else:
            dw = _app.get("plots_dock")
            if dock_is_alive(dw):
                dw.hide()

    plots_action.toggled.connect(_on_plots_toggled)
    view_menu.addAction(plots_action)
    _app["plots_action"] = plots_action


# ── Preferences menu ─────────────────────────────────────────────────────────

def create_preferences_menu(ctx: ViewerContext):
    """Build the Preferences menu on the napari menu bar.

    Reuses an existing Preferences menu if present, clearing stale actions
    from the previous dataset load before repopulating.
    """
    from qtpy.QtWidgets import QActionGroup, QMenu
    from qtpy.QtGui import QAction

    state = ctx.state
    menu_bar = ctx.viewer.window._qt_window.menuBar()

    # Reuse existing Preferences menu if present; clear stale actions from previous dataset
    prefs_menu = None
    for act in menu_bar.actions():
        if act.menu() and act.menu().title() == "Preferences":
            prefs_menu = act.menu()
            break
    if prefs_menu is None:
        prefs_menu = QMenu("Preferences", menu_bar)
        menu_bar.addMenu(prefs_menu)
    else:
        prefs_menu.clear()

    # Plot format. Every plot is written in each of the chosen formats, so the
    # choices are combinations rather than a single extension — PNG to look at,
    # PDF to publish, and the default writes both because asking a user to pick
    # in advance is what left half the figures in a format they had to redo.
    format_menu = prefs_menu.addMenu("Plot format")
    format_group = QActionGroup(format_menu)
    format_group.setExclusive(True)
    current = plot_formats(state)
    for label, formats in (("PNG + PDF", ["png", "pdf"]),
                           ("PNG", ["png"]),
                           ("PDF", ["pdf"]),
                           ("SVG", ["svg"])):
        act = QAction(label, format_group, checkable=True,
                      checked=(formats == current))
        act.setData(formats)
        format_menu.addAction(act)

    def _on_format_changed(action):
        state["plot_formats"] = list(action.data())
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

    # CPU cores — global budget for parallel analyses (currently CopyKAT).
    n_cpu = os.cpu_count() or 2
    current_cores = int(state.get("n_cores", max(1, n_cpu // 2)))
    # A small, machine-scaled set of choices: powers of two up to the core count,
    # plus explicit "half" and "all" (as labelled aliases), deduped and sorted.
    core_choices = sorted({c for c in (1, 2, 4, 8, 16, 32, 64) if c <= n_cpu}
                          | {max(1, n_cpu // 2), n_cpu})
    cores_menu = prefs_menu.addMenu("CPU cores")
    cores_group = QActionGroup(cores_menu)
    cores_group.setExclusive(True)
    for c in core_choices:
        if c == n_cpu:
            label = f"{c} (all)"
        elif c == max(1, n_cpu // 2):
            label = f"{c} (half)"
        else:
            label = str(c)
        act = QAction(label, cores_group, checkable=True, checked=(c == current_cores))
        act.setData(c)
        cores_menu.addAction(act)

    def _on_cores_changed(action):
        state["n_cores"] = int(action.data())
    cores_group.triggered.connect(_on_cores_changed)

    # Record code checkbox
    record_action = QAction("Record reproducible code", prefs_menu, checkable=True, checked=True)
    prefs_menu.addAction(record_action)

    def _on_record_toggled(checked):
        state["record_code"] = checked
        if checked:
            state["code_journal"].clear()
            state["code_journal_tags"].clear()
            state["prov_graph"] = ProvGraph()
            state["_legacy_counter"] = 0
            # Re-seed the preamble so a freshly-restarted recording is valid.
            ctx.record_preamble()
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

    # Continue from existing code file action
    continue_code_action = QAction("Continue from existing code file...", prefs_menu)
    prefs_menu.addAction(continue_code_action)

    def _on_continue_code():
        from pathlib import Path as _Path
        from qtpy.QtWidgets import QFileDialog as _QFD, QMessageBox as _QMB
        path, _ = _QFD.getOpenFileName(
            None, "Continue from Existing Code File",
            str(ctx.data_path), "Python Files (*.py)",
        )
        if not path:
            return
        reply = _QMB.warning(
            None, "Continue from Existing Code File",
            "Continuing from an existing code file may result in incomplete or "
            "duplicated analysis steps. Proceed?",
            _QMB.Yes | _QMB.Cancel,
        )
        if reply != _QMB.Yes:
            return
        with open(path, 'r') as f:
            content = f.read()
        state["code_journal"] = [content.rstrip()]
        state["code_journal_tags"].clear()
        if "from spatialdata_io import xenium" in content:
            state["code_journal_tags"].add("preamble")
        if "sc.pp.normalize_total" in content:
            state["code_journal_tags"].add("normalize")
        for key in ctx.clusterings:
            if f"# Add clustering: {key}" in content:
                state["code_journal_tags"].add(f"clustering_{key}")
        state["code_file"] = _Path(path).name
        ctx.set_status(f"Continuing code recording in {_Path(path).name}")
    continue_code_action.triggered.connect(_on_continue_code)
