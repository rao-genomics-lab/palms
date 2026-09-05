"""Tools -> Preprocess: how the expression matrix is normalised.

``normalize`` is the one step with no owning tab. It is pulled in implicitly by
``ctx.ensure_normalized()`` from nine call sites --- Clustering, Rank Genes,
Markers, UMAP, Correlation, Co-occurrence, the annotation neighbourhood tab and
(via ``ensure_spatial_neighbors``) the spatial statistics --- so there was
nowhere to put a control for it, and its scaling target stayed a literal in the
template text. This tab is that home.

Its reach is the reason it sits in Tools beside QC rather than inside
Clustering: QC answers "which cells?", this answers "on what scale?", and both
change what every later analysis is about. The tab order reads filter ->
normalise because that is the order the steps run in, and the order the barrier
puts them in when a filter is in force.

Changing the setting does not recompute anything on its own. The next analysis
that needs ``adata_norm`` re-runs the step, which revises the ``normalize`` node
and flags whatever depended on it stale --- which is honest: those results were
computed on a differently scaled matrix.
"""
from __future__ import annotations

from magicgui.widgets import CheckBox, FloatSpinBox
from qtpy.QtWidgets import QLabel

from palms.tabs._helpers import (
    DEFAULT_TARGET_SUM, StatusProxy, make_tab, normalize_preview,
)
from palms.utils.step_templates import Preview
from palms.utils.viewer_context import ViewerContext

TEMPLATE_ID = "normalize"


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    status = StatusProxy(ctx.viewer)

    # ── Widgets ──────────────────────────────────────────────────────────
    # Labels read as English; the tooltip names the template parameter, the
    # convention the Clustering tab established.
    median_check = CheckBox(
        label="Median counts (scanpy default)", value=False,
        tooltip="Scale every cell to the median count across cells, by passing\n"
                "no target_sum at all — scanpy's own default, and what most\n"
                "published scanpy/squidpy pipelines use.\n\n"
                "Untick to scale to a fixed number instead.",
    )
    target_spin = FloatSpinBox(
        label="Counts per cell", min=100.0, max=1e6, step=1000.0,
        value=DEFAULT_TARGET_SUM,
        tooltip="Every cell is scaled to this many counts before log1p.\n"
                "10,000 is the long-standing convention and the viewer's\n"
                "historical behaviour.\n\n"
                "Template parameter: target_sum",
    )

    readout = QLabel("")
    readout.setWordWrap(True)
    readout.setTextFormat(_rich())

    scope = QLabel(
        "Applies to <code>adata_norm</code>, which Clustering, Rank Genes, "
        "Markers, UMAP, Correlation and the spatial statistics all read. "
        "<b>ROI DEG and inferCNV normalise their own copies and are not yet "
        "affected</b> — they still scale to 10,000 whatever is set here."
    )
    scope.setWordWrap(True)
    scope.setTextFormat(_rich())

    # ── Provider ─────────────────────────────────────────────────────────
    def _normalize_preview() -> Preview:
        """What the next expression-based analysis would normalise with.

        Block selection is delegated to ``_helpers.normalize_preview`` so this
        and ``ctx.ensure_normalized`` turn the same setting into the same
        blocks. There is no "Run" button here to pair it with — the step runs
        when something needs it — so the readout below is what calls this, and
        the pane in Tools -> Templates renders from the same call.
        """
        return normalize_preview(_target_sum())

    ctx.state.setdefault("template_preview", {})[TEMPLATE_ID] = _normalize_preview

    # ── Readout ──────────────────────────────────────────────────────────
    def _target_sum():
        """The setting as the template wants it: a float, or None for median."""
        return None if median_check.value else float(target_spin.value)

    def _refresh_readout():
        target_spin.enabled = not median_check.value
        state["normalize_target_sum"] = _target_sum()
        blocks, params, _ = _normalize_preview()
        call = ("sc.pp.normalize_total(adata_norm)" if "scale.median" in blocks
                else f"sc.pp.normalize_total(adata_norm, "
                     f"target_sum={params['target_sum']!r})")
        readout.setText(f"Records: <code>{call}</code>")

    for widget in (median_check, target_spin):
        widget.changed.connect(lambda *_: _on_change())

    def _on_change():
        """Report what a change costs, since nothing recomputes here.

        Silent would be wrong in both directions: a user who expects an
        immediate re-clustering gets none, and a user who does not expect their
        existing results to go stale gets that anyway on the next run.
        """
        _refresh_readout()
        if state.get("_norm_src_id") is not None:
            status.value = (
                "Normalisation changed — the next analysis re-runs it, and "
                "results computed on the old scaling will show as stale."
            )

    # ── restore_session ──────────────────────────────────────────────────
    def _restore_session(session):
        """Display only; ``app.py`` has already seeded the state key.

        Same division as the QC tab: this handler runs after every other tab's,
        and one of them reaching ``ensure_normalized`` first would normalise on
        whatever the default was.
        """
        target = session.get("normalize_target_sum",
                             state.get("normalize_target_sum",
                                       DEFAULT_TARGET_SUM))
        median_check.value = target is None
        if target is not None:
            target_spin.value = float(target)
        _refresh_readout()

    _refresh_readout()

    widget = make_tab(median_check, target_spin, readout, scope)
    return widget, {"restore_session": _restore_session}


def _rich():
    from qtpy.QtCore import Qt

    return Qt.TextFormat.RichText
