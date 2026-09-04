"""Tab 2: Transcript overlay — multi-gene points display, and density binning.

The **points overlay** is display: it reads the per-gene feather index
(`utils/transcript_index.py`, built by `palms-preprocess`) because drawing a
gene has to feel instant, and nothing about which points are on screen belongs
in a notebook.

The **density heatmap** is analysis, and is recorded. It reads
`sdata.points['transcripts']` instead of that index — the index is a viewer
artifact a replayed notebook would not have — which costs a few seconds the
first time a gene is binned and nothing thereafter, since the fetch is its own
step and the bin-size knob only re-runs the histogram.

The **density preview** is display again, and is the reason both statements
above can stay true. Hunting for a gene meant paying the recorded fetch for
every candidate, so the preview bins the feather index instead and redraws in
about a tenth of a second. It draws into its own layer, says so, and records
nothing; `Compute Density` is still the only thing that runs and records the
analysis. It is safe because the two routes were measured to agree exactly —
row for row and bin for bin — and
`tests/test_transcript_density_step.py::test_the_feather_preview_bins_identically_to_the_recorded_step`
keeps them agreeing. It executes the same `transcripts.density` template text
rather than a second histogram, which is the whole point: a hand-rolled copy
would be the unwatched half of a drift pair.
"""

from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np
from magicgui.widgets import ComboBox, CheckBox, PushButton, Slider
from qtpy.QtCore import QTimer
from qtpy.QtWidgets import QListWidget, QHBoxLayout, QWidget, QLabel
from napari.qt.threading import thread_worker
from palms.tabs._helpers import make_tab, StatusProxy
from palms.utils.prov_graph import ARTIFACT, NOTE
from palms.utils.steps import Step, coerce
from palms.utils.step_templates import (
    Preview, builtin_spec, step_template as _resolved,
)

GENE_TEMPLATE_ID = "transcripts.gene"
DENSITY_TEMPLATE_ID = "transcripts.density"

#: Said in the status line, and echoed in the preview layer's name. An
#: unlabelled preview *is* the provenance hazard, so the labelling is asserted
#: by tests/test_preview_never_records.py rather than left to care.
PREVIEW_LABEL = "PREVIEW (not recorded)"
PREVIEW_LAYER_NAME = "transcript_density (PREVIEW - not recorded)"
#: Long enough that dragging the bin-size slider fires once when it stops, not
#: on every tick; short enough to feel like the picture is following you.
PREVIEW_DEBOUNCE_MS = 250

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state

    transcript_gene_widget = ComboBox(
        label="Transcript gene",
        choices=ctx.gene_names,
        value=ctx.gene_names[0] if ctx.gene_names else None,
    )
    # Published on ctx so refresh_gene_choices can re-populate it: the
    # choices are bound to var_names once, here, and a gene filter or a
    # segmentation swap can shrink that panel underneath them.
    ctx.transcript_gene_widget = transcript_gene_widget
    add_gene_button = PushButton(label="Add Gene", enabled=True)
    remove_gene_button = PushButton(label="Remove Selected", enabled=True)
    clear_genes_button = PushButton(label="Clear All", enabled=True)

    gene_list_qt = QListWidget()
    gene_list_qt.setMaximumHeight(150)

    legend_label_qt = QLabel("")
    legend_label_qt.setWordWrap(True)

    transcript_check = CheckBox(label="Show transcripts", value=False)
    qv_slider = Slider(label="Min QV", min=0, max=40, value=20)
    apply_transcripts_button = PushButton(label="Apply Transcripts", enabled=True)

    status_label = StatusProxy(ctx.viewer)

    def _update_legend():
        genes = state["transcript_genes"]
        if not genes:
            legend_label_qt.setText("")
            return
        color_names = [
            "Yellow", "Cyan", "Magenta", "Orange", "Green",
            "Sky Blue", "Red", "Violet", "Pink", "Brown",
        ]
        parts = []
        for i, g in enumerate(genes):
            parts.append(f"{color_names[i % len(color_names)]}: {g}")
        legend_label_qt.setText(" | ".join(parts))

    def on_add_gene():
        gene = transcript_gene_widget.value
        genes = state["transcript_genes"]
        if gene in genes:
            status_label.value = f"'{gene}' already in list"
            return
        if len(genes) >= 10:
            status_label.value = "Maximum 10 genes reached"
            return
        genes.append(gene)
        gene_list_qt.addItem(gene)
        _update_legend()

    def on_remove_gene():
        selected = gene_list_qt.currentRow()
        if selected >= 0:
            gene_list_qt.takeItem(selected)
            state["transcript_genes"].pop(selected)
            _update_legend()

    def on_clear_genes():
        state["transcript_genes"].clear()
        gene_list_qt.clear()
        _update_legend()

    def _on_transcripts_ready(result, _gen):
        if ctx.dataset_generation != _gen:
            return  # dataset reloaded while worker ran
        points, colors = result
        ctx.transcript_layer.data = points
        ctx.transcript_layer.face_color = colors
        ctx.transcript_layer.visible = True
        genes = state["transcript_genes"]
        status_label.value = (
            f"Transcripts: {', '.join(genes)} ({len(points):,} spots total)"
        )
        ctx.record_node(
            "viewer:transcript_overlay",
            f"\n# Display transcript overlay (viewer)\n"
            f"# genes={genes}, qv_threshold={qv_slider.value}, {len(points):,} spots total",
            deps=["preamble"], kind=NOTE, label="Transcript overlay",
        )
        apply_transcripts_button.enabled = True

    def on_apply_transcripts():
        genes = state["transcript_genes"]
        if transcript_check.value and genes:
            status_label.value = f"Loading transcripts for {len(genes)} gene(s)..."
            apply_transcripts_button.enabled = False
            gen = ctx.dataset_generation

            @thread_worker
            def fetch():
                return ctx.transcript_loader.get_multi_gene_points(genes)

            worker = fetch()
            worker.returned.connect(lambda result: _on_transcripts_ready(result, gen))
            worker.start()
        elif transcript_check.value and not genes:
            status_label.value = "No genes in list — add genes first"
        else:
            ctx.transcript_layer.visible = False
            status_label.value = "Transcripts hidden"

    add_gene_button.clicked.connect(on_add_gene)
    remove_gene_button.clicked.connect(on_remove_gene)
    clear_genes_button.clicked.connect(on_clear_genes)
    apply_transcripts_button.clicked.connect(on_apply_transcripts)

    # Button row
    btn_row = QWidget()
    btn_layout = QHBoxLayout()
    btn_layout.setContentsMargins(0, 0, 0, 0)
    btn_layout.addWidget(add_gene_button.native)
    btn_layout.addWidget(remove_gene_button.native)
    btn_layout.addWidget(clear_genes_button.native)
    btn_row.setLayout(btn_layout)

    # ── Transcript Density section ────────────────────────────────────────────
    density_separator = QLabel("── Transcript Density ──")

    density_gene_widget = ComboBox(
        label="Gene",
        choices=ctx.gene_names,
        value=ctx.gene_names[0] if ctx.gene_names else None,
    )
    # Published on ctx so refresh_gene_choices can re-populate it: the
    # choices are bound to var_names once, here, and a gene filter or a
    # segmentation swap can shrink that panel underneath them.
    ctx.transcript_density_gene_widget = density_gene_widget
    bin_size_slider = Slider(label="Bin size (µm)", min=10, max=500, value=50)
    cluster_filter_density_check = CheckBox(label="Filter by selected clusters", value=False)
    normalise_cells_check = CheckBox(label="Normalise by cells per bin", value=False)
    preview_check = CheckBox(
        label="Preview the density as I change settings",
        value=True,
        tooltip=(
            "Ticked: changing the gene, bin size, quality floor or filters redraws "
            "a density picture in about a tenth of a second, binned from the "
            "per-gene transcript index palms-preprocess built. It is a preview. It "
            "is not recorded, it is not in analysis.py, and it is not what the "
            "notebook replays.\n"
            "Unticked: nothing is drawn until you press Compute Density.\n"
            "Either way, Compute Density is the button that runs and records the "
            "analysis, and its result replaces the preview."
        ),
    )
    preview_caption_qt = QLabel(
        "The preview comes from the transcript index and is never recorded — it is "
        "there so you can find the right gene and bin size quickly. Press Compute "
        "Density for the real step: that is the one that lands in analysis.py and "
        "in the exported notebook.")
    preview_caption_qt.setWordWrap(True)
    preview_caption_qt.setStyleSheet("color: gray;")
    compute_density_button = PushButton(label="Compute Density")
    density_status = QLabel("")

    def _cluster_filter_params() -> dict:
        """The two filter params, or nothing at all.

        Read-only, like tab_roi's namesake: recording the clustering node and
        mirroring it onto obs are the run's job, not a side effect of drawing a
        preview pane.
        """
        if not cluster_filter_density_check.value:
            return {}
        key = (ctx.state.get("active_clustering_name")
               or getattr(ctx.clustering_widget, "value", None))
        if not key or key not in (ctx.clusterings or {}):
            return {}
        return {
            "clustering": key,
            "selected": sorted({str(i) for i in ctx.get_selected_cluster_ids()}),
        }

    # ── Preview providers ─────────────────────────────────────────────────────

    def _transcripts_gene_preview() -> Preview:
        return Preview(
            list(builtin_spec(GENE_TEMPLATE_ID).blocks),
            {"gene": density_gene_widget.value,
             # The QV slider is provenance-only for the point overlay, whose
             # feather files were filtered when palms-preprocess built them.
             # Here it is a real filter, applied as the transcripts are read.
             "min_qv": int(coerce(qv_slider.value))},
        )

    def _density_blocks(filtered: bool, normalised: bool) -> list[str]:
        """Both switches are blocks as well as params: "filter by cluster" adds
        a mask, "normalise" adds a second histogram."""
        blocks = ["head"]
        if filtered:
            blocks.append("filter")
        blocks.append("main")
        if normalised:
            blocks.append("normalise")
        return blocks

    def _density_preview() -> Preview:
        params = {
            "gene": density_gene_widget.value,
            "bin_size_um": float(coerce(bin_size_slider.value)),
            "pixel_size": float(coerce(ctx.pixel_size)),
        }
        filter_params = _cluster_filter_params()
        params.update(filter_params)
        return Preview(
            _density_blocks(bool(filter_params), normalise_cells_check.value),
            params,
        )

    ctx.state.setdefault("template_preview", {})[GENE_TEMPLATE_ID] = _transcripts_gene_preview
    ctx.state.setdefault("template_preview", {})[DENSITY_TEMPLATE_ID] = _density_preview

    # ── Density preview: display only, never recorded ────────────────────────
    # Everything below draws a picture and stops. It must not call run_step,
    # record_*, or ensure_* — a preview with a side effect is a preview that
    # changed the analysis, and tests/test_preview_never_records.py fails if one
    # appears here. The Step it builds carries a "preview:" id so that even a
    # bug that handed it to run_step could not overwrite a recorded node.

    def _paint_bins(layer, density, bin_size_um):
        """Shared layer plumbing for both the preview and the recorded result.

        The scale is the bin's size **in microns**, because that is the world
        every other layer lives in: `app.py` stamps each layer on insertion with
        `units.apply_to_layer`, which sets `layer.scale = pixel_size` — napari
        0.8 removed `ScaleBarOverlay.unit`, so the magnitude has to live in the
        scale (see utils/units.py). This used to assign the bin size in *image
        pixels* (`bin_size_um / pixel_size`), which overwrote that micron scale
        with a number 1/pixel_size too large: at 10 um bins on a 0.2125 um/px
        dataset the grid drew 4.7x oversized, spanning 34,120 um against a
        7,258 um image. Both layers had it; the preview only made it visible,
        by being the one on screen.
        """
        layer.data = density.astype(np.float32)
        layer.scale = [bin_size_um, bin_size_um]
        nonzero = density[density > 0]
        layer.contrast_limits = (
            [0, float(np.percentile(nonzero, 99))] if nonzero.size else [0, 1])
        layer.visible = True
        return nonzero

    def _hide_preview(reason: str = ""):
        """Take the preview off screen. A refusal must never leave a stale one.

        Only the *preview's own* status line is replaced: a recorded run's
        "Done — …" is the answer to a question the user actually asked, and
        turning the preview off is not a reason to take it away.
        """
        if ctx.transcript_bins_preview_layer is not None:
            ctx.transcript_bins_preview_layer.visible = False
        if reason or density_status.text().startswith(PREVIEW_LABEL):
            density_status.setStyleSheet("")
            density_status.setText(reason)

    def _preview_bindings():
        """The rendered step and its feather-derived points, or (None, reason).

        ``need_cell_id`` is read off the *rendered source about to be executed*
        rather than from which blocks are selected. Only the cluster filter
        reads ``cell_id`` today, but a user override that starts using it
        elsewhere is caught too, and the failure direction is a refused preview
        rather than a wrong one.
        """
        from palms.utils.transcript_index import points_for_preview

        loader = getattr(ctx, "transcript_loader", None)
        if loader is None or ctx.sdata is None:
            return None, None, ""
        gene = density_gene_widget.value
        if gene is None:
            return None, None, ""
        if "transcripts" not in getattr(ctx.sdata, "points", {}):
            return None, None, ""
        if "morphology_focus" not in getattr(ctx.sdata, "images", {}):
            return None, None, ""

        blocks, params, _ = _density_preview()
        if cluster_filter_density_check.value and "clustering" not in params:
            return None, None, "no preview — no clustering applied to filter by"
        key = params.get("clustering")
        if key and key not in getattr(ctx.adata, "obs", {}):
            # Mirroring the clustering onto obs is the *run's* job; doing it
            # here would make drawing a picture mutate the analysis.
            return None, None, (
                "no preview — press Compute Density once to apply the clustering")

        step = Step(
            id=f"preview:transcript_density:{gene}",
            **_resolved(DENSITY_TEMPLATE_ID, blocks),
            params=params, kind=ARTIFACT, outputs=["transcript_density"],
        )
        points, reason = points_for_preview(
            loader, ctx.sdata, gene, int(coerce(qv_slider.value)),
            need_cell_id="cell_id" in step.render())
        if points is None:
            return None, None, reason
        return step, points, ""

    def _on_preview_ready(density, token, generation, bin_size_um, normalised):
        if token != state.get("_preview_token") or generation != ctx.dataset_generation:
            return                      # superseded while it was in flight
        nonzero = _paint_bins(ctx.transcript_bins_preview_layer, density, bin_size_um)
        ctx.transcript_bins_layer.visible = False
        what = (f"{int(nonzero.size):,} non-zero bins, normalised by cells"
                if normalised else f"{int(density.sum()):,} transcripts binned")
        density_status.setStyleSheet("color: #b8860b;")
        density_status.setText(
            f"{PREVIEW_LABEL} — {what}. Press Compute Density to run and record it.")

    def _preview_start():
        if not preview_check.value or state.get("_density_running"):
            return
        step, points_or_none, reason = _preview_bindings()
        if step is None:
            _hide_preview(reason)
            return

        state["_preview_token"] = token = state.get("_preview_token", 0) + 1
        generation = ctx.dataset_generation
        bin_size_um = step.params["bin_size_um"]
        normalised = normalise_cells_check.value

        @thread_worker
        def _work():
            return ctx.preview_step(
                step, bindings={"transcript_points": points_or_none},
            )["transcript_density"]

        worker = _work()
        worker.returned.connect(
            lambda d: _on_preview_ready(d, token, generation, bin_size_um, normalised))
        worker.errored.connect(lambda e: _hide_preview(f"no preview ({e})"))
        worker.start()

    _preview_timer = QTimer()
    _preview_timer.setSingleShot(True)
    _preview_timer.setInterval(PREVIEW_DEBOUNCE_MS)
    _preview_timer.timeout.connect(_preview_start)

    # Note the prefix naming above: tests/test_tab_templates.py treats every
    # function whose name *ends* in ``_preview`` as a template-preview provider
    # and requires a direct call to it, which a signal connection is not. Hence
    # ``_preview_start`` rather than ``_start_preview``.
    def _schedule_preview_refresh(*_args):
        if not preview_check.value:
            _hide_preview()
            return
        _preview_timer.start()          # restart, so a drag fires once at the end

    def _on_density_ready(result):
        density, bin_size_um, normalised, needs_clustering = result
        compute_density_button.enabled = True
        state["_density_running"] = False
        # The recorded result supersedes the preview, and says so by taking the
        # screen back: the preview layer goes away rather than sitting under it.
        if ctx.transcript_bins_preview_layer is not None:
            ctx.transcript_bins_preview_layer.visible = False
        density_status.setStyleSheet("")
        if density is None:
            return
        nonzero = _paint_bins(ctx.transcript_bins_layer, density, bin_size_um)
        if normalised:
            density_status.setText(
                f"Done — {int(nonzero.size):,} non-zero bins, normalised by cells")
        else:
            density_status.setText(f"Done — {int(density.sum()):,} transcripts binned")

    def on_compute_density():
        gene = density_gene_widget.value
        if gene is None:
            return
        if "transcripts" not in getattr(ctx.sdata, "points", {}):
            density_status.setText("No transcripts element — cannot compute density")
            return
        if "morphology_focus" not in getattr(ctx.sdata, "images", {}):
            density_status.setText("No morphology image — cannot compute density")
            return

        gene_blocks, gene_params, _ = _transcripts_gene_preview()
        blocks, params, _ = _density_preview()
        if cluster_filter_density_check.value and "clustering" not in params:
            density_status.setText("No clustering applied — filter skipped")
            return

        # The clustering must exist as a node, and in adata.obs, before the step
        # can declare it as a dependency and read it. Both on the GUI thread.
        deps = [f"transcripts:{gene}"]
        clustering_key = params.get("clustering")
        if clustering_key:
            from palms.utils.gene_analysis import add_clustering_to_obs
            ctx.record_clustering(clustering_key)
            add_clustering_to_obs(ctx.adata, ctx.adata,
                                  ctx.clusterings[clustering_key], clustering_key)
            deps.append(f"clustering:{clustering_key}")

        gene_step = Step(
            id=f"transcripts:{gene}",
            **_resolved(GENE_TEMPLATE_ID, gene_blocks),
            params=gene_params,
            deps=[ctx.cell_root()],
            kind=ARTIFACT,
            label=f"Transcripts: {gene}",
            outputs=["transcript_points"],
        )
        density_step = Step(
            id=f"transcript_density:{gene}",
            **_resolved(DENSITY_TEMPLATE_ID, blocks),
            params=params,
            deps=deps,
            kind=ARTIFACT,
            label=f"Transcript density: {gene}",
            outputs=["transcript_density"],
        )

        # Fetching a gene's transcripts is the slow half — seconds, against
        # milliseconds for the histogram — so it is its own step and runs only
        # when the gene changes. Sliding the bin size re-runs the cheap half.
        fetch_key = (gene, gene_params["min_qv"], ctx.dataset_generation)
        need_fetch = state.get("_transcript_points_key") != fetch_key

        compute_density_button.enabled = False
        # Stops a queued preview from starting behind the real thing: they share
        # the executor lock, and a preview is never worth making a user wait.
        state["_density_running"] = True
        _preview_timer.stop()
        density_status.setStyleSheet("")
        density_status.setText(
            f"Fetching transcripts for {gene}..." if need_fetch
            else f"Computing density for {gene}...")
        normalised = normalise_cells_check.value
        bin_size_um = params["bin_size_um"]

        @thread_worker
        def _run():
            ctx.record_preamble()
            if need_fetch:
                ctx.run_step(gene_step)
                state["_transcript_points_key"] = fetch_key
            density = ctx.run_step(density_step)["transcript_density"]
            return density, bin_size_um, normalised, clustering_key

        def _failed(exc):
            compute_density_button.enabled = True
            state["_density_running"] = False
            state.pop("_transcript_points_key", None)
            density_status.setText(f"Density failed: {exc}")

        worker = _run()
        worker.returned.connect(_on_density_ready)
        worker.errored.connect(_failed)
        worker.start()

    compute_density_button.clicked.connect(on_compute_density)

    # The density controls become live: until now nothing here redrew without
    # the button. Min QV is included because it changes which rows are read.
    for _live in (density_gene_widget, bin_size_slider, qv_slider,
                  cluster_filter_density_check, normalise_cells_check,
                  preview_check):
        _live.changed.connect(_schedule_preview_refresh)

    widget = make_tab(
        transcript_gene_widget,
        transcript_check,
        qv_slider,
        apply_transcripts_button,
        gene_list_qt,
        btn_row,
        legend_label_qt,
        density_separator,
        density_gene_widget,
        bin_size_slider,
        cluster_filter_density_check,
        normalise_cells_check,
        preview_check,
        preview_caption_qt,
        compute_density_button,
        density_status,
    )
    return widget, {}
