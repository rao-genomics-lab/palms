"""Tab: Publish — export the dataset as Celldega DegaFiles.

A PALMS session already ends as a replayable notebook. This is the other half:
a *viewable* artifact. DegaFiles are WebP Deep Zoom pyramids plus Parquet vector
data, which Celldega's ``Landscape`` widget renders in a browser with no server
and no install, so a collaborator without the raw 10x output — or without a
Linux box to run a Xenium viewer on — can still look at the section.

The conversion is celldega's; ``utils/dega_export`` is the wrapper that keeps it
from writing into the raw output, and that module's docstring is where the
reasoning lives. This tab is the button, and three things about it are decided
by facts measured there rather than by taste:

* **The export is 322 s** on a 10.6 x 6.3 mm section, so it runs in a
  ``thread_worker`` with an elapsed readout. A blocking button would look like
  a hung viewer for five minutes.
* **celldega is an optional dependency and is checked lazily.** Importing it at
  build time would cost every launch a spatialdata-sized import for a tab most
  users never open, so nothing here touches it until the user asks — and when
  they do, ``dega_available`` is the single answer, used for the "Check" button
  and for the pre-flight inside the worker alike.
* **The destination is fixed at ``<data_path>/degafiles``.** Not a chooser: the
  recorded cell would then carry an absolute path, and a notebook that replays
  on another machine would either fail or write somewhere surprising. It is a
  deliverable to copy, not a mount point to write through.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from magicgui.widgets import ComboBox, PushButton, SpinBox
from napari.qt.threading import thread_worker
from qtpy.QtCore import QTimer
from qtpy.QtWidgets import QLabel, QTextEdit

from palms.tabs._helpers import (
    StatusProxy, attach_tqdm_progress, make_progress_bar, make_tab,
    qt_tqdm_context,
)
from palms.utils.dega_export import (
    IMAGE_TILE_LAYER, TILE_SIZE_UM, clear_staging, dega_available,
    NotExportable, degafiles_dir, is_cache_only, require_exportable,
    staging_dir,
)
from palms.utils.plot_output import recorded_paths
from palms.utils.prov_graph import TERMINAL
from palms.utils.steps import Step, StepError, coerce
from palms.utils.step_templates import (
    Preview, builtin_spec, step_template as _resolved,
)

TEMPLATE_ID = "export.degafiles"

#: The node id. One per dataset rather than one per run: a second export
#: supersedes the first in the same directory, so re-running revises the node
#: in place, which is what ``upsert`` already does for every other artifact.
NODE_ID = "export:degafiles"

if TYPE_CHECKING:
    from palms.utils.viewer_context import ViewerContext


def build_tab(ctx: ViewerContext) -> tuple:
    state = ctx.state
    status = StatusProxy(ctx.viewer) if ctx.viewer is not None else _NullStatus()

    intro = QLabel(
        "Write a Celldega DegaFile set — a self-contained, browser-viewable "
        "copy of this dataset — to <code>degafiles/</code> beside the data. "
        "Needs the optional <code>dega</code> extra."
    )
    intro.setWordWrap(True)
    intro.setTextFormat(_rich())

    layer_widget = ComboBox(
        label="Image layer", choices=["all", "dapi"], value=IMAGE_TILE_LAYER,
        tooltip="Which morphology channels to tile.\n\n"
                "all — DAPI, boundary, RNA and protein: the four\n"
                "morphology_focus images a Xenium bundle ships, and the four\n"
                "the Landscape viewer offers as toggles. Four pyramids, so\n"
                "roughly four times the runtime and the output size.\n\n"
                "dapi — the nuclear channel alone. Much faster, and enough\n"
                "when the image is only a backdrop for the cells.\n\n"
                "Template parameter: image_tile_layer",
    )
    tile_widget = SpinBox(
        label="Tile size (µm)", min=50, max=2000, step=50, value=TILE_SIZE_UM,
    )
    # celldega parallelises its tile generation over processes. One by default
    # because that is the configuration the 322 s measurement was taken in, and
    # a default that has never been run is not a default.
    workers_widget = SpinBox(label="Worker processes", min=1, max=16, value=1)

    destination = QLabel(f"<code>{degafiles_dir(ctx.data_path)}</code>")
    destination.setWordWrap(True)
    destination.setTextFormat(_rich())

    check_btn = PushButton(label="Check Celldega")
    publish_btn = PushButton(label="Publish DegaFiles")
    # A cache-only dataset is stated on arrival rather than on failure. Not
    # disabled, though: the button is what prints the full explanation, and a
    # dead control with no reason beside it is the worse of the two.
    if is_cache_only(ctx.data_path):
        intro.setText(
            intro.text() + "<br><br><b>This dataset has no raw 10x output.</b> "
            "It is a Crop Dataset export, whose SpatialData zarr <i>is</i> the "
            "data; celldega reads the original 10x bundle, so there is nothing "
            "here for it to convert. Publish the dataset this one came from."
        )
    clear_btn = PushButton(label="Clear Staging Files")

    progress = make_progress_bar()

    report = QTextEdit()
    report.setReadOnly(True)
    report.setFontFamily("monospace")
    report.setMinimumHeight(140)

    # ── Preview provider ─────────────────────────────────────────────────────
    def _publish_preview() -> Preview:
        """What "Publish DegaFiles" would run, with the widgets as they stand.

        ``out_dir`` is recorded *relative to the dataset* — the convention
        ``plot_output.recorded_paths`` sets for every other written artifact —
        so the replayed cell writes beside whatever data it was pointed at
        rather than at this machine's copy of it.
        """
        out = degafiles_dir(ctx.data_path)
        return Preview(
            list(builtin_spec(TEMPLATE_ID).blocks),
            {
                "out_dir": recorded_paths(ctx.data_path, [out])[0],
                "tile_size": coerce(tile_widget.value),
                "image_tile_layer": layer_widget.value,
                "max_workers": coerce(workers_widget.value),
            },
        )

    state.setdefault("template_preview", {})[TEMPLATE_ID] = _publish_preview

    # ── Availability ─────────────────────────────────────────────────────────
    def _exportable() -> str:
        """Empty when this dataset can be published, else why it cannot.

        Filesystem-only, so it costs nothing to call while building the tab —
        which is the point: a Crop Dataset export can *never* be published, and
        finding that out from celldega takes minutes and arrives as
        ``CalledProcessError: ['gzip', '-dk', 'cells.csv.gz']``.
        """
        try:
            require_exportable(ctx.data_path)
        except NotExportable as exc:
            return str(exc)
        return ""

    def _on_check():
        blocker = _exportable()
        ok, message = dega_available()
        parts = []
        if ok:
            import celldega
            parts.append(
                f"Celldega {getattr(celldega, '__version__', '?')} is importable, "
                f"with pyvips.")
        else:
            parts.append(message)
        if blocker:
            parts.append(blocker)
            status.value = "This dataset has no raw 10x output to publish."
        elif ok:
            parts.append("Ready to publish.")
            status.value = "Celldega is available."
        else:
            status.value = "Celldega is not installed."
        report.setPlainText("\n\n".join(parts))

    # ── Publish ──────────────────────────────────────────────────────────────
    def _on_publish():
        # Cheapest refusal first, and the one that can never be satisfied: a
        # dataset with no raw 10x output has nothing for celldega to read, so
        # asking whether celldega is installed would be beside the point.
        blocker = _exportable()
        if blocker:
            report.setPlainText(blocker)
            status.value = "This dataset has no raw 10x output to publish."
            return

        # The pre-flight is here and not in the worker so the remedy text lands
        # instantly; the worker's own require_dega() still guards the call, for
        # the notebook path where nobody pressed a button.
        ok, message = dega_available()
        if not ok:
            report.setPlainText(message)
            status.value = "Celldega is not installed."
            return

        blocks, params, _ = _publish_preview()
        step = Step(
            id=NODE_ID,
            **_resolved(TEMPLATE_ID, blocks),
            params=params,
            deps=["preamble"],
            kind=TERMINAL,
            label="Publish DegaFiles",
            outputs=["degafiles_path"],
        )

        publish_btn.enabled = False
        report.setPlainText(
            "Exporting… this takes several minutes — 322 s for a 10.6 x 6.3 mm "
            "section, most of it tiling the morphology image.\n"
            "The viewer stays usable; the raw data is not written to."
        )

        elapsed = QTimer()
        seconds = [0]

        def _tick():
            seconds[0] += 1
            status.value = f"Publishing DegaFiles… {seconds[0]} s elapsed"

        elapsed.timeout.connect(_tick)
        elapsed.start(1000)

        _post = [None]

        @thread_worker
        def _run():
            try:
                # celldega reports its tiling through tqdm; relaying it costs
                # nothing when a future version stops doing so.
                with qt_tqdm_context(_post[0], "DegaFile export: "):
                    return True, str(ctx.run_step(step)["degafiles_path"])
            except StepError as exc:
                return False, str(exc)

        def _on_done(result):
            elapsed.stop()
            publish_btn.enabled = True
            succeeded, detail = result
            if not succeeded:
                report.setPlainText(f"Export failed:\n{detail}")
                status.value = "DegaFile export failed."
                return
            report.setPlainText(
                f"Wrote DegaFiles to\n  {detail}\n"
                f"in {seconds[0]} s.\n\n"
                f"Open it with celldega's Landscape widget:\n"
                f"    from celldega.viz import landscape\n"
                f"    landscape(base_url='{detail}')\n\n"
                f"Staging files (safe to clear) are under\n"
                f"  {staging_dir(ctx.data_path).parent}"
            )
            status.value = f"Published DegaFiles to {detail}"

        def _on_error(exc):
            elapsed.stop()
            publish_btn.enabled = True
            report.setPlainText(f"Export failed:\n{exc}")
            status.value = "DegaFile export failed."

        worker = _run()
        _post[0], state["_dega_progress_timer"] = attach_tqdm_progress(
            worker, lambda m: setattr(status, "value", m),
            "DegaFile export: ", progress_bar=progress,
        )
        worker.returned.connect(_on_done)
        worker.errored.connect(_on_error)
        # Held so the elapsed timer is not collected while the worker runs.
        state["_dega_elapsed_timer"] = elapsed
        worker.start()

    # ── Staging ──────────────────────────────────────────────────────────────
    def _on_clear():
        """Drop the symlink farm and everything celldega extracted into it.

        Safe at any time and never touches the export: the farm is a build
        input, under ``viewer_cache/``, which the Dataset tab already treats as
        the viewer's own. Clearing it only makes the next export slower, since
        celldega skips an archive it finds already unpacked.
        """
        staging = staging_dir(ctx.data_path).parent
        if not staging.exists():
            report.setPlainText(f"No staging files to clear ({staging}).")
            return
        clear_staging(ctx.data_path)
        report.setPlainText(f"Cleared {staging}.")
        status.value = "Cleared DegaFile staging files."

    check_btn.changed.connect(_on_check)
    publish_btn.changed.connect(_on_publish)
    clear_btn.changed.connect(_on_clear)

    def _restore_session(session):
        pass  # Nothing here is session state; the export is a file on disk.

    scroll = make_tab(
        intro,
        layer_widget,
        tile_widget,
        workers_widget,
        QLabel("Destination:"),
        destination,
        check_btn,
        publish_btn,
        progress,
        report,
        clear_btn,
    )
    return scroll, {"restore_session": _restore_session}


def _rich():
    from qtpy.QtCore import Qt
    return Qt.TextFormat.RichText


class _NullStatus:
    """Stand-in for the napari status bar when there is no viewer (tests)."""

    value = ""
