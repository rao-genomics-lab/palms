# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Xenium Viewer** is a napari-based spatial transcriptomics viewer for Xenium 3.x output data — an open alternative to the commercial Xenium Explorer, which has no Linux build. It visualizes high-resolution spatial gene expression data with cell-level resolution. Supported on **Linux, macOS and WSL2**; Linux is the primary development platform and native Windows is not supported.

## Running the Viewer

```bash
# First-time setup. install.sh is `conda env create -f environment.yml` plus, on
# Linux/WSL only, `conda env update -f environment-linux.yml`. That overlay carries
# libglx-devel, which is linux-only and made environment.yml unsolvable on macOS
# while it lived there — conda env files have no platform selectors, hence a script.
./scripts/install.sh
conda activate xenium_viewer

# Optional one-time transcript preprocessing per dataset (~30-60 min)
# Without this, transcript loading falls back to scanning transcripts.parquet (~5s/gene instead of <100ms)
xenium-preprocess /path/to/xenium/output/

# Optional: build the SpatialData zarr cache without the GUI (the viewer builds
# it on first launch anyway). --check reports on an existing one and writes nothing.
xenium-build-cache /path/to/xenium/output/

# Rename or move a dataset. Nothing stores a dataset *name*, but the provenance
# graph records absolute paths; this rewrites them. --dry-run writes nothing.
xenium-rename-dataset /path/to/xenium/output/ new_name
xenium-rename-dataset /path/that/was/moved/ --repair

# Launch viewer (file dialog opens if no path given)
xenium-viewer [/path/to/xenium/output/]

# Launch without SpatialData zarr cache
xenium-viewer /path/to/xenium/output/ --no-cache
```

The package is installed as `xenium-viewer` (PyPI name) / `xenium_viewer` (import name) via `pip install -e .` (handled automatically by `environment.yml`). Console scripts: `xenium-viewer`, `xenium-preprocess`, `xenium-build-cache`, `xenium-rename-dataset`, `xenium-fetch-references`, `xenium-build-custom-segmentation`. You can also run `python -m xenium_viewer ...`.

There is a `pytest` suite in `tests/` (~320 tests) covering pure logic (provenance graph,
step templates, CopyKAT subsampling, registration math, LLM parsing, notebook export) and
the **zarr/SpatialData persistence paths** — crash-safe writes with simulated interrupted
writes, cache verify/repair, session save, loader cache policy and sidecar locations.
Several are *source guards* that fail if a fixed bug is reintroduced (calling
`delete_element_from_disk` outside `zarr_safe.py`, `rmtree`ing the live cache, writing a
sidecar into the store root, printing a warning instead of logging it).

`tests/test_notebook_replay.py` is the end-to-end one: it runs the real `Step` templates
over a synthetic AnnData, exports the provenance graph as a real `.ipynb`, executes it in
a **clean kernel** (`nbclient`, `allow_errors=False`), and requires ARI **exactly 1.0** on
every clustering plus identical top-N ranked genes. It takes ~35 s — most of it importing
scanpy/squidpy in the kernel — which is why the fixture is module-scoped. Its one
substitution is the `preamble` node (h5ad instead of `xenium(data_path)`, since CI has no
dataset); a test asserts that stays the only one. See "Verifying the claim" below.

Run with plain `pytest` from the repo root (`[tool.pytest.ini_options]` sets
`pythonpath = ["src"]`, so no install is needed). **No environment variables are
needed** — `tests/conftest.py` selects `QT_QPA_PLATFORM=offscreen` and
`MPLBACKEND=Agg` itself when there is no `DISPLAY`, because Qt does not *fail*
without a platform plugin, it aborts the process: the symptom is
`Aborted (core dumped)` with no test name, which is what CI reported for months.
An explicit `QT_QPA_PLATFORM` or a real display still wins, so a desktop run is
unchanged. On a desktop the platform stays unset, which is why
`reporting._headless()` does not rely on it alone: it also checks
`QThread.loopLevel()`, since a `pytest` process has a `QApplication` but no event
loop, and a modal `exec_()` there blocks forever with nothing to dismiss it.
GitHub Actions
(`.github/workflows/ci.yml`) runs the suite in the full conda env plus a fast `ruff`
error-only lint gate (`--select E9,F63,F7,F82`) on every push/PR. The napari GUI proper
and the spatial-analysis tabs still have no automated coverage — that testing remains
manual/exploratory.

## Architecture

### Entry Points & Load Sequence

1. **`src/xenium_viewer/app.py`** — Main entry point (~1300 lines). Validates data dir, orchestrates loading, builds napari viewer with layers, instantiates all managers, creates `ViewerContext`, builds all tab widgets, then restores session. The `main()` function is the `xenium-viewer` console-script entry point.
2. **`src/xenium_viewer/loader.py`** — Loads SpatialData from Xenium output. Uses a zarr cache (`sdata_cached.zarr/`) for 60–70% faster subsequent launches; staleness comes from a content hash in `.xv_manifest.json` (see "Cache safety"). Public API: `load_sdata`, `load_umap`, `load_clusterings`, `get_label_to_obs_mapping`. Its `main()` is the `xenium-build-cache` console script — the same load, headless, so the cold read can happen over ssh. `load_sdata(on_stale=)` is how a caller with no GUI authorises a rebuild: `_stale_preference` checks it *before* any branch that would prompt, because two of those branches never reach a dialog. Default `None` keeps the GUI behaviour exactly as it was, and `'keep'` deliberately does not satisfy `_ask_corrupt_cache` — a cache that will not open cannot be kept.
3. **`src/xenium_viewer/preprocess.py`** — One-time step that splits the transcript parquet into ~480 per-gene feather files for fast per-gene loading (~100ms vs 4–5s). The `main()` function is the `xenium-preprocess` console-script entry point.

### Central State: ViewerContext

**`src/xenium_viewer/utils/viewer_context.py`** — `ViewerContext` dataclass is the single shared state object passed to all tab modules. It holds:
- Core data: `viewer`, `adata`, `sdata`, `clusterings`, pixel size
- Layer references: `cell_labels_layer`, `transcript_layer`, `roi_layer`
- Manager objects: `CellColorManager`, `TranscriptLoader`, `UMAPViewer`
- Mutable state dicts: `state` (general), `he_state` (H&E registration), `arms_state` (ARMS overlay)
- Shared callables: `record_node()` / `record_code()`, `set_status()`, `refresh_clustering_choices()`

### Tab Modules (`src/xenium_viewer/tabs/`)

Each tab follows a consistent pattern:

```python
def build_tab(ctx: ViewerContext) -> tuple[QWidget, dict]:
    # Build UI, register callbacks that access ctx
    # Return (tab_widget, exports_dict)
    # exports_dict may contain "restore_session" callable
```

The tabs are grouped under Cells / Genes / Spatial / Images / Tools. **Tools → Cache**
(`tabs/tab_cache.py`) exposes the cache health check and repair actions described
under "Cache safety" below: verify, re-consolidate, recover from a backup, and a
force rebuild that moves the old cache aside rather than deleting it.

**Tools → Dataset** (`tabs/tab_dataset.py`) is the per-item view the Cache tab's
whole-store report cannot give: a `QTreeWidget` of everything on disk with sizes, and
checkbox deletion of the parts the viewer created. See "Deleting components" below.

**Anything that changes the zarr behind the viewer's back must call
`ctx.reload_dataset()`** (bound by `app.py`, rebuilt on every dataset load). The live
`SpatialData`, the napari layers and every tab's widgets are built from disk once at load
time and never re-read it, so recovered or externally-written elements stay invisible
until the whole dataset is reloaded. It is synchronous and tears down every widget —
never call it from inside a worker that holds a reference to a tab.

**A `build_tab` must return a scroll-wrapped widget** — `make_tab()`, or `scrollable()`
for a tab that builds its own root. The panel is a `QTabWidget` of `QTabWidget`s, and a
stacked widget's minimum size is the **maximum over all its pages, hidden ones included**,
so one unwrapped page becomes the floor for the whole Xenium Controls dock and the
separator stops moving — even if nobody ever opens that page. A `QScrollArea` reports a
fixed ~68px minimum whatever it holds, which is what keeps the floor low. The same applies
to a button row placed *outside* a tab's scroll area: buttons cannot shrink below their
labels, so use `toolbar_row()`, which stays pinned but scrolls sideways.
`tests/test_control_panel_size.py` measures it; both rules were once broken (Notebook and
Templates, pinning the dock at 536×534). **The same rule applies to the Plots dock**,
which is a second dock with the same failure mode — `utils/plots_panel.py` wraps its
gallery in `scrollable()` and its buttons in `toolbar_row()` for exactly that reason.

`src/xenium_viewer/tabs/_helpers.py` contains shared utilities (e.g., `StatusProxy`, `make_tab()`).

### Figures (`utils/plot_output.py`, `utils/plots_panel.py`, `utils/fig_render.py`)

**Every figure goes through `ctx.show_plot(fig, stem, title=)`.** It saves under
`<data_path>/plots/` in each format the user chose (`plot_output.plot_formats`, default
PNG **and** PDF), appends the figure to the **Plots** dock, reveals the dock, and returns
the paths written. `save=False` is for a figure its own Step template already wrote —
the template must stay the thing that writes the file, or the recorded code stops being
the code that ran; pass the `paths` it used so the card can still say where it went.

Four rules, each of which was a real defect (`tests/test_plot_consistency.py` is the
source guard, in the idiom of `test_persistence_safety.py`):

- **No tab may call `plt.show`.** Display draws into the dock's own `FigureCanvasQTAgg`
  (`fig_render.open_figure_window`), so it does not depend on which backend pyplot
  happens to be pointed at. Eighteen sites used four different display modes, two of them
  *blocking*.
- **No tab may call `matplotlib.use`.** It is process-wide and was never restored:
  `tab_umap` and `tab_marker_genes` set `Agg` inside a worker, so after a user saved a
  UMAP or ran a marker plot, every later `plt.show` in the session became a silent no-op
  and figures from unrelated tabs stopped appearing. `cnv_copykat_worker.py` still calls
  it, correctly — a detached process with no GUI.
- **No `savefig` outside `plot_output.save_figure`** (the Plots dock's own "Save as…"
  excepted — that is a copy to a place the *user* picked). A figure factory returns a
  Figure and lets the caller write it; `render_dag` and `make_cnv_heatmap` both follow it.
- **`plots_dir(data_path)` is the one definition of `<data_path>/plots`.** Five modules
  built it by hand.

A recorded `plot:*` node must name the file that was **actually written**, via
`ctx.recorded_plot_paths(paths)` (relative to `data_path`) and `recorded_save_code`,
which also emits the `mkdir` a replayed notebook needs — `savefig` does not create its
parent. Every hand-written node used to record a bare `"dotplot.svg"` that matched
nothing on disk. Stems are keyed by clustering/gene (`safe_stem`) so a second run does
not silently overwrite the first.

Batch outputs (the three pairwise-volcano generators) stay out of the gallery — an N×N
run is fifty figures — but default their directory to `batch_dir(data_path, ...)` and
honour the format setting.

### Key Utilities (`src/xenium_viewer/utils/`)

| Module | Purpose |
|---|---|
| `coloring.py` | `CellColorManager` with `DirectLabelColormap` for O(nonzero) raster colorization |
| `gene_analysis.py` | Rank genes, normalization, Leiden clustering |
| `spatial_analysis.py` | Squidpy-based spatial analysis (neighborhood enrichment, co-occurrence, L-R) |
| `registration.py` | Landmark-based similarity affine registration for H&E/ARMS |
| `transcript_index.py` | Per-gene feather loader |
| `session.py` | Zarr-based session persistence (ROIs, H&E/ARMS registration, clusterings, DEG results, provenance graph) |
| `prov_graph.py` | Provenance DAG for reproducible code — nodes/deps, upsert+staleness, topo-sort, cells/script/mermaid/dot rendering |
| `step_templates/` | The text of every analysis template, as `builtin/*.tmpl`. `spec.py` (TemplateSpec/BlockSpec/ParamSpec), `loader.py` (parse + `builtin_*`), `namespace.py` (`EXECUTOR_BASE_NAMES`) |
| `notebook_export.py` | Build/write/read `.ipynb` (nbformat) from the graph |
| `dag_view.py` | Matplotlib+networkx render of the provenance DAG |
| `umap_widget.py` | Separate linked napari window for UMAP scatter |

### Session Persistence

Stored in `sdata_cached.zarr/viewer_session/` as zarr arrays and parquet files (plus the code provenance graph as a JSON attr). Auto-saved on relevant actions and restored at startup. Supports zarr v2 and v3.

### Cache safety (`utils/zarr_safe.py`, `utils/cache_repair.py`)

**Never write to the zarr store directly.** Every element write goes through
`safe_write_element` / `safe_delete_element`, and plain zarr groups through
`safe_group_update`. These stage the new value in `<cache>/.xv_staging/`, journal the
operation, then swap it in with two `os.rename` calls; the previous copy is moved to
`<cache>/.xv_trash/` rather than deleted. `delete_element_from_disk` unlinks the live
element *before* its replacement exists, so any interruption left the store invalid —
`tests/test_persistence_safety.py` fails if it is called outside `zarr_safe.py`.
`recover_pending()` finishes or unwinds an interrupted swap at startup.

`cache_repair.verify()` is read-only and parses the root `zarr.json` with `json.loads`
(not `zarr.open`), so it works on a store too broken to open; `repair()` only renames,
clears debris and re-consolidates. `loader.load_sdata` repairs before asking, and
**never deletes a cache** — every destructive branch is a rename (guarded by a test).

Cache freshness comes from `<cache>/.xv_manifest.json`, a sha256 of `experiment.xenium`.
Caches without one fall back to an mtime comparison, treated as *uncertain* — it prompts
rather than rebuilding, since copying a dataset changes mtime without changing content.

**Sidecar analysis outputs** (`adata_norm_cache.h5ad`, `adata_cnv_cache_*.h5ad`,
`roi_deg_cache.parquet`, `arms_tile_deg_cache.parquet`, `cnv_*_result.json`) live in
`<data_path>/viewer_cache/`, *not* in the zarr store root — files there make zarr's
hierarchy walk warn on every consolidation, and a rebuild would delete them. Write via
`adata_persistence.sidecar_write_path()`; read via `find_sidecar()` / `glob_sidecars()`,
which fall back to the legacy in-store location for existing datasets.

### Renaming a dataset (`scripts/rename_dataset.py`)

**Nothing stores a dataset name.** Not `adata`, not `sdata`, not the cache manifest;
the folder name reaches only `viewer.title`. Cache freshness is a content hash of
`experiment.xenium`, and every cache/sidecar path is re-derived from `data_path` at
launch — verified: no `zarr.json`/`.zattrs`/`.zgroup`/`.zmetadata` under a real store
contains an absolute path. So a renamed dataset opens fine and does not rebuild.

What a bare `mv` breaks is the **absolute paths recorded in the provenance graph** —
the `preamble` node's `data_path = Path(r"…")` and each `clustering:<key>` node's
`read_csv` — plus the session's `he_path`/`arms_*_path` attrs and the (write-only)
CopyKAT result JSONs. It does not stay quiet: `app.py` re-emits the preamble for the
current `data_path` on every launch, so `upsert` flags **every descendant stale** and
the first launch after a manual rename marks the whole notebook ⚠ for nothing.
`xenium-rename-dataset` repairs all of it before that launch; `--repair` infers the
old path from the recorded preamble.

Two rules carry the safety, both under test:
- **Substitution is path-prefix only** (`/d/foo` must not match `/d/foobar` or
  `/d/foo.bak`), which is what leaves a file living *outside* the dataset directory
  untouched — that file did not move.
- **Every write swaps an inode**: `safe_group_update` / `atomic_json` for the store,
  and `_atomic_replace` for `analysis.py` / the notebook. Writing those two in place
  mutated a `cp -al` snapshot of the dataset through the hardlink — found by probing
  the real tool against a real dataset, and pinned by
  `test_derived_outputs_are_replaced_not_written_in_place`.

### Deleting components (`utils/store_inventory.py`, `tabs/tab_dataset.py`)

`store_inventory` is the model behind Tools → Dataset: five `Section`s of `Node`s
(raw output / cache elements + table contents / session state / derived caches /
backups & trash), each with a size, a detail and an honest `recoverable` value.
Filesystem-only and Qt-free like `cache_repair.verify`, so it reports on a store too
broken to open — and **read-only**, guarded by a test that greps it for every mutating
call. Two rules carry the safety:

- **`assert_deletable(path, roots, kind=)` is the single choke point.** A path may be
  removed only if it resolves inside a `deletable_roots()` directory — `sdata_cached.zarr`,
  `viewer_cache/`, `transcript_cache/`, or a `sdata_cached_*_*.zarr` backup. The raw 10x
  output is in none of them. A root that *is* or *contains* the dataset directory is
  refused outright (a `cache_path` of `.` would otherwise make every raw file deletable),
  and symlinks are resolved before the containment test but never followed when sizing or
  deleting. `tests/test_store_inventory.py` asserts the **property** over every node the
  inventory produces, not a list of remembered cases.
- **Unrecognised defaults to not deletable.** An unknown entry in the dataset directory is
  raw output; an unknown element or obs column is left alone. The `loader._USER_*`
  allow-lists decide only the yes cases. `CORE_ELEMENTS` (`tables/table`, both label
  rasters, `morphology_focus`, `points/transcripts`) are listed with sizes but blocked;
  so are `prov_graph*.json` and the `prov_graph` session attr.

Table contents (`OBS`/`UNS`/`OBSM`) carry `path=None` on purpose: deleting a column is a
rewrite of the whole table, so there is no path for an executor to `unlink`.

**A clustering is more than one obs column.** A Leiden run leaves the bare `<key>` (the
recorded step's `adata.obs[$key] = …`, needed so the notebook reproduces it) *and*
`clustering_<key>` (`save_clustering_to_adata`), plus `cluster_labels_<key>` once clusters
are named. `_clustering_twin_of` pairs the bare column with its prefixed one so all of
them cascade together — otherwise "delete this clustering" left an identical copy. The
pairing requires the `clustering_<name>` column to exist and the bare name not to be in
`_STRUCTURAL_OBS`, so no Xenium column becomes selectable.

**In the tree, a blocked row is dimmed, never `setDisabled(True)`.** Qt propagates a
disabled item down its whole subtree, so disabling `group:tables` and the core
`tables/table` also greyed out every clustering inside them — i.e. exactly what the tab
exists to delete. `tests/test_tab_dataset.py` checks the *effective* state through
ancestors, because a row's own flags do not tell you whether a user can tick it.

The executor (`tab_dataset._apply_deletion`) applies by kind in `Plan` order — **table
edits first**, backups last. Three things it must not skip, each of which was a real bug:
`_persist_table` runs **once** per batch; a deleted `clustering_*` column is also popped
from `ctx.clusterings` / `state["custom_clusterings"]` / `state["cluster_labels"]`
(`refresh_clustering_choices` reads that dict, *not* `adata.obs`); and a deleted session
node clears its `store_inventory.SESSION_MEMORY` mirror, or `save_session` writes it back
at exit. Failures are per node and named in the report — a partly applied batch has to be
describable, since several of these are irreversible. `_remove_tree` is the only function
that touches the filesystem, enforced by a test that parses the module.

### Code Recording (provenance DAG)

User actions are recorded as a **provenance graph** (`utils/prov_graph.py`), the
single source of truth for reproducible code — `ctx.state["prov_graph"]`. Each step
is a node with a stable `id` (the artifact it produces), its `code`, its `deps`
(parent node ids), and a `kind` (`setup` / `artifact` / `terminal` / `note`). The
notebook is *derived* from the graph by topological sort, so it always respects
dependencies regardless of the order actions were taken — even across sessions.

**`note` is the "not code, and not meant to be" kind.** A node whose cell is a
comment replays as a silent no-op: `allow_errors=False` sees a cell that ran, so a
missing step and a viewer-state annotation looked identical to every consumer.
`NOTE` declares the second case — it renders as *markdown* in the notebook, keeps
its comment in `analysis.py`, is labelled in the Notebook tab, and
`scripts/verify_notebook.py` counts it apart from the comment-only punch list. Use
it only for state with no notebook equivalent (canvas background, overlays,
crop-export). A comment-only node of any other kind is a defect, and
`tests/test_recorded_code_is_code.py` parses every `record_node` call site and
fails on one — `viewer:transcript_density` is the single listed exception.

- **Preferred: `ctx.run_step(Step(...))`** (`utils/steps.py`). A `Step` is a node id, a
  `string.Template` of plain scverse source, and a dict of literal `params`. `run_step`
  renders the template **once** and hands that same string both to `exec` and to the
  provenance graph — so the code the GUI runs *is* the code the notebook records, by
  construction rather than by discipline. **The invariant to enforce in review: a tab
  callback may never call an analysis function with a widget value; it may only build a
  `params` dict.** Params are validated with `ast.literal_eval(repr(v)) == v` (use
  `steps.coerce()` at the widget boundary for numpy scalars); templates use `$name` so
  `{...}` literals survive; execution is serialised and proceeds statement-by-statement
  so progress can be reported without changing the compiled source. Failures raise
  `StepError` naming the step, and nothing is recorded for a step that did not succeed.

  **Template text lives in `utils/step_templates/builtin/*.tmpl`, not in the tab modules.**
  A template is an ordered dict of *named blocks*; the `.tmpl` header declares the contract
  (`params`, `requires`, `outputs`, `assemblies`, `frozen-blocks`) and everything structural
  is a comment, so the file is valid Python. **The call site owns which blocks are selected;
  the registry owns their text** — the branch structure *is* what the widgets mean
  (`use_hvg or do_scale` → the `pca` block), so selection stays in Python. The tab modules'
  private constants (`_leiden_template`, `_NORMALIZE_TEMPLATE`, …) still exist, bound through
  `builtin_text`/`builtin_assemble`, which read only shipped files via `importlib.resources`
  and cannot see an override path — that is what keeps the six template-pinning test modules
  passing unchanged. `tests/test_template_registry.py` runs the `check_step` lint over every
  template × every declared assembly (40 renderings), which is where the five hand-written
  `check_step` calls became a registry-wide gate.

  **Users can override a template**, per user, in `~/.config/xenium-viewer/templates/*.tmpl`
  — **resolved per block**, so blocks the user did not touch keep tracking the shipped
  template and still receive upstream fixes. `loader.resolve()` never raises and never
  returns nothing: an invalid override is skipped, the builtin is used, and the problems ride
  along on the `ResolvedTemplate` for the GUI to surface (plus a once-per-session napari
  warning via `reporting.report_template_rejected`). Call sites use
  `step_template(id, blocks)`, which returns the text **and** its provenance stamp together
  — a stamp fetched separately could describe a different resolution than the text it labels.
  `validate.py` is the gate; the check that matters most is that a **required param the
  template no longer mentions is a hard stop**, because that template runs, succeeds, and
  silently ignores the user's setting. Two off switches: `--no-user-templates` and
  `XENIUM_VIEWER_TEMPLATE_PATH` (emptied by `tests/conftest.py`, so a dev's own overrides
  never change what the suite asserts). Saving derives its destination from the same search
  path reading uses, so a write cannot land where the reader does not look.

  **Upgrades.** `overrides.json` records the hash of the *shipped* text each overridden block
  was forked from. The `dpkg` "unmodified conffile" case needs no logic — an untouched block
  is simply absent from the file, so it already tracks upstream. Only a block the user changed
  can conflict: `ResolvedTemplate.stale_blocks` / `.needs_review` flag it, the edit still
  applies, and the tab offers a two-way diff plus "Take new default" for the moved blocks only.
  A conflict that no longer validates is deactivated instead of flagged. Missing/corrupt
  manifest ⇒ not flagged, deliberately: prompting on every override after every upgrade with
  nothing specific to point at trains people to dismiss the warning.

  **Downstream visibility.** `notebook_export.customisation_banner` prepends one markdown cell
  when any node's `template_origin != builtin` (in `notebook_export.graph_to_cells`, *not*
  `prov_graph`'s, so `analysis.py` and the replay test's verbatim-code-cell property are both
  unaffected). `scripts/verify_notebook.py::template_provenance` adds a `templates` report
  section and a `stock_templates` bool.

  **Tools → Templates** (`tabs/tab_templates.py`) shows the contract, **Default (read-only)
  beside Yours (editable)**, a live preview of the exact string that would be `exec`'d, plus
  Validate / Save / Revert and a problems list. **Save never refuses**; activation is what is
  gated, and an invalid file is simply rejected by the resolver.

  **The preview is rendered from the owning tab's `Preview(blocks, params, note="")`**
  (`step_templates/spec.py`), registered in `ctx.state["template_preview"]` — and *the tab's
  own callback runs from the same call* (`_leiden_preview` in `tab_clustering.py` is the
  worked example). **Blocks travel with the params**: block selection lives at the call site
  by design, so a params-only provider left the pane pinned to `assemblies[0]` — the numbers
  tracked the widgets while the code shape did not. `note` names a value that cannot come
  from a widget yet (a save-dialog path renders as the filename the dialog would propose and
  the header says so), and a provider must stay **read-only** — no `makedirs`, no
  `record_clustering` — since drawing a pane must not have side effects.

  Twelve of fourteen templates have a provider. The two exemptions are declared in
  `tests/test_tab_templates.py::_NO_PROVIDER` with their reasons: `normalize` takes no params,
  and `spatial_neighbors` takes `k` from whichever tab called `ensure_spatial_neighbors`, so
  it uses the `# sample-params:` header field instead. Four gates, each verified against the
  defect it describes: every template has a provider or an exemption; every provider is
  *called* by its own tab; every provider **answers** when invoked against a live tab (a
  raising provider is otherwise invisible — `_preview` catches it and shows sample values,
  exactly as it should for a half-built tab); and every provider selects a declared assembly.

  Five rules that exist because each was once broken:
  **(a) A template may only reference `EXECUTOR_BASE_NAMES`** (`utils/step_templates/namespace.py`:
  `sc sq sd pd np plt Path data_path sdata adata`) plus names a declared dependency binds.
  `_get_executor` calls `check_base_namespace` on the dict it built, so the set validation
  checks against and the set execution provides cannot drift — a name in one and not the
  other passes validation and then fails as a `NameError` on replay, in a clean kernel.
  **(b) Results come back through declared `outputs`, not by reading `ctx.adata` afterwards.**
  Reading back worked only while the executor namespace and `ctx.adata` were the same object;
  `StepExecutor` raises if the template does not bind a declared output, so a template edit
  that stops producing a result fails loudly instead of returning stale state. Leiden binds
  `leiden_labels` for exactly this reason.
  **(c) No template may carry a `$token` that is not a declared param.**
  `Template.substitute` (not `safe_substitute`) is the only thing checking a template against
  its params, so stripping a fake token with `str.replace` before `Step` sees it hides the
  template from that check. Two did (`$n_suffix`, `$dpi_kwarg`); both are now whole-line block
  variants, and `tests/test_template_placeholders.py` is a source guard against the idiom
  returning.
  **(d) A run site gets its template from `step_template`, never `builtin_assemble`.**
  `builtin_*` cannot see an override path — which is what makes the pinning tests immune to a
  developer's own config, and what silently disabled customisation for `genes.marker_plot`,
  the one call site that used it. That template could be edited, validated and saved with no
  effect, and its nodes carried no `template_id` for the notebook banner or `stock_templates`
  to notice. `tests/test_tab_templates.py::test_every_step_resolves_user_overrides` parses
  every `Step(...)` in every tab; `tab_templates._preview` is the one exemption, since it
  renders a spec it has already resolved.
  **(e) A template may not hand-roll what a library API already does.** Templates are read as
  the explanation of the analysis, so prefer `scanpy` / `squidpy` / `spatialdata` / `anndata`
  calls over manual numpy/pandas/shapely equivalents — `sc.get.obs_df` over indexing `.X` and
  densifying, `sc.tl.rank_genes_groups(key_added=)` over a results dict, `sd.polygon_query`
  over a point-in-polygon loop, `shapely.make_valid` over `buffer(0)` (which silently *deletes*
  a lobe of a self-intersecting polygon rather than repairing it). Hand-rolled code is allowed
  where no API covers it — the Welch + BH block in `roi.expression`, the frozen `arrow_shim`
  block — and says so in a comment. **A manual coordinate transform is the specific smell**: a
  `spatialdata` transformation *declares* the frame instead of applying it by hand, which is
  what keeps the notebook's coordinate conventions honest. The ROI templates were the worked
  example: a `[:, ::-1]` flip undone by the `[:, 1], [:, 0]` read on the very next line, so the
  comment asserted a convention the code reversed. One caveat found while applying the rule:
  `sc.pp.calculate_qc_metrics(inplace=True)` overwrites Xenium's own `obs['total_counts']`,
  which `store_inventory._STRUCTURAL_OBS` lists — a template must not mutate `adata` to reach
  an API, so `genes.correlation` uses `inplace=False`.

  `ProvNode`/`Step` also carry `template_id` / `template_origin` / `template_hash`.
  `code` still records what ran, so replay is unaffected; the fields let a reader tell a stock
  run from a customised or hand-edited one. They are excluded from `upsert`'s staleness
  comparison on purpose — they describe where the *same* code came from.

  Migrated so far: **Leiden clustering** (`tab_clustering.py`), **normalize**
  (`ctx.ensure_normalized()`, which binds `adata_norm` and replaces the old
  `record_normalize` + `get_normalized_adata` pair), **rank genes**
  (`tab_gene_analysis.py`), **spatial neighbours**
  (`ctx.ensure_spatial_neighbors(k)`, which builds the graph on `adata_norm` and
  replaces `record_spatial_neighbors`), **neighbourhood enrichment**, **co-occurrence**,
  **ligand-receptor**, **marker-gene plots**, **gene correlation**, **ROI DEG +
  `rois`**, **ROI expression + its CSV export**, and **inferCNV**. Every expression-based step consumes `adata_norm` and
  declares `deps=["normalize"]` — never an implicit reliance on `adata` having been
  normalised in place, which is what made the DAG lie before. Call
  `ctx.ensure_normalized()` (idempotent) before `ctx.run_step()` in any such step.
  One documented exception: the `preamble` node records
  `xenium(data_path)` while the viewer reaches the same objects via the zarr cache.
  A second documented exception: **CopyKAT** (`cnv:copykat`) stays on `record_node`
  because it runs detached in the `xenium_viewer_copykat` env — no in-process step can be
  the code that ran, so its cell says in-line that it is a reconstruction. `run_cnv_pipeline`
  is now the CopyKAT path only; the inferCNV template must stay in sync with it
  (`tests/test_cnv_step.py` pins both).
  Still unmigrated: the **annotation-neighbourhood** and **annotation-distance** tabs.
  Their synthetic virtual cells are sampled from a napari shapes layer the notebook has no
  access to — resolving that needs E3's spatialdata shapes. Their *figures* are now
  recorded, as `viewer:annot_nhood_plot` / `viewer:annot_distance_plot`, and deliberately
  as **`NOTE`**: a `TERMINAL` calling `plt.gcf().savefig(...)` with no preceding plot call
  replays as a silent no-op writing an empty figure, which is worse than recording
  nothing. **Plot/export terminals** across the migrated tabs are still on `record_node` —
  that is fine where the cell is real code (`sc.pl.*`, `to_csv`); what is not fine is a
  terminal whose cell is prose, which the source guard above now catches.
  Also migrated: **UMAP plots** (`umap.plot`, both by gene and by cluster).
- **Legacy: `ctx.record_node(id, code, deps=..., kind=..., label=..., params=...)`**
  in tab callbacks — still used by the not-yet-migrated tabs, and the reason the recorded
  and executed code could drift. Re-running a step (same `id`) revises its node in place
  and flags descendants stale; a missing dependency errors at record time.
  `ctx.record_code(code, tag)` remains as a thin backward-compat shim.
  Helper recorders: `record_preamble` (`preamble`
  node, defines `data_path` — and calls `record_environment`, which emits the
  `environment` node: package versions as a comment block plus seeds and
  `sc.logging.print_header()`, sorting first and with **no dependents**, so a version
  change is readable without flagging every result stale), `record_clustering`
  (`clustering:<key>`), `record_spatial_neighbors`. Identity conventions: `clustering:<col>`, `rank_genes:<key>`,
  `nhood:<key>`, `cooccur:<key>`, `ligrec:<key>`, `annotation:<col>`, `rois`, `roi_deg`,
  `cnv:<backend>` (`cnv:infercnv` / `cnv:copykat`); terminals `plot:*`
  (incl. `plot:cnv_heatmap:<backend>:<key>`) / `export:*` / `he:*` / `arms:*`;
  notes `viewer:*` and `crop_export:*`.
  A recorder that fails (missing dep, cycle) degrades to appending the snippet and
  reports through `reporting.report_recording_failure` — logged with a traceback and
  surfaced as a napari warning naming the node, never `warnings.warn`.
- **Outputs**: a flat `analysis.py` (derived, stable filename) written live, and
  `analysis_notebook.ipynb` (via `utils/notebook_export.py` + nbformat) on session save /
  the Notebook tab's "Export .ipynb". The notebook is code-only and replays from the raw
  Xenium output.
- **Notebook tab** (`tabs/tab_notebook.py`) renders the graph as topo-ordered cells with a
  ⚠ stale badge, an editable free-form cell area, and a "Show DAG" button (`utils/dag_view.py`,
  matplotlib+networkx). `graph_to_mermaid` / `graph_to_dot` give diagram text.
- **Persistence**: the graph is written to `<data_path>/viewer_cache/prov_graph.json`
  **on every recorded step** (`_helpers._save_prov_graph`) *and* serialized into
  `sdata_cached.zarr/viewer_session/` by `save_session` (dataset switch / exit).
  The sidecar wins on load (`app._load_prov_graph_items`), because the attr is
  behind whenever the viewer is still open or was killed. Persist the graph
  wherever an artifact is persisted — the artifacts are written eagerly, and a
  store holding results whose code is missing is the failure this pairing exists
  to prevent.
- **Every clustering needs a `clustering:<key>` node**, because analysis tabs
  declare `deps=["clustering:<key>"]`. The *producer* records it (Leiden, CNV,
  Novae, import); `ctx.record_clustering()` is only a backstop, and it records a
  `read_csv` **only when that CSV actually exists**. A producer that persists a
  column via `save_clustering_to_adata` without recording a node is a test
  failure (`tests/test_clustering_recording.py`).

### Verifying the claim (notebook replay)

`run_step` makes the recorded code equal the executed code *by construction*. Whether
replaying it reproduces the result is a separate, empirical question — so it is run, at
two tiers:

- **`tests/test_notebook_replay.py`** — the CI gate, described under "Running the Viewer".
- **`scripts/verify_notebook.py <dataset> --out report.json`** — the real measurement.
  Reads the graph out of `<cache>/viewer_session` attrs (zarr only — no GUI, no napari,
  no SpatialData load), exports the notebook, executes it against the **raw Xenium
  output** with per-cell timing, and compares the replayed `adata.obs` against the
  `clustering_*` columns and `uns['rank_genes_groups']` the viewer persisted. The JSON
  report carries per-clustering ARI, cluster counts, top-N gene agreement, wall-clock,
  package versions — and **the ids of every comment-only node**, which execute fine and
  do nothing, so `allow_errors=False` can never catch them. That list is the remaining
  recording work (Phase 0.3), enumerated by measurement; `note_nodes` is reported beside
  it and is *not* work — those are the declared viewer-state nodes. `--dry-run` produces
  both in seconds without replaying.

Both go through `notebook_export.execute_notebook()`, which runs the notebook in a
throwaway kernelspec pointing at `sys.executable`. Do not switch it to the installed
`python3` kernelspec: on a conda box that belongs to whichever env registered it last —
routinely base, which has no scanpy — and the resulting failure would look like a
reproducibility defect rather than a kernel-discovery one.

## Key Dependencies

- **napari** + **PyQt5/qtpy** + **magicgui** — UI framework
- **scanpy** / **anndata** — single-cell analysis
- **squidpy** — spatial transcriptomics analysis
- **spatialdata** / **spatialdata_io** — spatial data container and Xenium loader.
  `spatialdata` is bound as **`sd`** in the template namespace (`EXECUTOR_BASE_NAMES`) and
  imported by the recorded notebook preamble, so a template can call `sd.polygon_query` and
  `sd.transformations.*` directly — see rule (e) above.
- **zarr** — caching and session persistence
- **dask** — lazy array loading
- **tifffile**, **opencv**, **scikit-image** — image processing
- **insitucnv** (the `insituCNV-copykat` fork, in `environment.yml`) + **infercnvpy** — CNV
  inference. `utils/cnv_analysis.py` drives both backends via `run_cnv_pipeline(..., backend=)`.
  The fork is installed from an *unpinned* git URL (tracks master) in `environment.yml`,
  `environment-copykat.yml`, and the `cnv` extra — the same distribution in all three,
  because the two conda envs share no site-packages. Keep the fork's dependency bounds
  loose: upper pins there make the main env's pip section unsolvable against `-e .`.
  **inferCNV** runs in the main env. **CopyKAT** needs **rpy2 + R 4.3 + the `copykat` R package**,
  whose stack requires **python 3.11** — incompatible with the main env's python 3.12. So CopyKAT
  runs in a **second conda env** (`environment-copykat.yml` → `xenium_viewer_copykat`): the viewer
  resolves that env's python (`_resolve_copykat_python` in `tabs/tab_cnv.py`; override via
  `XENIUM_COPYKAT_ENV`/`XENIUM_COPYKAT_PYTHON`) and launches the detached worker
  (`src/xenium_viewer/cnv_copykat_worker.py`) there, passing the viewer source on `PYTHONPATH`.
  The detached process survives the GUI closing; the GitHub-only copykat R package auto-installs
  on first run (`src/xenium_viewer/install_copykat.py::ensure_copykat_installed`).

## Pending upstream deprecations (act before upgrading)

- ✅ **`sq.gr.spatial_neighbors` (removed in squidpy 1.9) — done**, issue #19. Now
  `sq.gr.spatial_neighbors_knn`, in `utils/spatial_analysis.py` and the
  `spatial_neighbors` template, which had to change together (the template is what the
  notebook replays). Measured identical graph, so no saved result was invalidated.
  `docs/squidpy-spatial-neighbors-migration.md` keeps the measurement, and the four
  things the original assessment got wrong — chief among them that Co-occurrence was
  never a dependent: it calls `sq.gr.co_occurrence`, which computes its own radii.
- ✅ **napari drops the PyQt5 backend in fall 2026 — done**, issue #15.
  `environment.yml` is `pyqt6` + `napari>=0.8`, `pyproject.toml` is `PyQt6`, and all 98
  unscoped enum sites and 8 `.exec_()` calls are in their scoped/Qt6 form.
  `docs/pyqt6-migration.md` keeps the measurement and the three things the original
  assessment got wrong. Two are worth carrying here, because both would mislead again:
  **conda-forge ships PyQt6 as `pyqt6`, not as a 6.x of `pyqt`** (that package stops at
  5.15, so `pyqt=6` fails as if PyQt6 did not exist), and **`QT_API=pyside6 pytest` never
  smoke-tested strict enums** — PySide6 runs in forgiveness mode and accepts every
  unscoped form, with qtpy supplying `exec_` on top. The same qtpy shim
  (`enums_compat.promote_enums`) means the pre-migration code would have run unchanged
  under PyQt6, so the enum edits were hygiene, not the blocker; the dependency solve was.
  CI asserts the resolved backend rather than adding a second leg.

- **`ScaleBarOverlay.unit` is removed in napari 0.9.0** — **no action needed**, recorded
  so nobody re-derives it. The viewer no longer touches that attribute: it is already a
  deprecated no-op (its setter warns and does nothing; napari's own test asserts it reads
  back `None`). It matters only as the reason `utils/units.py` exists. The unit used to be
  a `pint.Quantity` on the overlay, so `scale_bar.unit = "0.2125 um"` carried its magnitude
  and converted the bar at display level with nothing else changed; it is now derived from
  the layers as a `pint.Unit` × 1, i.e. magnitude forced to one. That is why the pixel size
  has to live in `layer.scale`, and therefore why every stored pixel affine needs converting
  at the napari boundary. `tests/test_units.py` pins the removed levers — **when 0.9.0 drops
  the attribute, re-check whether a display-level route has come back**, because if it has,
  `utils/units.py` and its ~11 call sites can be deleted rather than maintained.

## Known Compatibility Patches

- **ICE/X11 disconnect** — handled at startup of `src/xenium_viewer/app.py`, gated on
  `sys.platform.startswith('linux')` (macOS has no session manager, so clearing
  `SESSION_MANAGER` there would be an unexplained edit to the user's environment).
- **Missing `libglx-devel`** — `utils/gl_check.py`, called from `app.py` *before*
  `import napari`, because that import is what aborts. It lives in its own module so it can
  be tested without importing napari. It only reports; it does **not** repair: preloading the
  env's `libGLX.so.0` with `RTLD_GLOBAL` does not work, because PyOpenGL still `dlopen`s the
  host's *unversioned* `libGLX.so` as a separate mapping — only the unversioned name existing
  inside the env fixes it. It fires only when both halves of the collision are present (env
  lacks the name, host has a copy); warning on the missing package alone would fire on every
  correctly-working box with no `libglx-dev` installed. **Both checks are for the unversioned
  name**: `ctypes.util.find_library('GLX')` returns `libGLX.so.0`, which every working box
  has, so it cannot answer this question — hence the globbed path list.
- **pandas 3.0 PyArrow strings** — `_convert_arrow_strings()` in `src/xenium_viewer/loader.py`
- **NumPy 2.0** — `np.NAN` fallback for omnipath compatibility
- **matplotlib 3.9 `cm.get_cmap` removal** — `_patch_matplotlib_cm_compat()` in
  `src/xenium_viewer/utils/cnv_analysis.py`. Now a no-op against the pinned
  `insituCNV-copykat` fork (fixed there); retained for pre-existing environments
  and upstream InSituCNV.

## Version History

See `CHANGELOG.md`. The codebase was refactored from a 4295-line monolith into 11 modular tabs in March 2026.
