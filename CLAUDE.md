# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Xenium Viewer** is a napari-based spatial transcriptomics viewer for Xenium 3.x output data — a Linux equivalent to the commercial Xenium Explorer. It visualizes high-resolution spatial gene expression data with cell-level resolution.

## Running the Viewer

```bash
# First-time setup
conda env create -f environment.yml
conda activate xenium_viewer

# Optional one-time transcript preprocessing per dataset (~30-60 min)
# Without this, transcript loading falls back to scanning transcripts.parquet (~5s/gene instead of <100ms)
xenium-preprocess /path/to/xenium/output/

# Launch viewer (file dialog opens if no path given)
xenium-viewer [/path/to/xenium/output/]

# Launch without SpatialData zarr cache
xenium-viewer /path/to/xenium/output/ --no-cache
```

The package is installed as `xenium-viewer` (PyPI name) / `xenium_viewer` (import name) via `pip install -e .` (handled automatically by `environment.yml`). Console scripts: `xenium-viewer`, `xenium-preprocess`, `xenium-fetch-references`, `xenium-build-custom-segmentation`. You can also run `python -m xenium_viewer ...`.

There is a small `pytest` suite in `tests/` covering the codebase's *pure* logic
(provenance graph, CopyKAT subsampling, registration math, LLM prompt/response parsing,
patch-size inference, notebook export). Run it with `pytest` from the repo root
(`[tool.pytest.ini_options]` sets `pythonpath = ["src"]`, so no install is needed).
GitHub Actions (`.github/workflows/ci.yml`) runs the suite in the full conda env plus a
fast `ruff` error-only lint gate on every push/PR. The GUI, spatial-analysis, and
zarr/SpatialData persistence paths have no automated coverage — that testing remains
manual/exploratory.

## Architecture

### Entry Points & Load Sequence

1. **`src/xenium_viewer/app.py`** — Main entry point (~1300 lines). Validates data dir, orchestrates loading, builds napari viewer with layers, instantiates all managers, creates `ViewerContext`, builds all tab widgets, then restores session. The `main()` function is the `xenium-viewer` console-script entry point.
2. **`src/xenium_viewer/loader.py`** — Loads SpatialData from Xenium output. Uses a zarr cache (`sdata_cached.zarr/`) for 60–70% faster subsequent launches; staleness detected via `experiment.xenium` mtime. Public API: `load_sdata`, `load_umap`, `load_clusterings`, `get_label_to_obs_mapping`.
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

`src/xenium_viewer/tabs/_helpers.py` contains shared utilities (e.g., `StatusProxy`, `make_tab()`).

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

### Code Recording (provenance DAG)

User actions are recorded as a **provenance graph** (`utils/prov_graph.py`), the
single source of truth for reproducible code — `ctx.state["prov_graph"]`. Each step
is a node with a stable `id` (the artifact it produces), its `code`, its `deps`
(parent node ids), and a `kind` (`setup` / `artifact` / `terminal`). The notebook is
*derived* from the graph by topological sort, so it always respects dependencies
regardless of the order actions were taken — even across sessions.

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
  Migrated so far: **Leiden clustering** (`tab_clustering.py`), **normalize**
  (`ctx.ensure_normalized()`, which binds `adata_norm` and replaces the old
  `record_normalize` + `get_normalized_adata` pair), **rank genes**
  (`tab_gene_analysis.py`), **spatial neighbours**
  (`ctx.ensure_spatial_neighbors(k)`, which builds the graph on `adata_norm` and
  replaces `record_spatial_neighbors`), **neighbourhood enrichment**, **co-occurrence**,
  **ligand-receptor**, **marker-gene plots**, **gene correlation**, **ROI DEG +
  `rois`**, and **inferCNV**. Every expression-based step consumes `adata_norm` and
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
  Still unmigrated: the **annotation-neighbourhood** tab (records nothing; its synthetic
  virtual cells are sampled from a napari shapes layer the notebook has no access to —
  resolving that needs E3's spatialdata shapes). **Plot/export terminals** across the
  migrated tabs are still on `record_node`; the terminal-node policy is E4.
- **Legacy: `ctx.record_node(id, code, deps=..., kind=..., label=..., params=...)`**
  in tab callbacks — still used by the not-yet-migrated tabs, and the reason the recorded
  and executed code could drift. Re-running a step (same `id`) revises its node in place
  and flags descendants stale; a missing dependency errors at record time.
  `ctx.record_code(code, tag)` remains as a thin backward-compat shim.
  Helper recorders: `record_preamble` (`preamble`
  node, defines `data_path`), `record_clustering` (`clustering:<key>`),
  `record_spatial_neighbors`. Identity conventions: `clustering:<col>`, `rank_genes:<key>`,
  `nhood:<key>`, `cooccur:<key>`, `ligrec:<key>`, `annotation:<col>`, `rois`, `roi_deg`,
  `cnv:<backend>` (`cnv:infercnv` / `cnv:copykat`); terminals `plot:*`
  (incl. `plot:cnv_heatmap:<backend>:<key>`) / `export:*` / `viewer:*` / `he:*` / `arms:*`.
- **Outputs**: a flat `analysis.py` (derived, stable filename) written live, and
  `analysis_notebook.ipynb` (via `utils/notebook_export.py` + nbformat) on session save /
  the Notebook tab's "Export .ipynb". The notebook is code-only and replays from the raw
  Xenium output.
- **Notebook tab** (`tabs/tab_notebook.py`) renders the graph as topo-ordered cells with a
  ⚠ stale badge, an editable free-form cell area, and a "Show DAG" button (`utils/dag_view.py`,
  matplotlib+networkx). `graph_to_mermaid` / `graph_to_dot` give diagram text.
- **Persistence**: the graph is serialized into `sdata_cached.zarr/viewer_session/` and
  restored at startup, so a multi-session analysis accumulates into one notebook.

## Key Dependencies

- **napari** + **PyQt5/qtpy** + **magicgui** — UI framework
- **scanpy** / **anndata** — single-cell analysis
- **squidpy** — spatial transcriptomics analysis
- **spatialdata** / **spatialdata_io** — spatial data container and Xenium loader
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

## Known Compatibility Patches

- **ICE/X11 disconnect** — handled at startup of `src/xenium_viewer/app.py`
- **pandas 3.0 PyArrow strings** — `_convert_arrow_strings()` in `src/xenium_viewer/loader.py`
- **NumPy 2.0** — `np.NAN` fallback for omnipath compatibility
- **matplotlib 3.9 `cm.get_cmap` removal** — `_patch_matplotlib_cm_compat()` in
  `src/xenium_viewer/utils/cnv_analysis.py`. Now a no-op against the pinned
  `insituCNV-copykat` fork (fixed there); retained for pre-existing environments
  and upstream InSituCNV.

## Version History

See `CHANGELOG.md`. The codebase was refactored from a 4295-line monolith into 11 modular tabs in March 2026.
