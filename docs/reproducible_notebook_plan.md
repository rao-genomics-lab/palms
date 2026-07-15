# Reproducible Notebook Capture for Xenium Viewer — Provenance DAG

## Context

Today, every viewer session writes its own throwaway `code_<timestamp>.py` in the data
dir. A workflow that spans multiple sessions (clustering on Monday, neighborhood analysis
Thursday, re-clustering next week) is scattered across partial files, none complete.
Worse, much of the recorded code does not actually run: plot-save snippets reference an
undefined `fig`, ROI/ARMS DEG reference an undefined `data_path`, the preamble is only
emitted on clustering-derived paths, and whole tabs (Cell Coloring, Transcripts, ROI
polygon drawing, annotations) record only `#` comments instead of executable code.

The goal: a **single, durable, accumulating Jupyter notebook** that captures the full
analysis and **replays from the raw Xenium output**, so a user can open it in Jupyter/
VSCode outside the viewer and continue their analysis. Confirmed decisions:
- **Replay from raw** — the notebook re-runs every step from `xenium(<path>)`; it does not
  depend on `sdata_cached.zarr` state, so completeness + determinism are mandatory.
- **Persist & accumulate** across sessions into one coherent notebook.
- **Fix reproducibility** so the exported notebook runs top-to-bottom.
- **Model provenance as a DAG (chosen architecture).** Rather than an append-only journal,
  record the analysis as a dependency graph of artifacts. The notebook is *derived* from
  the graph, not from insertion order. This was chosen over the simpler append/upsert
  journal because it catches problems earlier (a missing dependency is an error at record
  time, not a `NameError` at replay time), guarantees correct ordering, and enables
  automatic staleness detection — the "re-run clustering in session 3" problem is solved
  by construction. Accepted cost: more upfront infrastructure and a per-call-site burden
  to declare dependencies.

## What already exists (reuse, don't rebuild)

- **In-app Notebook tab** — `src/xenium_viewer/tabs/tab_notebook.py` renders the code
  journal as live, editable Jupyter-style cells; `record_code()` auto-appends via the
  `state["_notebook_sync_fn"]` hook. Syntax highlighting + inline output.
- **Execution engine** — `src/xenium_viewer/utils/notebook_engine.py` runs cells in
  napari's in-process IPython kernel (`adata`/`sdata`/`ctx` in scope in-app).
- **Recorders** — `_helpers.py:238-318`: `record_code(code, tag)` plus `record_preamble`,
  `record_normalize`, `record_clustering`, `record_spatial_neighbors`. The
  `record_clustering → record_normalize → record_preamble` chain is already an implicit
  hardcoded dependency path — it becomes the first real edges in the graph.
- **"Continue from existing code file…"** menu (`_helpers.py:692-725`) — manual prototype
  of cross-session persistence, superseded by graph persistence.

## Architecture: the provenance DAG

### Node model

Replace `code_journal: list[str]` + `code_journal_tags: set` with a **graph**:
`dict[node_id, ProvNode]`. A `ProvNode` holds:

- **`id`** — stable identity of the *artifact produced* (see scheme below). Upsert key.
- **`code`** — self-contained snippet that regenerates the artifact given its deps.
- **`deps`** — list of parent node ids it consumes.
- **`kind`** — `setup` (preamble/normalize), `artifact` (reusable state: clustering,
  neighbors, DEG, nhood), or `terminal` (side-effect only: plot/export/viewer-only; no
  dependents, prunable).
- **`meta`** — optional params dict, session id / timestamp, human `label` (for markdown
  headers), and a `stale` flag.

Edges are implicit from `deps`. The graph is acyclic (validate on insert).

### Recording API

New `ctx.record_node(id, code, deps=(), kind="artifact", label=None, params=None)`
(replaces `record_code`; keep a thin `record_code` shim during migration):
- **Upsert by id** — if `id` is new, insert; if it exists with different `code`, replace
  (revise) and propagate staleness to descendants; if identical, no-op.
- **Validate deps** — every dep id must already be registerable; a missing/unknown dep is
  surfaced as an error at record time (the early-catch benefit that motivated the DAG).

Rewrite the helper recorders as node registrations that declare deps:
- `record_preamble()` → `preamble`, deps `[]`.
- `record_normalize()` → `normalize`, deps `[preamble]`.
- clustering → `clustering:<obs_column>`, deps `[normalize]` (or a custom-preprocess node).
- expression neighbors → `neighbors`, deps `[normalize]`; spatial neighbors →
  `spatial_neighbors`, deps `[preamble]`.
- nhood → `nhood:<cluster_key>`, deps `[clustering:<cluster_key>, spatial_neighbors]`;
  co-occurrence → `cooccur:<cluster_key>`; ligrec → `ligrec:<cluster_key>`.
- rank genes → `rank_genes:<groupby>`, deps `[clustering:<groupby>]`.

### Notebook derivation (`graph_to_notebook`)

1. **Topological sort** of the graph, respecting deps; tie-break by (`kind` order: setup →
   artifact → terminal, then insertion order) for stable, readable output. This is what
   guarantees preamble → normalize → clustering → neighbors → nhood regardless of the
   wall-clock order actions were taken across sessions.
2. Emit one code cell per node; optional markdown header from `label`.
3. **Prune** (setting): drop `terminal` nodes and/or graph branches not reachable from a
   "kept" output, so abandoned experiments don't bloat the notebook.
4. Code-only cells, no stored outputs (they regenerate on execution).

### Staleness

On upsert-with-changed-code, walk descendants (reverse edges) and set `meta.stale=True`.
Surface as a badge in the Notebook tab ("⚠ input changed — recompute"); advisory, does not
block. This is the capability the append/upsert model could not provide.

### Identity scheme (the crux — get this right)

Identity = the artifact produced; **all parameters live in the code body, never in the id**,
so re-running the same artifact with new params is a revise (upsert), not a silent drop or
duplicate. Conventions:
- Clustering: `clustering:<obs_column>` (e.g. `clustering:leiden_r1.0`). Note current Leiden
  key `leiden_r{resolution}` (`tab_clustering.py:48`) encodes only resolution — keep the
  column as identity but ensure changing n_neighbors/n_pcs/HVG/scale revises the same node.
- `neighbors`, `spatial_neighbors`, `rank_genes:<groupby>`, `nhood:<key>`,
  `cooccur:<key>`, `ligrec:<key>`, `rois`, `roi_deg:<group_def>`,
  `annotation:<src_col>:<tgt_col>`, `cnv:<params>`.
- Terminals: `plot:<kind>:<subject>`, `export:<kind>:<subject>` — dep on the artifact they
  consume, `kind="terminal"`.

### Non-artifact actions

Plots, CSV exports, and viewer-only actions register as `terminal` nodes depending on the
artifact they consume, with no dependents. They sort next to their producer and can be
filtered from export. Viewer-only visual state (cell coloring, background, transcript
overlay) → either a real standalone plotting terminal (e.g. `sq.pl.spatial_scatter`) or
omitted; not on the replay critical path.

## Implementation steps

### 1. New `src/xenium_viewer/utils/prov_graph.py`
`ProvNode` dataclass; graph container with `upsert(node)`, cycle check, `mark_stale`
(descendant propagation), `topo_sort()`, `prune()`, and `graph_to_notebook(graph) ->
nbformat.NotebookNode`. Pure, unit-testable module with no Qt/napari deps.

### 2. Rewrite recorder API in `_helpers.py:238-318`
Replace `record_code`/`record_preamble`/`record_normalize`/`record_clustering`/
`record_spatial_neighbors` with `record_node`-based versions that declare deps and upsert
into the graph. Extend the preamble (`_helpers.py:261`) with `from pathlib import Path` and
`data_path = Path("<ctx.data_path>")` (fixes ROI/ARMS DEG). Emit the `preamble` node at
startup when recording is enabled. Keep a `record_code(code, tag)` compatibility shim that
maps to a `terminal`/opaque node so un-migrated call sites keep working during rollout.

### 3. Migrate all call sites + fold in reproducibility fixes
Migrate each `ctx.record_code(str, tag)` (~50 sites) to `record_node(id, code, deps, kind)`.
As each is migrated, fix its reproducibility defect so the emitted node is executable:
- **Undefined `fig`** → `fig = plt.gcf()` before `fig.savefig`: `tab_nhood.py:179`,
  `tab_co_occurrence.py:219`, `tab_ligrec.py:210`, `tab_gene_analysis.py:221`.
- **Determinism** → thread the actual `random_state`/seed into Leiden/neighbors/UMAP/
  `sc.tl.ingest` (`tab_clustering.py:84/87`, UMAP). Document residual scanpy/leidenalg
  version dependence.
- **Comment-only → real code:** ROI polygons → inline vertex arrays (mirror landmark
  pattern `tab_he_registration.py:353`); cluster label maps → `adata.obs[k] =
  adata.obs[src].map({...}).astype('category')` (`_helpers.py:476`); annotations
  (CellTypist/LLM/label-transfer, `tab_gene_analysis.py:427/493/676`) → real `.map()`;
  clustering import/export (`tab_clustering.py:187/216`) → real `pd.read_csv`/`to_csv`.
Representative tabs: `tab_clustering`, `tab_nhood`, `tab_co_occurrence`, `tab_ligrec`,
`tab_gene_analysis`, `tab_roi`, `tab_arms`, `tab_he_registration`, `tab_cnv`, `tab_novae`,
`tab_transcripts`, `tab_cell_coloring`.

### 4. Persistence (`session.py`)
Serialize the graph (nodes: id, code, deps, kind, meta) as JSON under the existing
`viewer_session` group in `save_session()`/`load_session()`. At startup in `_do_full_init`
(`app.py`), load the graph and seed `ctx.state["prov_graph"]`; new actions upsert into it.
Drop the per-session timestamped `code_file` (`app.py:834-837`); write stable `analysis.py`
+ `analysis_notebook.ipynb` sidecars. Cross-session accumulation + correct ordering is now
automatic — the session-3 re-run lands in its correct graph position, not appended last.

### 5. Notebook tab (`tab_notebook.py`)
Render cells from the **derived topo-ordered** cell list (via `graph_to_notebook`) instead
of the raw journal; keep the `_notebook_sync_fn` hook but have it re-derive on graph change.
Add stale badges, an "Export .ipynb" button, and a real exports dict with `restore_session`
+ `get_cells()` (currently returns `{}`). Reconcile manual cell edits back into the graph on
save/export so user-authored cells survive.

### 6. `.ipynb` export (`utils/notebook_export.py`, new)
`write_notebook(nb, path)` / `read_notebook(path)` wrappers around nbformat; the notebook
body comes from `graph_to_notebook`. Also written on session save. Add `nbformat` to
`environment.yml`.

### 7. Docs
Update `CLAUDE.md` (Code Recording section) and `CHANGELOG.md` to describe the provenance
graph, derived notebook, `.ipynb` export, staleness, and standalone replay (requires the
same conda env). Keep this file in sync.

## Files

| File | Change |
|---|---|
| `src/xenium_viewer/utils/prov_graph.py` | **New** — ProvNode, graph, upsert/staleness/topo-sort/prune, `graph_to_notebook` |
| `src/xenium_viewer/utils/notebook_export.py` | **New** — nbformat read/write around derived cells |
| `src/xenium_viewer/tabs/_helpers.py` | Replace recorder family with `record_node`; preamble `Path`/`data_path`; compat shim |
| `src/xenium_viewer/app.py` | Init/seed graph from session; emit preamble node; stable `analysis.py`/`.ipynb` |
| `src/xenium_viewer/utils/session.py` | Serialize/deserialize graph JSON in `viewer_session` |
| `src/xenium_viewer/tabs/tab_notebook.py` | Render derived cells; stale badges; Export .ipynb; `restore_session`+`get_cells` |
| `tab_clustering, tab_nhood, tab_co_occurrence, tab_ligrec, tab_gene_analysis, tab_roi, tab_arms, tab_he_registration, tab_cnv, tab_novae, tab_transcripts, tab_cell_coloring` | Migrate call sites to `record_node`; fix `fig`/determinism/comment-only |
| `environment.yml` | Add `nbformat` |
| `CLAUDE.md`, `CHANGELOG.md`, `docs/reproducible_notebook_plan.md` | Docs |

## Rollout / risk

- The `record_code` compat shim lets migration proceed tab-by-tab without breaking
  recording; un-migrated sites just produce opaque terminal nodes until converted.
- `prov_graph.py` is pure Python — build and unit-test the graph, topo-sort, staleness, and
  `graph_to_notebook` in isolation before wiring into the UI.
- Highest-risk piece is the identity scheme; nail down the id conventions (table above)
  first, since ids are the upsert keys and drive correctness.

## Verification

1. **Single-session capture.** Run clustering → rank genes → nhood → draw ROI → annotate.
   Confirm each registers a node with correct `deps`, the Notebook tab shows topo-ordered
   cells, and `analysis.py` + `analysis_notebook.ipynb` are written.
2. **3-session acid test (the motivating case).** Session 1: cluster. Session 2: nhood.
   Session 3: re-cluster *same column, new params* → confirm the clustering node is
   **revised in place** (not dropped, not appended after nhood), `nhood` is flagged
   **stale**, and the derived notebook stays correctly ordered and self-consistent. Then
   re-cluster under a *new* key → confirm a new independent node and that the unused branch
   is prunable from export.
3. **Standalone replay.** `jupyter nbconvert --to notebook --execute analysis_notebook.ipynb`
   in a clean env **outside the viewer** against the raw dataset — runs top-to-bottom with no
   `NameError` (esp. the previously-broken `fig.savefig`/`data_path` lines) and reproduces
   the clustering / DEG / nhood.
4. **Ordering invariance.** Perform the same steps in different wall-clock orders across
   runs; confirm the exported notebook is identical (topo sort makes order deterministic).
5. **Determinism spot-check.** Recorded Leiden/UMAP carry the seed; replayed labels match
   in-viewer labels (modulo documented version dependence).
