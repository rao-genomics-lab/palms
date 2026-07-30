# Changelog

## [Unreleased] — 2026-07-30

### Added
- **Tools → Templates: see what an analysis button will run, before running it.**
  The source a step executes was only ever recoverable *after* the fact, from the Notebook
  tab. The text itself lived in 14 private module constants — seven of them assembled by a
  private function keyed on booleans — so there was no way to ask what a step takes, what it
  binds, or what it would run with the settings currently on screen.

  Template text now lives in `utils/step_templates/builtin/*.tmpl`: a header declaring the
  contract (`params`, `requires`, `outputs`, `assemblies`) and one or more named blocks.
  Everything structural is a comment, so a `.tmpl` is valid Python and the Notebook tab's
  syntax highlighter works on it unchanged. **The call site still owns which blocks are
  selected** — the branch structure is what the widgets mean — while the registry owns their
  text.

  The new tab shows the contract, the shipped source per block, and a live preview of the
  exact string that would be `exec`'d. The preview is not a reconstruction: it calls
  `Step.render()`, the same method the executor calls, with the owning tab's real widget
  values where the tab registers a provider. Leiden's `_leiden_params()` is now the single
  expression of "the current settings" that both the run and the preview read.

  The migration is provably text-preserving: all 47 template variants were compared against
  the previous constants and 43 are byte-identical, the other 4 differing only by the
  `leiden_labels` line added below. `tests/test_template_registry.py` then applies the
  `check_step` lint — which had existed since the Step system landed but was only ever called
  from five hand-written tests — to every template in every declared assembly, 40 renderings:
  each must parse, read only names the executor guarantees or the template declares, bind every
  output it claims, and never reach back into `xenium_viewer`. It immediately found one real
  contract error: `markers` was declared required, but `sc.pl.correlation_matrix` ignores
  `var_names`, so that assembly never referenced it.

### Changed
- **Groundwork for user-configurable analysis templates (Phase 0).** No user-visible
  feature yet; four preparatory changes that each stand on their own, ahead of moving the
  14 private template constants into a registry.

  - `EXECUTOR_BASE_NAMES` (`utils/step_templates/namespace.py`) now declares the names a
    template may reach for without binding them. `_get_executor` builds its dict and
    validates it against that declaration, so the set template validation checks against
    and the set the executor actually provides cannot drift — a name added to one and not
    the other would pass validation and then fail as a `NameError` on replay, in a clean
    kernel, long after the fact.

  - **Two templates no longer smuggle a fake `$token` past `Template.substitute`.** The
    gene-correlation tail carried `$n_suffix` and the marker-plot tail `$dpi_kwarg`, each
    stripped by `str.replace` *before* `Step` saw the text. Neither could ever become a
    real param — `$n_suffix` sits inside an f-string, where `repr('')` renders as `''` and
    breaks the literal — so they were load-bearing punctuation that looked exactly like
    parameters, and they hid those templates from the one check that validates a template
    against its declared params. Both are now whole-line block variants. All 26 assembled
    variants are byte-identical to before; a source guard fails if the idiom returns.

  - **The Leiden step declares its output.** It previously declared none, and the tab read
    the labels back off `ctx.adata.obs`, which worked only while the executor namespace and
    `ctx.adata` were the same object — an invariant maintained by hand and invisible to
    anyone editing the template. The template now binds `leiden_labels` and the tab reads
    the returned dict, so a template that stops producing labels raises `StepError` instead
    of silently handing back whatever obs column was left from a previous run.

  - `ProvNode` and `Step` gained `template_id` / `template_origin` / `template_hash`.
    `code` still records what actually ran, so replay is unaffected; the fields let a
    *reader* tell a stock run from a customised one, which rendered source cannot show.
    Old `prov_graph.json` sidecars and zarr session attrs load unchanged (every field but
    `id`/`code` is read with a default) and register as `builtin`. Template metadata is
    deliberately excluded from the staleness comparison — it describes where the same code
    came from, so changing it must not flag downstream results.

### Fixed
- **`pytest` hung forever on any machine with a display.** The suite stalled in
  `test_persistence_safety.py::test_a_failed_persist_leaves_the_table_readable` and never
  returned. `reporting._headless()` decided "nobody can answer a dialog" from
  `QT_QPA_PLATFORM` alone — but `conftest.py` deliberately leaves that unset when there is
  a `DISPLAY`, so on a developer desktop every signal it looked at was clear while the
  process still had no Qt event loop. `QMessageBox.exec_()` starts a *nested* event loop
  and returns when something dismisses the dialog; with no window manager and no user,
  nothing ever did. It needed any earlier test to have created the `QApplication`, which is
  why the file passed in isolation and why CI (no `DISPLAY` → `offscreen` → modal already
  suppressed) never saw it. The bug outlived the `_headless()` docstring that describes
  exactly this failure.

  `_headless()` is now the OR of two independent signals: the platform check, plus
  `_event_loop_running()` — `QThread.loopLevel() > 0` on the main thread, which is readable
  from a worker and is the direct answer to the question that actually matters ("is anyone
  processing events", not "is there a screen"). An unexpected Qt error assumes a loop *is*
  running: a needless dialog beats silence about a save that did not happen. Three tests
  cover it, including one that runs a real `exec_()` to prove the check does not simply
  always say "no loop" — which would have suppressed every modal in the real GUI.

- **The Notebook tab presented hand-edited cells as the recorded provenance.**
  `_reconcile_edits` folded a user's edit of a graph cell into `node.code` at export time
  **without executing it**, so the exported notebook showed never-run source in the place
  where every other node guarantees the code that produced the artifact. The node is now
  marked `hand-edited` and its `template_hash` cleared. The logic moved to a module-level
  `reconcile_edits(graph, cells)` so it is tested against a real graph rather than
  reproduced in the test.

- **CI could never run a Qt test: `pytest` aborted with a core dump.** The workflow set
  `MPLBACKEND` but not `QT_QPA_PLATFORM`, and on a runner with no display Qt does not fail
  when it cannot load the `xcb` platform plugin — it calls `abort()`. The whole run died at
  the first use of the `qapp` fixture with no test name and no traceback, only
  `Aborted (core dumped)`. Present since `test_tab_cache.py` was added and invisible
  because CI runs on PRs and pushes to `main`, and this work has been on feature branches.

  `tests/conftest.py` now chooses `offscreen` (and `Agg`) itself whenever there is no
  `DISPLAY`, so bare `pytest` works headless for CI, ssh sessions and containers alike; an
  explicit `QT_QPA_PLATFORM` or a real display still wins. `ci.yml` also sets it, to keep
  the runner honest about what it needs. The documented
  `env -u DISPLAY QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg` incantation is no longer
  required — a requirement that has to be remembered by hand is one CI will forget.

### Added
- **Tools → Dataset: see what the dataset holds on disk, and delete the parts the
  viewer created.** A dataset accumulates viewer data in four places nobody can see —
  `sdata_cached.zarr/`, `viewer_cache/`, `transcript_cache/` and sibling backup stores.
  Until now the only visibility was a comma-joined `Elements (N): …` line in the Cache
  tab's report, with no sizes at all, and the only way to remove anything was to delete
  a registered image or patch overlay from the tab that created it. There was no way to
  drop a clustering you no longer wanted.

  The new tab is a tree with a size per row: the original 10x output (read-only), every
  cache element, every obs/uns/obsm key inside the tables, session state, the derived
  caches, and the backups and trash — where the gigabytes usually are. Ticked rows are
  deleted through one executor, with a confirmation listing every path, the bytes
  reclaimed and a `⚠ not recoverable` block. On the reference dataset it scans in 2 s and
  its section totals match `du`.

  The safety property is structural, not a promise: `store_inventory.deletable_roots()`
  names the four directories the viewer created, `assert_deletable` refuses anything that
  does not resolve inside one of them, and a test asserts that over *every* node the
  inventory produces. Structural elements (`tables/table`, both label rasters,
  `morphology_focus`, `points/transcripts`) are listed with their sizes but cannot be
  selected — deleting the table bricks the dataset and the others break Crop Export and
  Segmentation-revert. Anything unrecognised defaults to not deletable, so an unfamiliar
  vendor file shows up read-only instead of becoming selectable. `prov_graph.json` and its
  dated backups are blocked: the sidecar wins over the session attr on load, so deleting
  it would silently lose every step recorded since the last save.

  Deleting a clustering also drops it from `ctx.clusterings`, which is what the combos
  actually read — `refresh_clustering_choices` never looks at `adata.obs`, so without that
  the column was gone from disk and still colouring cells. Deleting session state clears
  its in-memory mirror, or `save_session` writes it straight back at exit. `transcript_cache/`
  is offered (a `xenium-preprocess` re-run brings it back) and says so on the row.

  **Deleting a clustering means the whole clustering.** A Leiden run leaves *two* obs
  columns of the same data — the recorded step writes `adata.obs[$key]` so the notebook
  reproduces it, and `save_clustering_to_adata` writes `clustering_<key>` for the viewer —
  plus `cluster_labels_<key>` once you name any clusters. Ticking the clustering now takes
  all three, listed in the confirmation, instead of leaving an identical copy behind. The
  bare twin is only ever paired when its `clustering_<name>` exists and the name is not
  one of the Xenium table's own columns, so no raw column becomes selectable.

## [Unreleased] — 2026-07-29

### Fixed (found by replaying a real session)
- **The exported notebook could not get past the dotplot.** `plot:dotplot:<key>`
  recorded `sc.pl.rank_genes_groups_dotplot(adata, …)`, but the rank-genes step writes
  its result to `adata_norm`; the replay died at cell 29 of 39 with
  `KeyError: 'rank_genes_groups'`. `plot:rank_panel:<key>` had the same bug, and
  `plot:volcano:<key>` had its silent form — `run_pairwise_deg(adata, …)` runs happily
  against **raw counts** and returns different genes. All three were written when the
  viewer normalised `adata` in place and were never updated when `normalize` moved to
  binding `adata_norm` — exactly the drift `run_step` prevents for migrated steps.
  Guarded twice now: statically, by a check that no recorded cell passes `adata` to a
  rank-genes consumer, and by executing the recorded dotplot string against a ranked
  `adata_norm`.

- **The notebook kept only the last ranking.** Each rank-genes cell rebinds `rank_df`,
  and scanpy overwrites `uns['rank_genes_groups']` in place, so a session that ranked
  two clusterings exported markers for one of them. `_RANK_GENES_TEMPLATE` now also
  writes `rank_results[<groupby>]`, and the verification dumps one tagged frame per
  clustering.

- **The verification's own comparison was wrong in two ways**, both of which cost a
  ten-minute replay to discover:
  - It compared the viewer's stored ranking against whichever `rank_df` the notebook
    bound last. On a session with two clusterings that meant igraph's 31 groups against
    leidenalg's 28 — reported as every group "diverged" when nothing had: right genes,
    group numbering from a different clustering. `compare_rank_genes` now selects the
    replayed ranking by name and reports `different_groupby` when the notebook never
    ranked what the viewer stored.
  - It filtered unlabelled cells by comparing the *rendered* label to `"nan"`. Under
    pandas 3 a null in a categorical renders `<NA>`, so it passed through into sklearn
    (`ValueError: Input contains NaN`) — hit on a CNV clustering computed on a subset
    and reindexed onto the whole table. Now masked on `.notna()`, with the labelled and
    unlabelled counts in the report so partial coverage reads as what it is.

  Also: `--graph PATH` replays a given `prov_graph.json` against a dataset, which is
  what makes a corrected recording measurable before the user happens to repeat the
  action that recorded it.

  Result on the reference dataset (19 nodes, 63,355 cells, 327.8 s): ARI **1.0** with
  identical labels on all three clusterings — `leiden_igraph_r1.0` (31), 
  `leiden_leidenalg_r1.0` (28), `cnv_leiden_res0.2` (27 over 12,157 labelled cells) —
  and top-10 ranked genes identical in all 31 groups.

### Added
- **The notebook now records what it was run with.** A replay only reproduces a result
  against the same software, and the recorded code named the functions but never the
  versions that answered the call — so a disagreement gave no way to separate a real
  difference from a scanpy upgrade. An **`environment`** node now opens every exported
  notebook: the versions present when the analysis was recorded, as a comment block, plus
  `random.seed(0)` / `np.random.seed(0)` and `sc.logging.print_header()`, so the replay's
  own versions print directly beneath the recorded ones. It is deliberately *not* an
  assertion — a version mismatch is information, not a failure.

  It has no dependents by design: an environment change is something to read, not a
  reason to flag every clustering and DEG table in the session stale. Re-opening a
  dataset in an unchanged environment leaves the node alone rather than rewriting its
  timestamp. The CI replay test executes the cell in a clean kernel, since a version pin
  that raises would be worse than none. (`utils/environment.py`, `tabs/_helpers.py`)

- **Comment-only nodes are now either code or declared as notes.** A recorded node whose
  cell is a comment replays as a silent no-op — `allow_errors=False` sees a cell that ran
  fine, the notebook "passes", and the step it claims to document is simply absent. Some
  of those nodes were real gaps; others were viewer state (the canvas background, an
  overlay) that has no notebook equivalent at all. Both looked identical to every
  consumer, so the Tier-2 report's punch list was mostly display state with the real
  defects buried in it.

  - A new node kind, **`NOTE`**, declares "viewer state, no code equivalent". It renders
    as markdown in the notebook (marked as such), keeps its comment in the flat
    `analysis.py`, is labelled in the Notebook tab, and `verify_notebook.py` counts it
    separately from the punch list. The canvas background, the cluster size filter, the
    UMAP window, patch and transcript overlays, and crop-export are now notes.
  - **ROI expression is real code.** `roi_expression:<gene>` — the one node the first
    Tier-2 run against a real dataset flagged — recorded two lines saying the per-region
    means were "shown in the viewer". It now runs as a `Step`: shapely point-in-polygon
    membership (the same idiom as the ROI DEG step), per-region count/mean/median/std/
    min/max, and pairwise Welch's t-tests with Benjamini-Hochberg correction via
    `scipy.stats.false_discovery_control`. The tab formats its text from that step's
    outputs instead of computing them itself, and `export:roi_expression` is now the
    `to_csv` that writes the file rather than a comment saying one was written.
  - **H&E/ARMS registration nodes carry their data**: the flips bind
    `he_flip_vertical`/`he_flip_horizontal`, the coarse alignment records the affine
    matrix it computed (previously discarded, with only its scale printed), and saving
    landmarks records the `save_landmarks(...)` call with the points inlined.
  - **A source guard** (`tests/test_recorded_code_is_code.py`) parses every
    `ctx.record_node` call site and fails if one records prose where the notebook needs a
    statement. One known gap remains, listed with its reason: `viewer:transcript_density`
    computes a 2-D histogram and needs the transcript loader expressed as plain code
    first.

- **Recorder failures are now visible.** `record_node` degrades rather than aborting when
  provenance bookkeeping fails — a bug there must never lose the user's analysis — but it
  announced the degradation with `warnings.warn`, which in a GUI process goes to a
  terminal nobody reads, and only once per unique message under Python's default filter.
  What was left behind is exactly the failure this work exists to prevent: a result on
  screen with no cell that produces it. `reporting.report_recording_failure` now logs with
  a traceback, keeps the failure for the session tally, and shows a non-modal napari
  warning naming the node — no dialog, since the analysis itself succeeded.

- **Notebook replay verification — the reproducibility claim, measured.** Until now
  nothing executed an exported notebook and compared its results to the viewer's. The
  step executor makes the recorded code *be* the executed code by construction, but
  whether replaying that code from the raw output reproduces the result is an empirical
  question, and it was unanswered. Two artifacts now answer it:

  - **`tests/test_notebook_replay.py`** (CI gate, hermetic, ~35 s). Runs the real
    `Step` templates — Leiden, normalize, rank genes, spatial neighbours, neighbourhood
    enrichment — over a synthetic AnnData (`replay_adata` in `tests/conftest.py`), exports
    the provenance graph as a real `.ipynb`, and executes it in a **clean kernel** with
    `allow_errors=False`. Adjusted Rand index must be exactly **1.0** and the labels
    identical (ARI alone is blind to relabelling); top-N ranked gene names must match in
    order; nhood z-scores must be `allclose`. Further tests assert that the notebook's
    cells are the recorded node sources *verbatim*, so a passing replay cannot be the
    exporter quietly fixing something up. One documented substitution: the `preamble`
    node reads an h5ad instead of `spatialdata_io.xenium(data_path)`, which CI has no
    dataset for — the same preamble exception already documented, and a test asserts it
    stays the only one.
  - **`scripts/verify_notebook.py`** (evidence, run against a real dataset). Reads the
    provenance graph straight out of `<cache>/viewer_session` attrs — no GUI, no napari —
    replays it against the raw Xenium output with per-cell timing, and emits a JSON
    report: per-clustering ARI and cluster counts, top-N gene agreement, wall-clock,
    package versions, and **the ids of every comment-only node the notebook silently
    skipped**. Those nodes execute fine and do nothing, so no amount of `allow_errors`
    catches them; naming them turns the remaining recording work into a measurement.
    `--dry-run` produces that list in seconds without replaying. A replay that *fails*
    reports the **node id** that broke, not just nbclient's cell index, along with the
    timings of every cell that did run.

  Supporting: `notebook_export.execute_notebook()` runs a notebook in a throwaway
  kernelspec pointing at `sys.executable`, because the installed `python3` kernelspec
  belongs to whichever environment registered it last — routinely the conda base env,
  which has no scanpy. `nbclient`/`ipykernel` added to `environment.yml` and to a new
  `test` extra; ruff's CI gate now covers `scripts/` too.

- **Choosable Leiden flavour in the Clustering tab.** `sc.tl.leiden` has two backends —
  `igraph` (fast) and `leidenalg` (scanpy's historical default, optimising the
  RBConfiguration objective rather than igraph's modularity). The viewer hard-coded
  `flavor='igraph'`, so a partition from an existing scanpy pipeline could not be
  reproduced. A **flavor** dropdown now selects between them, `igraph` remaining the
  default, alongside an **n_iterations** spinbox that resets to the selected backend's
  default (`2` for igraph, `-1` — iterate to convergence — for leidenalg). `directed` is
  derived from the flavour rather than exposed, because scanpy raises on
  `directed=True` under igraph. All four values are written literally into the recorded
  step, so the notebook shows exactly which backend produced the labels; the replay test
  now runs over both flavours, since each is seeded from `random_state`.
  (`tabs/tab_clustering.py`, `tests/test_clustering_step.py`)

  **Result keys now carry the flavour**: `leiden_igraph_r1.0` / `leiden_leidenalg_r1.0`,
  so both backends can be run at one resolution and compared instead of one silently
  overwriting the other. Clusterings computed before this change keep their older
  `leiden_r{resolution}` keys and still load and appear in every dropdown, but a new run
  at that resolution writes the new key *alongside* the old one rather than replacing it.

### Fixed
- **inferCNV failed under pandas 3 with `ArrowInvalid: only handle 1-dimensional
  arrays`.** infercnvpy's `_running_mean` slices the gene list with a **2-D** index
  array. Under pandas 3 `var.index.values` is an `ArrowStringArray`, which routes that
  to pyarrow's `take()` — it accepts only 1-D indices. The recorded template already
  carried a shim converting `.obs`/`.var` strings to object dtype, and it did nothing:
  **AnnData re-infers string dtypes when a frame is assigned back**, so with
  `future.infer_string` at its pandas-3 default the Arrow array landed straight back
  where it started. The option has to be off across the *assignment*, not just the
  conversion.

  This is precisely the drift `utils/steps.py` exists to prevent, one level down: the
  in-process helper `_convert_adata_arrow_strings` had the option toggle, the template's
  hand-inlined copy of it did not — so the CopyKAT path (which calls the helper) worked
  while inferCNV (which runs the template) died. `tests/test_cnv_step.py` now *executes*
  the real shim and asserts the 2-D indexing it exists to enable, restores the global
  option it changes, and pins it against the helper.

- **Saving the CNV heatmap crashed when the run had few windows.**
  `chromosome_heatmap(dendrogram=True)` ends in `sc.tl.dendrogram`, which represents
  cells with `pd.DataFrame(_choose_representation(...))`. Above `settings.N_PCS` (50)
  columns that representation is a PCA — dense, fine. At or below it, it is `.X` itself,
  and `pd.DataFrame(csr_matrix)` does not densify: it builds a one-column *object* frame
  of 1×n row matrices, so the `.groupby().mean()` that follows dies with
  `TypeError: agg function failed [how->mean,dtype->object]`. The heatmap therefore
  worked on a wide CNV matrix and crashed on a narrow one. `make_cnv_heatmap` now
  densifies a narrow sparse `X_cnv` for the duration of the plot and restores the
  original afterwards, so the live session object is unchanged.

- **A CNV run on a non-human panel produced a result instead of an error.** InSituCNV's
  default gene-position reference is the infercnvpy Maynard 2020 table, which is human.
  A mouse panel matches only the symbols spelled identically in both nomenclatures — **8
  of 5006** on the dataset that surfaced this (`C2`, `C3`, `C6`, `C7`, `F3`, `F8`, `F9`,
  `H19`). The pipeline ran, clustered those 8 genes into 5 windows and reported CNV
  clusters; the first sign of trouble came several steps later, as the heatmap crash
  above. `run_cnv_pipeline` now refuses to continue when under 5% of the panel has
  coordinates, and says so — naming the counts, and adding that the symbols look like
  mouse nomenclature when the casing suggests it. Supplying a non-human annotation is
  not yet wired up (`prepare_cnv_input` accepts `gene_reference`/`gene_reference_path`);
  until it is, CNV inference is human-only.

- **Exported notebooks died on any viewer-derived clustering.** `record_clustering` is
  the backstop that gives a clustering a `clustering:<key>` node so analysis tabs can
  declare `deps=[...]` on it. It recorded `pd.read_csv(".../analysis/clustering/<key>/
  clusters.csv")` **whether or not that file existed** — true only for the clusterings
  10x ships, false for every one the viewer derives. The first replay of a real session
  against its own dataset died there with `FileNotFoundError`, three cells in.

  The producers now record the code that actually made the column: the CNV tab publishes
  its `cnv_leiden_res*` labels from `cnv_clusters` (inferCNV) or the subsampled
  `adata_copykat` (CopyKAT), each propagated CopyKAT column gets its own
  `clustering:<col>_propagated` node carrying the overlay that used to be a loop hidden
  inside the extrapolation cell, and Novae records `clustering:novae_domains` — under its
  old `novae` id nothing could depend on it, and the recorded cell bound
  `novae_domain` while the viewer stored `novae_domains`. A source guard fails if a new
  producer persists a clustering without recording a node for it.

  The backstop itself, now reached only for columns from a session recorded before its
  producer recorded code, emits a *reload* from the viewer's cache and says so in-line —
  the CopyKAT precedent for code that cannot be the code that ran. It is tested by
  executing it, not by reading it.

- **The provenance graph reached disk only at exit.** Artifacts were persisted eagerly
  (`save_clustering_to_adata` writes the column immediately) but the graph explaining
  them only in `save_session`, which runs on a dataset switch or viewer exit. Measured on
  a real session: the store held a 16-minute-old three-node graph while the same table
  already carried two Leiden clusterings, two rank-genes results, an ROI DEG and a
  neighbourhood enrichment. Anything reading the store mid-session — the verification
  script, the next launch after a crash — saw results with no code behind them.

  The graph is now written to `viewer_cache/prov_graph.json` on every recorded step (one
  small atomic write; updating the zarr group would mean copying every parquet under
  `viewer_session/`). `save_session` still writes the attr, and the sidecar takes
  precedence on load and in `scripts/verify_notebook.py`, which reports which source it
  used. Datasets without the sidecar restore from the attr exactly as before.

  Three guards, because the first version of this ate a graph twice: writes are gated on
  `state["prov_graph_restored"]`, which `app.py` sets *after* the session is restored —
  tabs seed a preamble node while the viewer is still being built, and persisting that
  replaced a 13-node graph with a one-node stub which the next launch then preferred, so
  the DAG came up showing only "Setup & data loading". And on load, a sidecar with
  *fewer* nodes than the attr loses: it is written on every step and the attr only at
  exit, so it is never legitimately smaller, and nothing in the GUI removes nodes.
  Finally `save_session` refuses to shrink the stored graph at all — a viewer that came
  up empty used to write that emptiness over the last remaining copy on exit, which is
  how the 13-node analysis was lost a second time.

- **A write-failure dialog could block a process with nobody at the keyboard.**
  `reporting._surface` guarded the modal only on `QApplication.instance() is None`. A
  test run, a headless script or CI creates an instance with no event loop and no user,
  and `QMessageBox.exec_()` then blocks forever with nothing able to dismiss it — it hung
  the test suite for an hour, silently, as soon as a new fixture created the
  QApplication before the tests that deliberately inject write failures. The modal is now
  suppressed when `QT_QPA_PLATFORM` is `offscreen`/`minimal`/`vnc`; the log entry, the
  failure tally and the non-modal notification are unaffected.

- **Crash-safe zarr cache writes.** The viewer persisted elements with
  `delete_element_from_disk` followed by `write_element`. That is not a metadata
  operation — spatialdata does `del root[element_type][element_name]`, which recursively
  unlinks, so the bytes were gone before the replacement started being written. Its own
  docstring warns "data loss may occur if the execution is interrupted during writing."
  `_persist_table` ran it on *every* analysis action, so every clustering, DEG run and
  label edit opened a window in which a kill, a full disk or any exception left the store
  structurally invalid — and the loader then discarded the whole cache (30 GB on the
  dataset that prompted this, of which the table is 320 MB).

  `utils/zarr_safe.py` replaces it with stage-then-swap: the new element is written to a
  throwaway sibling store under `.xv_staging` (live store untouched, so a failure there
  costs nothing), then journalled and swapped in with two `os.rename` calls. The previous
  copy moves to `.xv_trash` rather than being deleted. `recover_pending()` finishes or
  unwinds an interrupted swap at startup, inferring the phase from the filesystem.
  All 16 delete-then-write call sites now use it, and the store lock covers every writer
  rather than 3 of ~20. (`utils/zarr_safe.py`, `utils/adata_persistence.py`,
  `tabs/tab_he_registration.py`, `tabs/tab_arms.py`, `tabs/tab_external_images.py`,
  `tabs/tab_patch_overlays.py`)

- **Four save functions erased data when a layer was transiently empty.**
  `save_rois_to_sdata`, `save_annotations_to_sdata`, `save_landmarks_to_sdata` and
  `save_arms_tiles_to_sdata` deleted the stored element *before* checking whether there
  was anything to write. A napari layer that was empty for any reason — mid-teardown, a
  missed snapshot — silently wiped the persisted ROIs, annotations or tiles.

- **Session save destroyed the session it was writing.** `save_session` called
  `create_group("viewer_session", overwrite=True)` and only wrote the replacement ~110
  lines later; any exception in between (most plausibly a non-serialisable value reaching
  `attrs.update`) left an empty group and one printed warning, after the user had closed
  the window. Attrs are now built and JSON-validated before the group is touched, and the
  write goes through `safe_group_update`. (`utils/session.py`)

- **Startup migrations re-armed on every launch.** `save_session` rebuilt attrs from
  scratch, wiping the four `migrated_*` markers each clean exit — so migrations re-ran at
  the next launch, including two that themselves rewrote the whole cell table. Attrs are
  now preserve-by-default. (`utils/session.py`)

- **Caches were discarded too readily.** Three paths did it:
  - an unreadable cache was renamed aside — or `rmtree`'d if the rename failed — with no
    repair attempt and no check for user data. It now runs `verify`/`repair` first and
    escalates to restoring a missing element from its backup.
  - staleness compared `experiment.xenium`'s mtime against the cache *directory* mtime,
    which only moves when a direct child is added or removed — so `rsync`/`cp -p`/a
    re-download condemned a good cache. Caches now carry `.xv_manifest.json` with a
    sha256 of the source. Pre-existing caches keep the mtime check as an *uncertain*
    hint that prompts rather than rebuilding.
  - the sidecar list omitted `adata_cnv_cache_*.h5ad`, so a cache whose only user data
    was a multi-hour CopyKAT run reported "no user data" and was rebuilt with **no dialog
    at all**.

  Also: `_restore_user_elements` is now a deny-list (CNV scores, `cnv_runs` and
  `rank_genes_groupby` were dropped even on "Rebuild and restore my data"); the rebuild
  stages and renames rather than overwriting in place and `rmtree`ing on failure, after a
  free-space check; and no-dialog no longer means "rebuild" — stale keeps, unopenable
  raises `CacheLoadAborted` for a clean exit. (`loader.py`, `app.py`)

- **The Marker Genes correlation-matrix button never worked** — it called
  `sc.tl.correlation_matrix`, which does not exist in scanpy. See the 2026-07-28 entry.

### Added
- **Recovery from a corrupt cache no longer needs the corrupt cache to open.** The first
  version of "Recover from Backup" called `spatialdata.read_zarr` on the backup — which is
  by definition a store that failed to read — and died in `_read_table` on an unreadable
  table, taking the perfectly salvageable shapes, images and clusterings with it. Recovery
  is now filesystem-level: element directories are self-contained and obs columns are
  individual zarr arrays, so a broken root index or an unreadable table condemns neither.
  `cache_repair.salvageable_elements()` and `read_obs_columns()` do the reading;
  `zarr_safe.safe_import_element()` does the writing, journalled like any other swap.
  Verified against the reported cache: 12 obs columns (8 clusterings, cluster labels, 3 CNV
  scores), 6 landmark sets and 3 images recovered from a store spatialdata refuses to open.
  (`utils/cache_repair.py`, `utils/zarr_safe.py`, `tabs/tab_cache.py`)

- **External images and patch overlays were not counted as user data.** `_detect_user_data`
  matched a fixed list of element names, but these are named per file
  (`ext_<filename>`, `ext_<filename>_xenium_lm`, `patch_*`) — so a dataset with a registered
  PhenoCycler image and its landmarks reported "no user data" and could be rebuilt over
  without a prompt. Matching is now by prefix/suffix as well as exact name. (`loader.py`)

- **`errored` handlers indexed the exception as an exc_info triple.** napari emits the
  exception itself, so `exc_info[1]` raised `TypeError` and replaced the real error in the
  traceback — masking, among other things, the recovery failure above. Fixed in the Cache
  tab and at the three pre-existing sites in `tabs/tab_segmentation.py`, where it had been
  hiding async segmentation-save failures. (`tabs/tab_cache.py`, `tabs/tab_segmentation.py`)

- **Recovered registration was undone by the reload that followed it.** Reloading saves the
  current session first, as any dataset switch does, and `save_session` deletes the
  `he`/`arms` groups and rewrites them from `ctx.he_state` — which was still empty, so it
  erased the affine recovery had just written. The images came back but unaligned.
  Recovery now hydrates `ctx.he_state` / `ctx.arms_state` in memory as well as on disk, so
  the pre-reload save writes the recovered values instead of blanking them. Landmarks were
  never affected — they load from `sdata.shapes`, so importing those elements was already
  enough. (`tabs/tab_cache.py`)

- **External-image and patch-overlay UI state was blanked on every save with none
  loaded.** `_snapshot_layers` yields `[]` rather than `None` when nothing is loaded, so
  the "fall back to the previous value" branch never fired — losing saved contrast,
  opacity and affine after a recovery, and any time the attrs were written before restore
  had run. These now fall back on empty too, which is safe because restore is driven by
  the sdata elements with the attrs used only as decoration: an entry left behind for a
  removed image is never looked up. (`utils/session.py`)

- **Recovered data was invisible until the dataset was reopened.** Recovery writes
  elements straight into the zarr, so the live SpatialData, the napari layers and every
  tab's widgets knew nothing about them. `app.py`'s dataset-open path is now factored into
  `_load_dataset_into_viewer(path)` and exposed as `ctx.reload_dataset`, so the Cache tab
  offers a reload as soon as recovery finishes. Recovery also merges the backup's
  `viewer_session` — H&E/ARMS filenames, flips and affines — for keys the live session
  lacks: restoring `he_image` and its landmarks without that is half a job, since the
  element would exist while the session still said no H&E was loaded.
  (`app.py`, `tabs/tab_cache.py`, `utils/viewer_context.py`)

- **A Cache tab (Tools → Cache).** The cache was a black box: when it broke, the loader
  moved it aside and rebuilt from raw, and the only signal was a long wait. The tab shows
  size, free space, build manifest, write failures this session and the log path, and
  offers **Verify** (read-only), **Re-consolidate Metadata** (fixes the most common
  corruption without touching data), **Recover from Backup** (pull elements out of a
  `.xv_trash` copy or a previous cache the loader kept aside, including CopyKAT sidecars)
  and **Force Rebuild** (moves the cache aside and rebuilds on the next launch, so the
  running session is never left pointing at a freed store). All work runs in a
  `thread_worker` behind `store_lock`, so nothing here can race `_persist_table`, and a
  test fails if the tab ever gains an `rmtree`. (`tabs/tab_cache.py`, `app.py`)

- **`utils/reporting.py` — a per-dataset log and non-modal error surfacing.** Every write
  failure went to stdout, which a GUI user never reads, and the only dialog was for
  permission errors, shown **once per process** via a module-level flag — so the second
  failure and everything after it was invisible. When the cache was being corrupted, the
  warnings that would have explained it were lost. Now: a rotating
  `<data_path>/xenium_viewer.log` (2 MB × 3) started before anything can fail;
  `report_write_failure()` always logs with a traceback, shows a non-modal napari
  notification marshalled to the GUI thread via `ensure_main_thread` (the old dialog could
  be constructed from a `thread_worker`, a real Qt violation), and reserves a modal for
  permission and disk-full errors only — tracked per (dataset, error class), not per
  process. A running tally (`failure_summary()`) makes failures visible in aggregate
  without a popup per event. The ~20 `print("Warning: could not ...")` sites in the cache
  write paths now log, so the file captures them; a source guard fails if any come back.
  (`utils/reporting.py`, `utils/adata_persistence.py`, `utils/session.py`, `app.py`)

- **`utils/cache_repair.py`** — `verify()` (read-only; parses the root `zarr.json` with
  `json.loads` rather than `zarr.open`, so it reports on a store too broken to open) and
  `repair()` (idempotent; replays journals, clears debris, drops stray groups,
  re-consolidates, and at `FULL` restores an element from its `.xv_trash` backup).
  Replaces the ad-hoc block in `app.py`, which handled two hard-coded cases and assumed a
  nested consolidated-metadata layout — zarr 3.1 writes a flat one, so it could not have
  detected the case it was written for. (`utils/cache_repair.py`, `app.py`)

### Changed
- **Sidecar analysis outputs moved out of the zarr store root** into
  `<data_path>/viewer_cache/`. Files in the store root make zarr's hierarchy walk emit a
  `ZarrUserWarning` each — the most likely source of the reported "several warnings",
  since `app.py` called `consolidate_metadata` without the filter spatialdata itself uses.
  It also meant a cache rebuild deleted them, including hours of CopyKAT compute. Readers
  fall back to the old location, so existing datasets keep working and nothing is migrated
  eagerly. (`utils/adata_persistence.py`, `tabs/tab_cnv.py`, `loader.py`)

### Tests
- First coverage the zarr/persistence paths have ever had: `test_zarr_safe.py` (26,
  including interrupted-write simulation with both a recoverable `OSError` and a
  `KeyboardInterrupt` that bypasses cleanup the way a `kill -9` does),
  `test_persistence_safety.py` (9), `test_session_persistence.py` (14),
  `test_cache_repair.py` (20), `test_loader_policy.py` (16), `test_sidecar_location.py`
  (20), `test_reporting.py` (21), `test_tab_cache.py` (24). Plus source guards that fail if `delete_element_from_disk` is called outside
  `zarr_safe.py`, if `loader.py` `rmtree`s the live cache, or if a sidecar is written into
  the store root, or if a cache write path prints a warning instead of logging it, or if recovery opens a backup as a whole.
  289 tests pass.

## [Unreleased] — 2026-07-28

### Added
- **Step executor: the code the GUI runs is now literally the code the notebook
  records.** New `utils/steps.py` introduces `Step` (a provenance node id, a
  `string.Template` of plain scverse source, and a dict of literal `params`) plus
  `StepExecutor`, which renders the template **once** and hands that same string both
  to `exec` and to `ProvGraph.upsert`. This is the E1 infrastructure for closing the
  drift between executed and recorded code — the defect that let the GUI run
  `leidenalg.find_partition(..., seed=42, n_iterations=2)` while the notebook recorded
  `sc.tl.leiden(..., random_state=0)` (scanpy's default is `n_iterations=-1`), and let
  the GUI normalise with `target_sum=1e4` while the notebook recorded scanpy's median
  default. With one rendering there is no second expression to drift.

  The invariant that makes the guarantee hold, and which review must enforce: *a tab
  callback may never call an analysis function with a widget value — it may only build
  a `params` dict.*

  Design notes:
  - `string.Template` (`$name`) rather than `str.format`, so `{...}` dict literals and
    f-strings inside templates are left alone.
  - Params are substituted via `repr()` and validated with
    `ast.literal_eval(repr(v)) == v`, which rejects numpy scalars (whose NumPy-2 repr
    is `np.float64(1.0)` — not importable in a bare notebook, and not stable across
    versions), non-finite floats, and objects with a default `<... at 0x...>` repr.
    `coerce()` is provided for use at the widget boundary. Float noise such as
    `1.0000000000000002` round-trips exactly and is therefore allowed.
  - Execution is serialised behind an `RLock` (steps mutate shared namespace state) and
    proceeds one top-level statement at a time via `ast`, so long steps can report
    progress while the compiled source stays byte-identical to the recorded source;
    statement line numbers are preserved so tracebacks point into the recorded cell.
  - `compile(..., "<step:id>")` puts the step id in the traceback, and failures raise
    `StepError` naming the step and statement instead of being swallowed.
  - Recording happens only on success, so a failed step leaves no node claiming an
    artifact that does not exist.
  - `free_names()` / `check_step()` provide the template lint: a rendered template must
    reference only names the namespace guarantees. This is what makes the exactness
    guarantee auditable rather than merely asserted.
  - `params` finally becomes meaningful — every node now carries a machine-readable
    parameter record alongside its source (it was previously populated at exactly one
    call site and never rendered).

  `tests/test_steps.py` (36 tests) covers the executed-equals-recorded guarantee, param
  round-tripping for every literal type the GUI can produce, numpy rejection plus
  coercion, brace survival, free-variable analysis, per-statement progress, failure
  naming, and the inherited upsert/staleness semantics on re-run with changed params.
  No existing module imports `steps.py` yet — tabs migrate in E2.
  (`utils/steps.py`, `tests/test_steps.py`)

- **Leiden clustering migrated onto the step executor — the first analysis whose
  recorded code is the code that ran.** `tab_clustering.py` now builds a single `Step`
  from `_leiden_template(use_hvg, do_scale)` and calls `ctx.run_step()`; the separate
  `_record_leiden_code` branch that hand-wrote a parallel description of the pipeline is
  gone. `ctx.run_step` / `ctx.executor` are attached in `create_shared_helpers`, with a
  namespace seeded with `sc`/`sq`/`pd`/`np`/`plt`/`Path`/`data_path`/`sdata`/`adata` so
  steps operate on the same objects the viewer holds.

  Behaviour changes that follow, all of them fixes:
  - **Leiden now runs `sc.tl.leiden(flavor='igraph', n_iterations=2, random_state=0)`**
    instead of `leidenalg.find_partition(..., seed=42)` on a directed graph in a spawned
    subprocess. `flavor`, `n_iterations` and `random_state` are pinned explicitly because
    scanpy 1.12 emits a `FutureWarning` that the default backend will become `igraph`
    (which also requires `directed=False`) — leaving them implicit would silently change
    clusterings on a scanpy upgrade. **Cluster labels will differ from previous runs**;
    existing saved clusterings are untouched.
  - **Both preprocessing branches are now one step**, which starts from `adata_norm`
    (bound by the `normalize` step below) and works on its own `adata_leiden` copy, so
    neighbours/Leiden/HVG-subsetting don't mutate the shared normalised object. It
    declares `deps=["normalize"]`, so the DAG carries a real `normalize -> clustering`
    edge and the exported notebook normalises exactly once. PCA is recomputed only when
    HVG selection or scaling changed the feature space; otherwise `adata_norm`'s `X_pca`
    is what a recomputation would produce anyway.
  - Flat-journal/notebook-tab updates are bounced to the GUI thread via
    `superqt.utils.ensure_main_thread`, since steps execute in napari worker threads.

  `tests/test_clustering_step.py` (12 tests) runs the real `normalize` + Leiden pair on
  synthetic data across all four HVG/scale combinations, asserts recorded source ==
  executed source, and — the reproducibility claim in miniature — replays the whole
  topo-sorted graph in a clean namespace and checks the labels match exactly.
  (`tabs/tab_clustering.py`, `tabs/_helpers.py`, `utils/viewer_context.py`,
  `tests/test_clustering_step.py`, `CLAUDE.md`)

- **`normalize` and rank genes migrated onto the step executor — the second known
  divergence closed.** The viewer normalised with `target_sum=1e4`
  (`gene_analysis.get_normalized_adata`) while the recorded cell said
  `sc.pp.normalize_total(adata)`, i.e. scanpy's *median* default. Different `X`, so
  different PCA, neighbours, clusters and DEG. `ctx.ensure_normalized()` replaces the
  old `record_normalize` + `get_normalized_adata` pair with a single step that records
  the `target_sum` it uses.

  Two structural fixes come with it:
  - **`normalize` binds `adata_norm = adata.copy()` instead of mutating `adata`.** The
    old node was `kind=SETUP`, so it sorted ahead of every artifact node and silently
    log-normalised the object other cells then copied — an implicit, invisible edge that
    was wrong whenever the node happened to be absent. Consumers now name `adata_norm`
    explicitly and declare `deps=["normalize"]`, which is what puts the edge in the DAG.
  - **Rank genes now ranks on `adata_norm`, not on raw `adata`.** The recorded cell had
    been `sc.tl.rank_genes_groups(adata, ...)`, correct only when a `normalize` SETUP node
    happened to be in the graph and had mutated `adata` first. The step declares
    `deps=["normalize", "clustering:<key>"]` and copies the clustering from `adata.obs`
    onto `adata_norm.obs`.

  `record_clustering`'s CSV-import node now depends on `preamble` rather than `normalize` —
  it reads a CSV into `.obs` and never needed normalised values. `ctx.record_normalize` is
  gone from `ViewerContext`, replaced by `ctx.ensure_normalized`.

  `tests/test_normalize_rank_genes_steps.py` (9 tests) asserts the step reproduces
  `get_normalized_adata`'s output to `rtol=1e-6`, that `adata` is left untouched, that
  rank-genes reads `adata_norm`, and that replaying the whole topo-sorted graph in a clean
  namespace yields an identical `rank_df`.
  (`tabs/_helpers.py`, `tabs/tab_gene_analysis.py`, `utils/viewer_context.py`,
  `tests/test_normalize_rank_genes_steps.py`, `CLAUDE.md`)

- **The spatial, marker, correlation and ROI tabs migrated onto the step executor.**
  Every one of them ran one expression and recorded a different one; each is now a single
  templated `Step`.

  - **Spatial neighbours** is now `ctx.ensure_spatial_neighbors(k)`, building the graph on
    `adata_norm`. The old node built it on `adata` while every consumer was handed the
    normalised copy — so a replayed notebook ran neighbourhood enrichment against an object
    with no `.obsp` graph on it at all. It also stopped rebuilding `obsm['spatial']` from
    `x_centroid`/`y_centroid`, columns the Xenium table does not have under those names.
  - **Neighbourhood enrichment** and **co-occurrence** run on `adata_norm` with that graph.
  - **Ligand-receptor**: the interaction-database checkboxes reached the notebook only as a
    `# interactions: OmniPath, LigRecExtra` prose comment, so a replay silently fell back to
    omnipath's defaults. They are now `InteractionDataset` members reconstructed by name in
    the recorded source, along with the `CellPhoneDB`-only restriction, `use_raw=False` and
    `copy=True`.
  - **Marker genes** recorded *nothing at all* despite being five plain scanpy calls. All
    five are now `plot:markers:<plot>:<key>` terminals carrying the marker dict and the
    cluster display labels as literals.
  - **Gene correlation**: the whole figure is one step, so the scatter the viewer shows is
    the scatter the notebook draws. The recorded cell previously omitted the annotation box,
    the *n* in the title, and the cluster filter entirely — a notebook that correlated all
    cells while the GUI showed a filtered subset. Expression is pulled with `sc.get.obs_df`
    instead of hand-indexing `.X` and calling `.toarray()`.
  - **ROI DEG** no longer records a call into `xenium_viewer.utils.gene_analysis`, so the
    notebook is standalone scverse code: the shapely point-in-polygon assignment and the
    scanpy DEG are written out in full. `rois` became a real step that binds `roi_polygons`.
    The cluster filter is recorded as an explicit `obs[key].isin([...])`.

  `tests/test_spatial_roi_steps.py` (28 tests) executes the real templates and replays the
  recorded graph in a clean namespace; ROI DEG is additionally checked frame-for-frame
  against `compute_roi_deg`, the implementation it replaces.
  (`tabs/_helpers.py`, `tabs/tab_nhood.py`, `tabs/tab_co_occurrence.py`, `tabs/tab_ligrec.py`,
  `tabs/tab_marker_genes.py`, `tabs/tab_gene_correlation.py`, `tabs/tab_roi.py`,
  `utils/viewer_context.py`, `utils/spatial_analysis.py`, `tests/test_spatial_roi_steps.py`)

- **The inferCNV backend migrated onto the step executor.** `cnv:infercnv` is now a
  templated step running in-process; the tab builds the result dict from its outputs
  instead of calling `run_cnv_pipeline`. Three drifts closed: the recorded cell normalised
  at scanpy's *median* default while the viewer used `target_sum=1e4`, and it dropped both
  `lfc_clip` (so a replay used infercnvpy's default, not the pipeline's 4.0) and
  `dendrogram=False`. The pandas-3 PyArrow-string conversion `infercnvpy` needs is inlined
  as plain pandas, so the notebook no longer imports `xenium_viewer` at all.

  **CopyKAT is deliberately not migrated.** It runs detached, in a second conda env
  (its R stack needs python 3.11), so no in-process step can be the code that ran. Its
  node stays on `ctx.record_node` and now says so in the cell itself — it is a
  reconstruction of that run, not executed source. Its recorded parameters were corrected
  the same way (`target_sum=1e4`, `dendrogram=False`).
  `run_cnv_pipeline` is now the CopyKAT path only; `tests/test_cnv_step.py` (11 tests)
  pins the parameters both sides carry so they cannot drift apart again.
  (`tabs/tab_cnv.py`, `utils/cnv_analysis.py`, `tests/test_cnv_step.py`)

  **Known remaining divergence:** the annotation-neighbourhood tab records nothing at all —
  its synthetic virtual cells are sampled from a napari shapes layer the notebook has no
  access to, which needs E3's spatialdata shapes to resolve. Plot/export **terminals**
  across the migrated tabs are still on `ctx.record_node`; their code strings were
  corrected to read the object the artifact now lives on, but the terminal-node policy
  itself is E4.

### Fixed
- **Re-running a clustering under an existing key kept the previous run's colors.**
  `CellColorManager.get_cluster_colors` caches on the series' `name`, which is the
  clustering key. Re-running Leiden at the same resolution (or re-importing the same
  file) replaces the series behind that key, but nothing dropped the cache entry — so
  the raster kept the old color array while the legend and cluster filter were rebuilt
  from the new assignment. On screen that reads as a clustering that was only *partially*
  overwritten: cells whose new cluster id happened to match their old one looked right,
  the rest did not. Compounding it, the Leiden tab named every series `"leiden"`
  regardless of resolution, so one cache entry served every run at every resolution.
  Producers now call `color_manager.invalidate_cluster_cache(key)` — a named method
  replacing the `_cluster_cache.pop()` the CNV and Novae tabs were already doing by
  hand — and the Leiden series is named for the key it is stored under.
  `tests/test_cluster_color_cache.py` (7 tests) covers the behaviour and adds a
  source-level guard that fails if any tab rebinds `ctx.clusterings[key]` without
  invalidating. (`utils/coloring.py`, `tabs/tab_clustering.py`, `tabs/tab_cnv.py`,
  `tabs/tab_novae.py`)
- **The Marker Genes correlation-matrix button never worked.** It called
  `sc.tl.correlation_matrix`, which does not exist in scanpy — the `AttributeError` was
  raised inside a worker thread where nothing surfaced it. `sc.pl.correlation_matrix` reads
  the matrix `sc.tl.dendrogram` computes, so that is what the step now runs. Surfaced by
  migrating the tab onto the executor, which reports step failures instead of swallowing
  them. (`tabs/tab_marker_genes.py`)
- **A dead `_get_adata_norm` in the Marker Genes tab** passed `add_clustering_to_obs` its
  arguments in the wrong order (`adata_orig` and `clustering_series` swapped). Removed.
  (`tabs/tab_marker_genes.py`)

### Removed
- **`spatial_analysis.run_ligrec` and `spatial_analysis.run_co_occurrence`.** One squidpy
  call each; the templates in the tabs now *are* that call. (`utils/spatial_analysis.py`)
- **`utils/leiden_worker.py`.** The spawned-subprocess Leiden existed for GUI
  responsiveness, but it was also the second expression of the algorithm that drifted from
  the recorded one. scanpy's `igraph` flavor is, per its own warning, "orders of magnitude
  faster" than `leidenalg`, which removes the motivation. (`utils/leiden_worker.py`,
  `cnv_copykat_worker.py` docstring reference)

### Changed
- **`.gitignore`**: added `manuscript/` (preprint drafts and planning notes, kept out of
  the public repo) and `data/` (untracked local datasets).

## [Unreleased] — 2026-07-19

### Fixed
- **`conda env create -f environment.yml` failed to solve because of `insitucnv`.**
  The `insituCNV-copykat` fork's published metadata carried stale upper pins
  (`anndata<0.12`, `pandas<3`). Since the pip section of `environment.yml` resolves
  `-e .` and `insitucnv` together, pip had to satisfy `anndata>=0.12` (this package's
  own requirement) *and* `anndata<0.12` at once — an unsatisfiable constraint, so the
  environment build died with `ResolutionImpossible`. The bounds were never real: the
  fork only uses stable AnnData APIs and has been exercised end-to-end against
  anndata 0.13 / pandas 3.0. Fixed at the source by dropping both upper bounds in the
  fork (which also drops its unused `seaborn` dependency and replaces the
  matplotlib-3.9-removed `matplotlib.cm.get_cmap` with `matplotlib.pyplot.get_cmap`).
  No dependency line in this repo changed — `environment.yml`, `environment-copykat.yml`,
  and the `cnv` extra all track the fork's master and now resolve unmodified. The
  `pip install --no-deps insitucnv` workaround is retired from the install docs, and
  `_patch_matplotlib_cm_compat()` is now a no-op guard kept only for pre-existing
  environments and upstream InSituCNV. (`README.md`, `docs/Installation.md`,
  `environment-copykat.yml`, `utils/cnv_analysis.py`, `CLAUDE.md`)
- **`mamba env create` prompted for GitHub credentials.** Separately from the pin
  conflict above, the `insituCNV-copykat` fork was a *private* repo installed over
  `https://`, so pip's clone stopped to ask for a username/password that could never
  work (GitHub dropped Git password auth in 2021). The fork is now **public**, so the
  pip section clones anonymously on any machine. This also retires the CI workaround
  that stripped the private dependency out of `environment.yml` — CI now builds from
  the unmodified env file, so the environment it tests matches the one users get.
  (`.github/workflows/ci.yml`)
- **CNV clusterings showed no cells when the "Filter by cluster" checkbox was on.**
  Selecting a CNV clustering (inferCNV `cnv_leiden_res*` or a CopyKAT `*_propagated`
  column) in Cells → Coloring with the cluster filter engaged blanked *every* cell,
  even though the cluster IDs were listed. CNV clusterings carry *string* categories
  (`'0','1','2'`, or `'tumor'/'normal'/'unknown'`) unlike ordinary Leiden's *integer*
  categories, so they take the `factorize` path whose `_cluster_raw_to_id` map is keyed
  by the raw strings — but `_repopulate_cluster_checkboxes` coerced the checkbox ids to
  `int`, so `translate_selected_ids_to_int` matched nothing, returned `[]`, and the
  `~np.isin(...)` mask removed all cells. The checkbox ids now keep their raw category
  type (sorted numerically for display only), matching the map. Reviewed and verified
  the CopyKAT extrapolation/propagation itself is correct — the "sometimes wrong
  propagation" was this same blanking artifact, not a propagation bug. Added
  `tests/test_cluster_filter.py`. (`tabs/_helpers.py`, `tests/test_cluster_filter.py`)

### Docs
- **Tracked TODO to migrate napari off the deprecated PyQt5 backend.** Napari warns
  that PyQt5 support is deprecated and will be removed in fall 2026. No code change
  yet (the viewer runs fine on PyQt5 until then); the codebase already routes all Qt
  access through `qtpy`, so the migration is small — only the backend pins plus a few
  Qt5-isms (8 unscoped enums, 7 `.exec_()` calls). Captured the migration checklist in
  `docs/pyqt6-migration.md` for a future session. (`docs/pyqt6-migration.md`)

### Added
- **Unit tests + continuous integration.** Added a `pytest` suite over the codebase's
  pure logic — `tests/test_prov_graph.py` (extended with Mermaid/DOT rendering),
  `test_cnv_subsample.py` (the CopyKAT budget-split invariant), `test_registration.py`
  (landmark affine + JSON round-trip), `test_gene_analysis.py` (LLM prompt/response
  parsing, label mapping), `test_patch_overlay_io.py` (patch-size/stride inference), and
  `test_notebook_export.py` (graph → `.ipynb` round-trip). A GitHub Actions workflow
  (`.github/workflows/ci.yml`) runs the suite in the full conda env (micromamba) plus a
  fast `ruff` error-only lint gate on every push/PR; README now shows CI / license /
  Python badges. Configured via `[tool.pytest.ini_options]` in `pyproject.toml`
  (`pythonpath = ["src"]`). (`tests/`, `.github/workflows/ci.yml`, `pyproject.toml`,
  `README.md`, `.gitignore`)
- **Global CPU-cores preference wired into CopyKAT.** A new **Preferences → CPU
  cores** submenu sets `ctx.state["n_cores"]` (default `max(1, os.cpu_count()//2)`,
  session-only like the other preferences). The CopyKAT path threads it through
  the launch params → `cnv_copykat_worker` → `run_cnv_pipeline(n_cores=)` →
  `run_copykat(n_cores=)` (CopyKAT's R `n.cores`, which speeds its `parallelDist`
  passes); inferCNV is unaffected. The value is echoed into
  `cnv_copykat_params.json` and the recorded `run_copykat(...)` code cell.
  (`app.py`, `tabs/_helpers.py`, `tabs/tab_cnv.py`, `cnv_copykat_worker.py`,
  `utils/cnv_analysis.py`)
- **Extrapolate CopyKAT calls to the whole dataset.** A new **"Extrapolate CopyKAT
  calls to all cells"** checkbox on the CNV tab. Because CopyKAT runs on a subsample,
  only those cells get a call; when enabled, after the run finishes the viewer
  extends the per-cell tumor/normal (`cnv_status`), `copykat_pred`, and CNV-subclone
  (`copykat_leiden_res*`) results to every cell: cells CopyKAT actually ran keep their
  real value, and each un-run cell is filled with the majority value among run cells in
  its reference-clustering group (the fork's `propagate_cnv_labels(method="cluster")`,
  with run cells' real values overlaid back on top so the true 0/1/2/3 subclones aren't
  collapsed to the dominant one). Adds a colorable, session-persisted `<col>_propagated`
  clustering for each and a `cnv:copykat_propagated` provenance node. These are copied cluster-level calls, not
  per-cell inferred CNV; groups with no sampled cell are labelled `unknown`. Requires
  the updated `insituCNV-copykat` fork (force-reinstalled in both envs).
  (`tabs/tab_cnv.py`, `cnv_copykat_worker.py`)

### Fixed
- **`NameError` on a successful Novae run (surfaced by the new CI lint gate).** In the
  Novae tab, `_on_novae_ready()` referenced `level` when recording the reproducible code
  and building the results summary, but `level` was only bound in `on_run_novae()`'s
  scope (unlike `species`/`n_domains`, which are re-read from their widgets there). Any
  completed Novae domain inference would raise `NameError: name 'level' is not defined`.
  Now `_on_novae_ready()` reads `level = level_slider.value` alongside the others.
  (`tabs/tab_novae.py`)
- **CopyKAT subsample starved the analyzed cells when the reference cluster was large.**
  `subsample_indices` kept *all* reference (baseline) cells first, so a reference
  population bigger than `max_cells` filled the entire subsample and **no analyzed cell
  was ever run** — every analyzed cluster came back with no CNV call and showed as
  `unknown` in the extrapolation (e.g. a 17.9k-cell reference cluster consumed the whole
  10k budget). The budget is now split: analyzed cells get priority for the slots while a
  modest reference baseline (25% of `max_cells`, ≥500) is reserved to seed CopyKAT's
  diploid baseline, and any unused budget tops the reference back up. Small-reference runs
  are unchanged. Existing CopyKAT caches produced before this fix must be re-run.
  (`utils/cnv_analysis.py`)
- **Stale CopyKAT "running" marker after a killed worker.** A SIGTERM/SIGKILL on the
  detached CopyKAT worker (e.g. from htop) skips its `finally` cleanup, leaving
  `plots/copykat_RUNNING.txt` behind so the next launch wrongly reported a job "in
  progress". The worker now records its **PID** in that marker, and the viewer
  distinguishes a genuinely-running detached worker from a dead one by checking
  `/proc/<pid>/cmdline` (with an `os.kill(pid, 0)` fallback; a reused/unrelated PID
  reads as interrupted). Stale markers are now auto-cleared — on restore and in the
  live poll loop — and the CNV tab reports the true state. Terminating the worker was
  already safe for the SpatialData zarr (the worker only writes standalone
  `.h5ad`/`.json`/plot files, never the zarr store). (`tabs/tab_cnv.py`,
  `cnv_copykat_worker.py`)

### Added
- **CopyKAT CNV backend (inferCNV / CopyKAT / both).** The CNV tab can now call
  copy-number with **inferCNV**, **CopyKAT**, or **both** (default), via the
  `insituCNV-copykat` fork's `run_copykat` (which writes the same
  `obsm["X_cnv"]`/`uns["cnv"]` keys, so clustering + heatmaps are shared). A
  per-backend registry keeps both results live at once — cluster columns are
  namespaced (`cnv_leiden_res*` vs `copykat_leiden_res*`), both are colorable,
  and the heatmap saver has a **backend** selector alongside the resolution one
  (`cnv_heatmap_<backend>_<res>.png/.pdf`, fork settings: dendrogram, ±0.4, dpi
  200). CopyKAT is slow (~2 h), so it runs on a random ≤10k-cell subsample
  (reference cells kept) as a **detached background process** that survives the
  GUI: closing the app mid-run prompts **Stop / Continue in background / Cancel**,
  and a finished background run writes `adata_cnv_cache_copykat.h5ad` +
  `cnv_copykat_result.json` + `plots/cnv_heatmap_copykat_*` +
  `plots/copykat_DONE.txt`, restored on next launch. New module
  `cnv_copykat_worker.py`. **CopyKAT runs in a second conda env** — its R stack
  (r-base 4.3 + rpy2 3.5.11) requires python 3.11, which is incompatible with the
  viewer's python 3.12 (scanpy≥1.12 / zarr≥3). Create it once with
  `conda env create -f environment-copykat.yml` (`xenium_viewer_copykat`); the
  viewer auto-detects it (override via `XENIUM_COPYKAT_ENV` /
  `XENIUM_COPYKAT_PYTHON`) and launches the detached worker there, passing the
  viewer source on `PYTHONPATH` so that env needn't install the package. The
  GitHub-only copykat R package auto-installs on first run. inferCNV still runs
  in the main env. (`tabs/tab_cnv.py`, `utils/cnv_analysis.py`,
  `utils/adata_persistence.py`, `app.py`, `install_copykat.py`,
  `cnv_copykat_worker.py`, `environment.yml`, `environment-copykat.yml`,
  `pyproject.toml`)
- **Limit CNV analysis to specific cell types.** A new "Cell types to analyze
  (CNV subclones)" checkbox grid in the CNV tab (drawn from the same
  clustering/annotation column as the reference selector, so it shows your
  cell-type labels) lets you restrict inference to chosen cell types. Before
  running, the AnnData is subset to the selected types **plus** the reference
  population (inferCNV needs the reference as its baseline), so the CNV profile,
  score, subclone clustering, and chromosome heatmap only cover the cells of
  interest — immune/stromal cells you don't care about are excluded entirely.
  Leaving all boxes checked (the default) analyzes the whole tissue as before.
  The analyzed cell-type set joins the CNV profile signature (so it interacts
  correctly with per-resolution heatmap accumulation), persists across sessions
  (`cnv_run_info["analyze_categories"]`), and the `cnv` provenance node emits the
  leading `adata = adata[...].copy()` subset step. (`tabs/tab_cnv.py`,
  `utils/cnv_analysis.py`, `utils/adata_persistence.py`)

## [Unreleased] — 2026-07-17

### Added
- **Per-resolution CNV chromosome heatmaps.** The CNV tab now *remembers* every
  leiden resolution run under the same core parameters (reference, neighbors,
  smoothing, window, step) instead of overwriting the previous run. A new
  **"Heatmap resolution"** selector lets you save the chromosome heatmap for any
  accumulated resolution; each writes its own `plots/cnv_heatmap_<key>.png/.pdf`
  and a distinct `plot:cnv_heatmap:<key>` provenance terminal, so comparing
  resolutions no longer clobbers earlier heatmaps. Accumulation is scoped to a
  shared CNV profile — changing a core parameter starts a fresh profile (with a
  status note), since prior resolutions' heatmaps require the earlier profile.
  The accumulated resolution list persists across sessions
  (`cnv_run_info["cluster_keys"]`) and the `cnv` provenance node emits the full
  `cluster_cnv_resolutions(adata, [...])` list. (`tabs/tab_cnv.py`,
  `utils/adata_persistence.py`)

## [Unreleased] — 2026-07-16

### Added
- **Reproducible analysis as a provenance DAG.** User actions are now recorded
  as a dependency graph of artifacts (`utils/prov_graph.py`) rather than an
  append-only script. Each step is a node keyed by the artifact it produces,
  with its code and dependencies; the notebook is *derived* from the graph by
  topological sort. Re-running a step revises its node in place and flags
  descendants **stale** (instead of being silently dropped or appended out of
  order), and a missing dependency errors at record time rather than as a
  `NameError` on replay. New recorder API `ctx.record_node(id, code, deps, kind,
  label, params)`; `ctx.record_code(code, tag)` kept as a compat shim.
- **Cross-session accumulation.** The graph is serialized into
  `sdata_cached.zarr/viewer_session/` and restored at startup, so a workflow
  spanning several sessions builds one coherent notebook. Output filename is a
  stable `analysis.py` (was a per-launch `code_<timestamp>.py`).
- **Jupyter export.** `utils/notebook_export.py` (nbformat) writes
  `analysis_notebook.ipynb` — code-only, dependency-ordered, replayable from the
  raw Xenium output — on session save and via the Notebook tab's "Export .ipynb"
  button. Added `nbformat` to `environment.yml`.
- **Notebook tab overhaul** (`tabs/tab_notebook.py`): cells are derived from the
  graph with a ⚠ stale badge, re-running a step updates its cell in place, and a
  **"Show DAG"** button renders the graph (`utils/dag_view.py`, matplotlib +
  networkx) to `plots/provenance_dag.png`. `graph_to_mermaid` / `graph_to_dot`
  provide diagram text.

### Fixed
- Recorded code now actually replays: undefined `fig` in every plot-save snippet
  (`fig = plt.gcf()`), undefined `data_path` (now defined in the preamble), the
  preamble is always emitted, and durable comment-only actions became real code —
  ROI polygons (inlined vertex arrays instead of a cache-only
  `load_rois_from_sdata`), cell-type annotations (CellTypist/LLM/label-transfer
  now emit a real `.map()` producing a `<key>_annotated` column), clustering
  import/export, and the pairwise-volcano loops.
- Leiden HVG/scale branch now copies labels back onto `adata.obs` (previously
  left only on `adata_leiden`); CNV takes raw counts from `sdata['table'].X`
  rather than an already-normalized `adata.X`. `random_state` threaded into
  recorded Leiden/UMAP for determinism.

## [Unreleased] — 2026-07-15 (b)

### Changed
- **CNV tab defaults now match InSituCNV's reference notebook** —
  `smoothing_neighbors` 30→20, `window_size` 10→60, `step` 2→10 in
  `utils/cnv_analysis.py::run_cnv_pipeline()` and the corresponding
  `tabs/tab_cnv.py` spin boxes. The original lower window/step defaults
  were based on a mistaken assumption that a 60-gene window would produce
  zero CNV windows for chromosomes with fewer genes on a small Xenium
  panel. Verified against infercnvpy's actual implementation
  (`_running_mean()`): windowing is computed independently per
  chromosome, and a chromosome with fewer genes than `window_size` isn't
  dropped — it falls back to a single window averaging all of that
  chromosome's genes. So a larger window mainly trades sub-chromosomal
  resolution for a less noisy per-window estimate, which suits
  whole-chromosome/arm-level CNV signal; there was no failure mode being
  avoided by the smaller custom defaults. `n_neighbors`/`lfc_clip` already
  matched the notebook and are unchanged. Added a tooltip and an inline
  hint label to the **CNV cluster resolution** field noting that this
  value (kept at 0.2) may need per-dataset tuning, since InSituCNV's own
  notebook evaluates multiple resolutions rather than recommending one.

## [Unreleased] — 2026-07-15

### Added
- **CNV inference tab (Genes → CNV)** — infers copy-number variation from
  expression data using the [InSituCNV](https://github.com/Moldia/InSituCNV)
  method (`insitucnv` + `infercnvpy`, new optional `cnv` extra:
  `pip install -e ".[cnv]"`). Pick an existing clustering, mark some of its
  categories as the "normal" reference population, and run CNV inference to
  get (a) a new colorable CNV-subclone clustering (`cnv_leiden_<res>`,
  registered into `ctx.clusterings` exactly like every other clustering —
  usable app-wide in Rank Genes, ROI DEG, etc.), (b) a continuous per-cell
  "CNV score" coloring mode, and (c) a chromosome heatmap plot. New
  `utils/cnv_analysis.py` (pure-logic pipeline wrapper), new
  `tabs/tab_cnv.py`, new `CellColorManager.get_continuous_colors()` for
  coloring by arbitrary continuous per-cell scores (not just gene
  expression), new `save_cnv_results_to_adata`/`load_cnv_results_from_adata`
  persistence pair in `utils/adata_persistence.py` (results survive session
  reload, including an `adata_cnv_cache.h5ad` sidecar so the heatmap can be
  regenerated without recomputation). Human genome reference only for now
  (infercnvpy's default GRCh38 gene-position table, auto-downloaded/cached).
  Window/step defaults are set lower than infercnvpy's bulk-RNA-seq defaults
  since Xenium panels are much smaller; both are user-adjustable, and the
  results panel reports how many genes mapped to the genome and how many
  windows were produced so users can judge result quality.

  Installing the `cnv` extra needs two known workarounds, both verified
  against a live end-to-end run: `insitucnv==0.1.0`'s PyPI metadata pins
  `anndata<0.12`/`pandas<3`, which conflicts with this app's own
  `anndata>=0.12` — install it with `pip install --no-deps insitucnv` after
  `pip install -e ".[cnv]"` fails (the pin is stale; insitucnv only uses
  stable AnnData APIs and imports/runs fine against `anndata` 0.13 /
  `pandas` 3.0 in practice). Separately, `run_cnv_pipeline()` now calls
  `_convert_adata_arrow_strings()` before `run_infercnv()` — infercnvpy does
  numpy-style fancy indexing on `var_names` that breaks on pandas 3.0's
  PyArrow-backed string dtype — and patches the removed
  `matplotlib.cm.get_cmap()` API that `insitucnv`'s own cluster-coloring
  code still calls (removed in matplotlib 3.9+, `insitucnv` not yet
  updated).

## [Unreleased] — 2026-07-14 (e)

### Added
- **Crop Dataset now carries over clusterings** — every clustering in
  `ctx.clusterings` (built-in Xenium ones like `graphclust`/`kmeans_*` and
  any custom/Leiden ones), subset to the exported cells, is written into the
  cropped dataset's `adata.obs` using the same `clustering_<name>` /
  `cluster_labels_<name>` column convention already used elsewhere in the app
  (`save_clustering_to_adata`/`save_cluster_labels_to_sdata`), so
  `load_custom_clusterings_from_adata` picks them up automatically on reopen
  — no `analysis/clustering/` folder needed, and no need to recompute
  clustering after cropping. New `_carry_over_clusterings()` helper in
  `utils/crop_export.py`. Verified against a real dataset (10 built-in
  clusterings + 1 synthetic custom one + cluster labels, all round-tripped
  with zero mismatches for the kept cells).

## [Unreleased] — 2026-07-14 (d)

### Fixed
- **Datasets with no clusterings crashed at startup** — a dataset with an empty
  clustering list (e.g. a Crop Dataset export, which has no raw Xenium
  `analysis/clustering/` folder) crashed while building the control panel:
  `ValueError: None is not a valid choice. must be in ()`, raised by every
  "Clustering" magicgui `ComboBox` when handed `value=None` against empty
  choices. Added a shared `combo_value_kwargs()` helper
  (`tabs/_helpers.py`) that omits the `value` kwarg when the choice list can't
  supply the requested index (magicgui then safely defaults to the first
  choice, or `None` when empty), and applied it to every clustering selector
  (Cell Coloring, Rank Genes, Markers, Lig-Rec, Nhood Enrich, Co-occur, Annot
  Nhood, Annot Dist) plus the Gene Correlation gene selectors. Such datasets
  now open normally, with the clustering dropdowns simply empty until a
  clustering is computed.

## [Unreleased] — 2026-07-14 (c)

### Fixed
- **Crop Dataset: exported datasets failed to open with `FileNotFoundError:
  ... projection.csv`** — `load_umap()`/`load_clusterings()` unconditionally
  read `analysis/umap/.../projection.csv` and `analysis/clustering/`, which a
  Crop Dataset export never has (no raw Xenium `analysis/` folder). Both now
  tolerate a missing `analysis/` folder (return `None`/`{}` instead of
  raising). `app.py`'s dataset loader also now reconstructs the UMAP tab's
  coordinates from `adata.obsm['X_umap']` when the raw CSV is absent but the
  table already carries embedded UMAP coordinates from the source dataset —
  so a cropped dataset's UMAP view still works if the source had one computed
  at crop time. Verified against a real cropped dataset via the full
  `app._load_dataset()` reload path.

## [Unreleased] — 2026-07-14 (b)

### Fixed
- **Crop Dataset: exported datasets failed to open with `--no-cache`** —
  `load_sdata()` unconditionally rebuilt from raw Xenium files whenever
  `use_cache=False`, but a Crop Dataset export has no raw files (by design —
  it's a lightweight, self-contained zarr package), so this crashed with
  `FileNotFoundError: ... cells.zarr.zip`. Now `load_sdata()` falls back to
  the zarr cache whenever raw Xenium files aren't present, regardless of the
  `use_cache`/`--no-cache` flag, since that's the only way such a directory
  can ever be loaded. Fixes already-exported crop datasets too, not just new
  ones, since the check is based on file presence rather than anything
  written at export time.

## [Unreleased] — 2026-07-14 (a)

### Fixed
- **Crop Dataset: orphan nuclei in `nucleus_labels`** — nucleus IDs are their own
  independent numbering, unrelated to `cell_id`/`cell_labels`, so an ID-based overlap
  check (added, then removed, in the previous fix) couldn't reliably mask them and
  the fallback left `nucleus_labels` completely unmasked. Now masked correctly via
  spatial overlap with the already-masked `cell_labels` crop: a nucleus is kept only
  if it occupies at least one pixel within a kept cell's footprint.

## [Unreleased] — 2026-07-13

### Added
- **Crop Dataset (Tools tab)** — draw one or more polygons in a new "Crop Regions"
  napari layer and export each as its own standalone, independently-openable
  xenium-viewer data directory (cropped morphology image, cell/nucleus labels,
  transcripts, and AnnData table). Images/labels are cropped to each polygon's pixel
  bounding box; cells and transcripts are filtered to the exact drawn polygon via
  true point-in-polygon tests. Output/name for each region is chosen via sequential
  folder-picker + name-prompt dialogs, run in a background thread with a progress
  dialog. New module `src/xenium_viewer/utils/crop_export.py`; new tab
  `src/xenium_viewer/tabs/tab_crop_dataset.py`.

## [Unreleased] — 2026-07-01

### Fixed
- **Permission-denied dialog gaps** — two write paths were missing coverage from the
  read-only zarr dialog: (1) `delete_element_from_disk` in the Patches tab only printed
  to console on failure; (2) `record_code` in `_helpers.py` had no exception handling at
  all when writing `code.py`. Both now call `_maybe_show_permission_dialog` on
  `PermissionError`/`OSError`.

## [Unreleased] — 2026-06-30

### Added
- **tab10 palette in Patches tab** — `tab10` (matplotlib's 10-colour categorical
  palette) is now available in the Patches tab palette dropdown alongside tab20,
  glasbey_dark, Set1, Set3, and ARMS.

## [Unreleased] — 2026-06-26 (d)

### Fixed
- **Leiden clustering UI freeze** — `sc.tl.leiden` held the Python GIL during graph
  construction and partition, causing the progress bar and status-bar spinner to freeze.
  The Leiden step now runs in a subprocess via `ProcessPoolExecutor` (spawn context),
  giving it its own GIL so the main process's Qt event loop stays responsive throughout.
  New module: `src/xenium_viewer/utils/leiden_worker.py`.

## [Unreleased] — 2026-06-26 (c)

### Added
- **Progress bar for long-running analyses** — an indeterminate `QProgressBar` now
  appears inside the control panel directly below the run button while an analysis is
  in progress (Leiden clustering, rank genes, L-R, neighbourhood enrichment,
  co-occurrence, ROI DEG, annotation nhood enrichment, Novae domains). The bar
  disappears automatically when the analysis finishes or errors out. The existing
  napari status-bar spinner/tqdm text is retained unchanged.

## [Unreleased] — 2026-06-26 (b)

### Added
- **Wiki screenshots** — all 52 documentation screenshot placeholders are now filled with
  actual PNGs captured from a running viewer instance. A new script
  `scripts/capture_screenshots.py` automates future recapture by programmatically
  navigating each tab and grabbing the control-panel and full-window views.
  One placeholder (`tutorial-clustering-step5.png`, the matplotlib dotplot window)
  remains a comment for manual capture.

## [Unreleased] — 2026-06-26

### Fixed
- **Zarr cache rebuild no longer silently destroys user data** — when `experiment.xenium`
  is newer than `sdata_cached.zarr` (e.g. after Xenium Explorer opens the dataset, or a
  backup/rsync resets timestamps), a Qt dialog now appears listing all user-generated data
  found in the cache (ROIs, H&E/ARMS landmarks and images, clusterings, analysis results,
  etc.) and offers three choices: *Rebuild and restore my data* (backs up the old cache,
  rebuilds from raw files, merges user elements back, then deletes the backup),
  *Rebuild without restoring* (previous behaviour), or *Keep existing cache* (skip
  rebuild and load from the stale cache as-is — safe when only metadata changed).
  When the cache is unreadable (corrupt), it is now preserved as
  `sdata_cached_corrupt_<timestamp>.zarr` instead of being silently deleted, so data can
  be recovered manually.

## [Unreleased] — 2026-06-25

### Fixed
- **ROI DEG region grouping** — `compute_roi_deg` was reading `uns['rank_genes_groups']`
  (the scanpy default key) instead of the key it had just written (`uns['wilcoxon']` /
  `uns['t-test']`). Because the subset AnnData is a copy of `ctx.adata`, it could
  inherit stale cell-type DEG results from a previous analysis, causing "Run ROI DEG"
  to display cluster-level rather than region-level results. Fixed by passing `key=method`
  to `sc.get.rank_genes_groups_df`.
- **Permission errors on read-only datasets** — write failures (clustering import, H&E
  registration, session save, etc.) previously silenced to the console only. A
  `QMessageBox` warning dialog now fires on the first `PermissionError` / read-only
  `OSError` of a session, telling the user to copy the dataset to a writable location
  or launch with `--no-cache`.

### Added
- **ARMS Overlay: save/load landmarks** — "Save Landmarks..." and "Load Landmarks..."
  buttons added to the ARMS Overlay tab (below "Compute Registration"), mirroring the
  H&E Registration tab. Landmarks and the computed affine are saved to a portable JSON
  file via the existing `save_landmarks` / `load_landmarks` API. The save button
  enables when ≥1 landmark pair is placed and disables on "Clear All".

### Changed
- **H&E Registration: removed Save Affine button** — the "Save Affine..." button has
  been removed. There was no corresponding "Load Affine..." button, making it a dead end.
  The affine is already auto-persisted to the sdata zarr cache on every registration.

## [Unreleased] — 2026-05-05

### Changed
- **Installable Python package** — the repo is now packaged with `pyproject.toml`
  (hatchling build backend) and shipped under `src/xenium_viewer/`. The three
  numbered entry-point scripts have been renamed:
  `00_preprocess_transcripts.py` → `xenium_viewer/preprocess.py`,
  `01_load_sdata.py` → `xenium_viewer/loader.py`,
  `02_xenium_viewer.py` → `xenium_viewer/app.py`. Cross-imports were rewritten
  from `from utils.X` / `from tabs.X` to `from xenium_viewer.utils.X` /
  `from xenium_viewer.tabs.X`.
- **Console-script entry points** — `xenium-viewer`, `xenium-preprocess`,
  `xenium-fetch-references`, and `xenium-build-custom-segmentation` are
  installed on `PATH`.
- **Conda environment file** — `environment.yml` reproduces a working env in
  one step (`conda env create -f environment.yml`); pypi-only deps and the
  package itself (`-e .`) are listed in the embedded `pip:` block. Optional
  extras (`celltypist`, `r`, `gpu`, `references`, `full`) are exposed via
  `[project.optional-dependencies]`.
- **Removed importlib hacks** — `app.py` no longer loads its sibling modules
  via `importlib.util.spec_from_file_location`; they are normal package
  imports. The `sys.path.insert` at module top is gone.

## [Unreleased] — 2026-04-18

### Changed
- **External Images tab: single composite layer** — multichannel images (e.g. 25-channel PhenoCycler)
  now display as a single RGB composite napari layer instead of one layer per channel. Per-channel
  visibility checkboxes, color buttons, and contrast range sliders are in the tab widget.
  Composite is built lazily via dask so large images render efficiently.
- **External Images tab: landmark registration** — external images can now be independently aligned
  to the Xenium image via landmark-based registration (same `compute_landmark_affine` as H&E tab).
  Includes flip V/H checkboxes. The "Apply transform from" dropdown remains as an alternative.
  Landmarks are persisted to sdata.shapes; affine to sdata transformations.

## [Unreleased] — 2026-04-17

### Fixed
- **Overlay affine persistence** — affine transforms for patch overlays and external images are now
  saved to the SpatialData object via `set_transformation()` / `write_transformations()` (same
  pattern as H&E registration). On restore, the affine is read back from sdata directly, so
  overlays are correctly positioned even before the source layer (e.g. H&E) finishes loading.
  A deferred-linking listener also re-establishes live affine mirroring once the source layer
  appears.
- **QCheckBox crash on close** — `_snapshot_layers` now reads `entry["hidden_cluster_ids"]` (a
  plain set maintained by the tab) instead of iterating Qt checkbox widgets that may already be
  destroyed during shutdown.

### Changed
- **ARMS palette for subclone predictions** — subclone prediction overlays now default to the
  ARMS palette (RColorBrewer Set1+Set2+Dark2, 24 colours) to match the R-based ARMS visualisation
  pipeline. Cluster-to-colour mapping is 0-normalised so 1-based genomic cluster IDs (1,2,3) map
  to Set1 red, blue, green — matching the R package exactly.

## [Unreleased] — 2026-04-16

### Added
- **External Images tab** — load arbitrary multichannel OME-TIFF/TIFF/SVS files (e.g. PhenoCycler)
  as lazy-loaded napari layers. Channel axis is detected from `tif.series[0].axes` / OME-XML so
  IF images are never misinterpreted as RGB. Each channel gets its own napari sub-layer with
  native contrast / colormap / visibility controls; the tab adds group opacity, show/hide-all,
  and an "Apply transform from…" selector that live-mirrors another layer's affine (updates when
  the source is re-registered). Pixels and affine are written to `sdata.images["ext_<slug>"]`;
  only per-entry UI state (opacity, affine source) lives in zarr attrs.
- **Patch Overlays tab** — load phikon patch-cluster outputs (folder with `patches/coordinates.npy`
  + `clustering/cluster_labels.npy`) and subclone-prediction CSVs (`x_coord`, `y_coord`,
  `predicted_genomic_cluster`, `morphology_cluster`, `prediction_confidence`). Patches render as
  polygon rectangles in a napari Shapes layer (scales to 10k+ patches). Controls for cluster
  column, palette (tab20 / glasbey_dark / Set1 / Set3 / ARMS), outline-only, edge width, opacity,
  confidence threshold, and per-cluster visibility. Patch size is inferred from folder / file name
  with a stride-based sanity check and confirmation dialog. Geometry, cluster columns, and
  confidence are stored in `sdata.shapes["patch_<slug>"]`; UI settings persist in zarr attrs.
  Subclone prediction overlays default to the ARMS palette (RColorBrewer Set1+Set2+Dark2) to
  match the R-based ARMS visualisation pipeline.
- **Affine mirroring helper** — `utils/affine_linking.py` lists layers with non-identity affine
  and wires up `events.affine` subscriptions so external overlays stay aligned when the source
  image is re-registered.

## [Unreleased] — 2026-04-13

### Added
- **Cluster labels persisted to SpatialData** — user-assigned cluster names (from the label editor,
  reference atlas annotation, and LLM annotation) are now saved as `adata.obs["cluster_labels_<key>"]`
  columns (e.g. `cluster_labels_leiden_r1.0`) inside `sdata.tables["table"]` immediately on
  assignment. Labels are loaded back on the next launch and merged with the session-attrs fallback
  (sdata wins on conflict). The obs columns are readable in any standalone Python session.

## [Unreleased] — 2026-04-12

### Changed
- **ROI DEG and ARMS tile DEG persistence migrated to SpatialData** — results are now saved as
  sidecar parquets inside the zarr cache (`roi_deg_cache.parquet`, `arms_tile_deg_cache.parquet`)
  immediately on computation rather than at session close. Restores automatically on relaunch.
  One-time migration copies old `viewer_session/*.parquet` files to the new location.
- **ARMS tile DEG code recording fixed** — generated code snippet in `code.py` is now fully
  executable: loads tile polygons from `sdata` via `load_arms_tiles_from_sdata`, reconstructs the
  ARMS registration affine (fine × flip), applies it to the tiles, filters by selected tile
  clusters, and optionally applies a Xenium cell cluster mask before calling `compute_arms_tile_deg`.

## [Unreleased] — 2026-04-10

### Added
- **Custom segmentation SpatialData persistence** — custom segmentations are now cached inside the
  SpatialData zarr store (`sdata.labels["custom_cell_labels"]` + `sdata.tables["custom_table"]`)
  after first load, so subsequent opens do not require re-selecting the h5ad file.
  - On "Load Custom Segmentation...": if a cached version is detected in sdata, a dialog offers to
    load from cache (fast) or select a new h5ad file.
  - Auto-restore on launch: if custom segmentation was active when the dataset was closed,
    `_restore_session` reloads it from sdata automatically.
  - All per-action saves (clustering, rank genes, nhood enrichment, etc.) now route to
    `sdata.tables["custom_table"]` when custom segmentation is active, so analysis results are
    preserved across sessions.
  - `session.py` now persists `segmentation_source` ("xenium" | "custom") in the viewer session attrs.
  - **"Update SpatialData on disk"** button in the Segmentation tab force-syncs current in-memory
    state to sdata (custom table + labels in custom mode; xenium table in xenium mode).

## [Unreleased] — 2026-03-26 (2)

### Added
- **Custom cell segmentation** — two-stage pipeline to replace the native Xenium segmentation
  with a custom one from a Seurat v5 RDS file.
  - `scripts/extract_seurat_segmentation.R` — Stage 1: extracts polygon vertices, count matrix,
    and cell metadata from a Seurat RDS via slot access (no Seurat installation required). Outputs
    `segmentation_polygons.csv`, `counts.mtx`, `genes.txt`, `barcodes.txt`, `cell_metadata.csv`.
  - `scripts/build_custom_segmentation.py` — Stage 2 (run in `xenium_viewer` conda env): reads
    Stage 1 outputs, validates coordinate system against `cells.parquet`, rasterizes 300K+ polygons
    using `rasterio.features.rasterize`, writes a multi-scale zarr label store
    (`custom_labels.zarr`) and an AnnData file (`custom_segmentation.h5ad`) with
    `obs['cell_id']` = integer label and `obsm['spatial']` = xy-µm centroids, plus a metadata JSON.
  - **Segmentation tab** — new viewer tab with "Load Custom Segmentation..." and "Revert to Xenium
    Segmentation" buttons. Load flow: select `custom_segmentation.h5ad` → zarr found automatically
    via JSON metadata → cell labels layer replaced, `ctx.adata`, `ctx.color_manager`,
    `ctx.centroids_yx`, `ctx.clusterings` all rebuilt from the new data, clustering ComboBoxes
    refreshed across all tabs, stale analysis caches cleared. Revert restores the native
    spatialdata labels and reloads the original clustering files.
- `segmentation_source: str` field on `ViewerContext` (values: `"xenium"` | `"custom"`).

## [Unreleased] — 2026-04-01

### Added
- **UMAP save** (`tabs/tab_umap.py`) — "Save UMAP Plot..." button saves a scanpy-native `sc.pl.umap`
  figure in PNG or SVG. Uses the current cell-coloring clustering; cluster labels (if set) are
  applied as category names so they appear in the legend and on-data annotation.
- **Marker Genes tab** (`tabs/tab_marker_genes.py`) — new tab accepting a JSON marker-genes dict
  (`{"Cell type": ["Gene1", "Gene2"]}`). Generates dotplot, heatmap, matrixplot, tracksplot, and
  correlation matrix using scanpy's native plotting functions with the chosen clustering. Format
  (PNG/SVG) is selectable. The marker dict is persisted to the zarr session and restored on startup.
- **Volcano plot (revised)** — `make_volcano_plot()` redesigned with EnhancedVolcano aesthetics:
  only genes passing **both** thresholds are coloured (up-regulated red `#DC0000`,
  down-regulated blue `#4DBBD5`); everything else is grey. x-axis is symmetric around zero with
  auto-computed 99th-percentile limits. Points outside the display range are shown as directional
  triangle markers (`>`, `<`, `^`) pinned to the axis edges. Default labelled genes increased from
  10 to 20. Uses seaborn `despine` for theme_classic look and `adjustText` for label repulsion
  (requires `pip install adjusttext`; falls back to plain text if not installed).
- **ROI DEG: volcano plots** — "Save Volcano Plot(s)..." button in the ROI tab generates pairwise
  volcano plots after running ROI differential expression. With 2 ROIs, 1 plot is saved; with N ROIs,
  all C(N, 2) pairwise plots are saved as `roi_volcano_Region_X_vs_Region_Y.png` (300 dpi).
  Thresholds: adjusted p-value < 0.01 and |log2FC| > 1. Significantly up-regulated genes are red,
  down-regulated are blue, with dashed threshold lines. Progress shown in status bar.
  - `compute_roi_deg()` now returns `(df, adata_norm)` tuple so the normalized subset AnnData is
    available for pairwise comparisons (matching the existing `compute_arms_tile_deg()` signature).


- **ARMS tiles: "Outline only" checkbox** — when checked, tiles are rendered as colored outlines
  (edge = cluster color, fill = transparent) instead of filled polygons. The checkbox is enabled
  once tiles are loaded and toggles the existing layer in place without reloading.
- **ARMS tiles: "Tile edge width" slider** — adjusts outline thickness live (range 1–100, default
  20). Enabled once tiles are loaded.

### Changed
- **ARMS tiles: ColorBrewer Set1+Set2 color palette** — replaced the hardcoded 8-color custom
  palette with a concatenation of ColorBrewer Set1 (9 colors) + Set2 (8 colors) = 17 distinct
  colors. Cluster filter checkboxes and the legend now show `C{id}` labels without color names.

## [Unreleased] — 2026-03-27 (3)

### Fixed
- **Annot. Nhood and Annot. Distance clustering dropdowns not refreshed on segmentation swap** —
  both tabs created their `clustering_widget` ComboBoxes without registering them on `ctx`, so
  `refresh_clustering_choices` skipped them. Fix: register as `ctx.annot_nhood_clustering_widget`
  and `ctx.annot_dist_clustering_widget`, add corresponding fields to `ViewerContext`, and include
  them in the refresh loop in `_helpers.py`.

## [Unreleased] — 2026-03-27 (2)

### Fixed
- **Custom segmentation: cluster coloring all-transparent** — obs columns (`cell_type`, etc.)
  added to `ctx.clusterings` from `new_adata.obs` were indexed by obs row numbers (`'0'`, `'1'`,
  …) rather than `cell_id` integer label values. `get_cluster_colors` calls
  `cluster_series.reindex(cell_ids)` using label values as the lookup key, so the index mismatch
  produced all-NaN alignments and every cell got RGBA `(0,0,0,0)`. Fix: wrap obs columns in a new
  `pd.Series` indexed by `adata.obs['cell_id'].values` before storing in `ctx.clusterings`.
  Gene expression coloring was unaffected because it uses `label_to_obs` directly.
- **Custom segmentation: `IndexError` in UMAP hover** — `UMAPViewer._valid` has shape
  `(n_original_cells,)` but `cluster_ids_per_obs` from a custom segmentation has a different
  length. Guard added: skip UMAP hover cluster IDs when sizes don't match.

## [Unreleased] — 2026-03-27

### Added
- **Annotation layer** — new napari Shapes layer "Annotations" for drawing named tissue regions
  (bone, adipocyte, vessel, etc.) with per-shape type labels displayed as text overlays.
  Annotations persist across sessions via `sdata.shapes['annotations']`.
- **Annotations tab** — manage annotation type labels (assign, colour-pick, delete), import from
  GeoJSON (compatible with QuPath exports), export to GeoJSON.
- **Annot. Nhood tab** — neighbourhood enrichment analysis that includes annotation regions as
  virtual cell types. Samples a configurable-density grid of virtual cells inside annotation
  polygons, builds an augmented AnnData, and runs the squidpy neighbourhood enrichment pipeline.
  Annotation types appear as rows/columns in the Z-score heatmap.
- **Annot. Distance tab** — for each real cell computes minimum distance (µm) to the boundary of
  a selected annotation type (using shapely's vectorised distance API). Visualises the distribution
  per cell-type cluster as violin, box, or strip plots. Cells can optionally be coloured on the
  canvas by their distance value using a Points layer.

## [Unreleased] — 2026-03-26

### Added
- **View menu: Show Minimap toggle** — checkable "Show Minimap" action appended to napari's native
  View menu. Enabled and checked when the minimap overlay is available; disabled (grayed out) when
  there is no morphology data. State is reset on dataset reload.

## [Unreleased] — 2026-03-25

### Changed
- **Phase 3 SpatialData storage refactoring** — H&E and ARMS landmark pairs, and ARMS tile
  geometries migrated from custom zarr arrays into native SpatialData shapes elements, making
  the sdata object self-sufficient for Python analysis.
  - **H&E landmarks** → `sdata.shapes['he_xenium_landmarks']` and `sdata.shapes['he_he_landmarks']`
    as GeoDataFrames of shapely Points (xy coords in native pixel spaces). Saved after each
    registration event; cleared on landmark clear.
  - **ARMS landmarks** → `sdata.shapes['arms_xenium_landmarks']` and `sdata.shapes['arms_he_landmarks']`
    — same pattern as H&E landmarks.
  - **ARMS tile polygons** → `sdata.shapes['arms_tiles']` as a GeoDataFrame of shapely Polygons
    with `tile_name` and `cluster_id` columns. Saved after GeoJSON+CSV load; restored at startup
    without requiring the original files.
  - UI-only state (flip toggles, coarse_affine, file paths) remains in zarr attrs as before.
  - `utils/session.py` no longer writes landmark zarr arrays; zarr load paths kept as migration
    fallback for existing datasets.
  - 4 new functions in `utils/adata_persistence.py`: `save_landmarks_to_sdata`,
    `load_landmarks_from_sdata`, `save_arms_tiles_to_sdata`, `load_arms_tiles_from_sdata`.
  - **Automatic migration**: on first launch after upgrade, `migrate_landmarks_to_sdata()`
    migrates existing landmarks into sdata.shapes. Sources checked in order: (1) zarr session
    arrays (snapshot-captured at close time), (2) `landmarks.json` in the dataset folder
    (saved via Save Landmarks button). Also re-parses GeoJSON/CSV tile files into
    sdata.shapes['arms_tiles'] if the files still exist. Marked with
    `migrated_landmarks_to_sdata` flag so it only runs once.
  - All shape GeoDataFrames now constructed via `ShapesModel.parse()` to ensure spatialdata
    compatibility (attaches coordinate transform metadata required by the zarr writer).
- **Phase 2 SpatialData storage refactoring** — rank genes results, normalized AnnData, and
  ROI polygons migrated out of custom `viewer_session/` files into native SpatialData locations.
  - **Rank genes** → `adata.uns['rank_genes_groups']` (scanpy native); DataFrame reconstructed
    via `sc.get.rank_genes_groups_df()` on restore. Removes `rank_genes.parquet`.
  - **Normalized AnnData** → `sdata.tables['adata_norm']` (zarr-backed); replaces the
    `rank_genes_adata_norm.h5ad` file. Required for dotplot and volcano plots after restore.
  - **ROI polygons** → `sdata.shapes['rois']` as a GeoDataFrame of shapely Polygons
    (xy coords); replaces per-polygon zarr arrays in `viewer_session/rois/`. Old zarr arrays
    serve as fallback for datasets not yet saved with the new code.
  - Rank genes results now enable Show/Export buttons at startup (loaded before `restore_fn`,
    same pattern as Phase 1 analyses).
  - Old `rank_genes.parquet` + `rank_genes_adata_norm.h5ad` files are auto-migrated and
    deleted on first launch with the new code.
  - `utils/session.py` no longer handles rank genes or ROI polygon storage.
- **Phase 1 SpatialData storage refactoring** — clusterings, nhood enrichment, co-occurrence,
  L-R interaction results, and UMAP coordinates are now stored in native AnnData locations
  (`adata.obs`, `adata.obsm["X_umap"]`, `adata.uns`) instead of custom zarr/parquet files in
  `viewer_session/`. Persisted via `sdata.delete_element_from_disk` + `sdata.write_element`.
  - Analysis results are now saved immediately after computation (crash-safe), rather than
    deferred to viewer exit.
  - Old `viewer_session/` data is automatically migrated on first launch with the new code.
  - New module: `utils/adata_persistence.py` centralizes all adata read/write/migration logic.
  - `utils/session.py` no longer handles clusterings, nhood, co-occurrence, or L-R data.

### Fixed
- **Index alignment for adata persistence** — `adata.obs.index` uses integer strings (`'0'`, `'1'`, ...)
  while clustering series and UMAP are indexed by cell barcode (`'aaaagflk-1'`, ...). Save/load now
  maps between the two via the `cell_id` column.
- **Show/Export buttons disabled at startup for persisted analyses** — nhood enrichment,
  co-occurrence, and L-R results were loaded from `adata.uns` *after* tab session restore ran,
  so the buttons stayed disabled despite results existing. Fixed by loading `adata.uns` analysis
  results before calling `restore_fn()`.

## [Unreleased] — 2026-03-23

### Added
- **Level slider in Novae tab** — exposes the `level` parameter (hierarchical tree level,
  default 7, range 1–15) for `assign_domains()`, giving finer control over domain granularity.
- **Console variable injection** — key variables (`adata`, `sdata`, `viewer`, `ctx`,
  `clusterings`, `color_manager`, `gene_names`, `data_path`) are now pushed into napari's
  built-in IPython console on dataset load, with a help message listing what's available.
  Variables are refreshed on dataset reload. Enables a GUI-to-code handoff workflow alongside
  the existing code recording feature.

### Fixed
- **Console button hidden by layer widgets** — capped layer controls and layer list dock
  widgets to 200px max height so the console toggle button stays visible.
- **Console not resizable when opened** — the Xenium Controls dock (QTabWidget with 13 tabs)
  had a ~970px minimum height from stacked widget content, leaving no room for the console.
  Fixed by wrapping tab content in `QScrollArea` inside `make_tab()`, dropping the dock's
  minimum height to ~117px. Also defers the `resizeDocks` call via `QTimer.singleShot(0)` to
  ensure it runs after Qt finishes the visibility layout pass.

## [Unreleased] — 2026-03-20

### Added
- **Spatial Domains tab (Novae)** — new tab that runs [Novae](https://mics-lab.github.io/novae/)
  zero-shot spatial domain inference. Select species (human/mouse), optionally specify N domains
  (0 = auto-detect), and click "Run Novae Domains". On completion the cell labels layer is
  automatically recolored by inferred domains, the `novae_domains` key is added to the clustering
  dropdowns, and results are persisted to the session cache so they are restored on relaunch.
  Requires `pip install novae`. The full pipeline is recorded to `code.py`.

## [Unreleased] — 2026-03-19

### Added
- **Session-scoped code file** — each viewer session now writes to a timestamped file
  (`code_YYYYMMDD_HHMMSS.py`) instead of overwriting `code.py`. Opening a second dataset
  creates a fresh timestamped file; prior session files are preserved.
- **"Continue from existing code file" in Preferences** — lets users resume a previous
  session's code journal. Opens a file dialog, warns about potential duplication, then
  re-detects preamble/normalize/clustering tags from the file content so subsequent actions
  do not emit duplicate boilerplate.
- **Cluster label mapping recorded as executable dict** — editing cluster labels via the
  label editor now emits a `cluster_labels_<key> = {...}` dict to the code journal instead
  of a generic comment, making the output directly usable in downstream plotting calls.

### Changed
- **Auto-save plots on Show** — clicking "Show" in the Neighborhood Enrichment,
  Co-occurrence, Ligand-Receptor, Gene Correlation, and Gene Analysis (dotplot) tabs now
  automatically saves the figure to `<data_dir>/plots/<stem>.<format>` and reports the
  path in the status bar. The separate "Save Plot" button has been removed from all tabs.
- **Default plot format is now SVG** — plots are saved as vector graphics out of the box;
  PNG remains available via Preferences → Plot Format.

### Added
- **Gene Correlation tab** — scatter plot of per-cell expression for any two selected
  genes, annotated with Pearson r and Spearman ρ (both with p-values). Optionally
  restricted to the current cluster filter selection. Plot can be saved in the
  user's preferred format (PNG/SVG).
- **Normalisation options for Gene Correlation** — new "Normalisation" ComboBox with three
  choices: *Raw counts* (default off), *Fraction of total* (library-size correction), and
  *Log1p(CPM)* (default; matches normalisation used in DEG / dotplot tabs). Axis labels and
  status-bar output reflect the selected mode. `code.py` emits the correct normalisation
  snippet for each choice.
- **Filter small clusters by cell count** — new "Min cluster size" slider (100–10,000)
  and "Filter Small Clusters" button in the Cell Colouring tab. Clicking the button
  auto-unchecks clusters with fewer than N cells and enables the cluster filter,
  propagating to all downstream analyses (cell colouring, ligand-receptor,
  neighborhood enrichment, co-occurrence).

### Fixed
- **Rank genes results lost on crash** — rank genes results are now auto-saved immediately
  after computation via `save_rank_genes_incremental` in `utils/session.py`. This writes
  `rank_genes.parquet`, `rank_genes_adata_norm.h5ad`, and the `rank_genes_groupby` key
  to the zarr session. On session restore, all downstream buttons (dotplot, rank plot,
  volcano, edit/reset labels, export, LLM annotate) are re-enabled and the results
  preview is repopulated. Previously, only export and volcano were re-enabled and
  dotplot/rank-plot required re-running rank genes. Skipped when `--no-cache` is active.
- **Custom clusterings lost on crash** — Leiden runs and CSV imports now write an
  incremental parquet + update `custom_clustering_names` in the zarr session immediately
  after creation (`save_clustering_incremental` in `utils/session.py`). Previously,
  clusterings were only persisted on clean exit; a force-kill or crash meant they were lost.
  Skipped automatically when `--no-cache` is active or the zarr store does not yet exist.
- **Duplicate File/Preferences menus on dataset reload** — `create_file_menu()` now appends
  Xenium actions to napari's native `viewer.window.file_menu` (called once at startup) instead
  of creating a parallel custom "File" QMenu on every dataset load. `create_preferences_menu()`
  now finds-or-creates a single "Preferences" menu and clears stale actions before repopulating,
  so N dataset loads no longer produce N Preferences menus. The now-redundant
  `_attach_bare_file_menu()` helper and its cleanup code in `_do_full_init()` have been removed.
- **Dotplot crash with `'NoneType' object has no attribute 'get_axes'`** — in a
  `@thread_worker` context `plt.gcf()` returns `None`; the `_relabel_axes` fallback
  now checks `.fig` and `.figure` attributes explicitly and skips relabelling rather
  than crashing when no figure object is found.
- **Dotplot/Rank Genes Plot crash with duplicate cluster labels** — when multiple clusters
  share the same user-assigned label (e.g. several clusters all named "Unknown"), the old
  code called `rename_categories()` which raised `ValueError: Categorical categories must
  be unique`. Fixed by removing the rename + re-run approach entirely: plots now use the
  original integer cluster IDs from `uns['rank_genes_groups']` and apply label substitution
  post-hoc via axis tick replacement (`_relabel_axes`).
- **Minimap persists on dataset reload** — when opening a new dataset via File > Open Dataset,
  the old `MinimapWidget` is now explicitly hidden and scheduled for deletion before the new
  dataset is loaded, preventing stale thumbnails from overlapping the new minimap.
- **Cluster label dialog sorts numerically** — the "Edit Cluster Labels" dialog in the
  Clustering and Gene Analysis tabs now sorts cluster IDs numerically (0, 1, 2, … 19) instead
  of lexicographically (0, 1, 10, 11, … 19, 2, 3); falls back to string sort for non-integer IDs.

### Changed
- **"Edit Cluster Labels..." moved to Cell Colouring tab** — the button is now in the Cell
  Colouring tab (below "Apply Cell Coloring"), where it is more naturally discovered alongside
  the clustering selector and cluster filter controls. Removed from the Clustering tab.

### Added
- **LLM-based cluster annotation** — new "LLM Annotation" section in Gene Analysis tab.
  After running Rank Genes, click "Annotate with LLM" to send top 10 marker genes per
  cluster to a locally installed AI CLI (Claude, Gemini, or Codex). The LLM returns cell
  type annotations as JSON, which are stored as cluster labels and shown in dotplots,
  rank gene plots, and hover text. Runs in a background thread to keep UI responsive.

- **Cluster ID shown on hover** — when hovering over cells with custom labels (from
  CellTypist, label transfer, or LLM annotation), the status bar now shows both the
  numeric cluster ID and the label name, e.g. "Cell 12345 — leiden: 3 (T cells)".

- **Reference dataset auto-discovery** — the label transfer section now discovers
  reference datasets from `.metadata.json` sidecar files in `reference_datasets/`
  instead of using a hardcoded dictionary. Each h5ad gets a companion JSON file with
  structured metadata (paper, authors, journal, year, platform, tissue, annotation
  columns, default column, annotation workflow). Selecting a dataset in the ComboBox
  shows paper/platform/tissue info in the status bar. New datasets can be added by
  placing an h5ad + metadata.json pair in `reference_datasets/`.

- **`scripts/fetch_reference_datasets.py`** — standalone download & conversion script
  for 4 additional public scRNA-seq reference datasets: Prostate Cell Atlas
  (Clatworthy/Tuong 2021, direct h5ad), HuPSA (Cheng et al. 2024, Seurat .rds via
  rpy2), GSE176031 (Song/Huang 2021, GEO TXT matrices), GSE181294 (Mei et al. 2023,
  GEO MTX/CSV). Run with `--all` or `--dataset <name>`. Supports `--force` re-download,
  `--keep-raw` to preserve staging files. Each dataset generates both an h5ad file and a
  `.metadata.json` sidecar.

- **Label transfer via sc.tl.ingest()** — new "Label Transfer" section in Gene Analysis
  tab. Select a reference scRNA-seq dataset (two preconfigured prostate cancer datasets
  included, or browse for any h5ad file), load it, pick an annotation column, and click
  "Run Label Transfer" to project cell type labels onto Xenium clusters. Pipeline: find
  common genes → subset both → preprocess reference (normalize, log1p, HVG, PCA, neighbors,
  UMAP) → preprocess Xenium (normalize, log1p) → `sc.tl.ingest()` → majority vote per
  cluster. Reference datasets are cached in memory after first load. Preconfigured datasets:
  Zhao et al. (222K cells, ct.main/ct.sub/ct.sub.epi) and eBioMedicine (21K cells,
  ID/ID_coarse/refinedID). New utility functions: `load_reference_h5ad()`,
  `get_annotation_columns()`, `run_label_transfer()` in `utils/gene_analysis.py`.

- **CellTypist automated cell type annotation** — new "CellTypist Annotation" section
  in the Gene Analysis tab. Select a CellTypist model from a dropdown, click "Annotate
  with CellTypist" to run per-cell predictions, then majority-vote assigns cell type
  labels to each cluster. Labels are stored in `state["cluster_labels"]` (same system
  as manual labels), so dotplots, rank genes plots, hover info, and exports all pick
  up the new names automatically. Includes "Download Models" button for one-time model
  download, and a **confidence threshold slider** (default 0.5) that filters out
  low-confidence per-cell predictions before the majority vote — prevents immune-focused
  models from labeling non-immune clusters as "T cells". Status bar shows how many cells
  passed the threshold. Gracefully degrades when celltypist is not installed (widgets
  shown but disabled). New `run_celltypist_annotation()` utility in
  `utils/gene_analysis.py`.

- **Reset Labels button** in Gene Analysis tab — quickly reverts cluster labels back to
  original numeric IDs after CellTypist annotation or manual editing. Removes the entry
  from `state["cluster_labels"]` for the selected clustering.

- **Complete code recording coverage** — all analysis actions now record parameters to
  `code.py` via `ctx.record_code()`:
  - **Transcripts**: QV threshold in overlay recording; new density heatmap recording
    (gene, bin_size, normalise, cluster_filter)
  - **Ligand-receptor**: interaction database selection and CellPhoneDB flag
  - **ROI DEG**: cluster filter status and selected clusters
  - **UMAP**: display action recorded as comment
  - **Cell coloring**: background color toggle
  - **H&E registration**: flip state changes
  - **ARMS overlay**: flip state changes; Xenium cluster filter in tile DEG
  - **Co-occurrence**: cluster subset in run recording; filter_targets in plot recording
  - **Cluster labels**: label editing recorded in `_helpers.py`

### Fixed
- **Progress feedback not showing** — two bugs prevented tqdm progress from appearing in
  the status bar during nhood/L-R/co-occurrence analyses:
  - QTimer objects were garbage-collected immediately after `worker.start()` because no
    Python reference was retained; fixed by storing timers in `state['_spinner_timer']`
    and `state['_progress_timer']`
  - `qt_tqdm_context` was not patching `tqdm.std.tqdm` (only `tqdm.tqdm` and
    `tqdm.auto.tqdm`); squidpy's `parallelize()` imports from `tqdm.std` when
    `ipywidgets` is absent, so the monkey-patch had no effect; fixed by also patching
    `tqdm.std.tqdm`
  - `attach_tqdm_progress` now returns `(post_fn, timer)` instead of just `post_fn`

### Added
- **Progress feedback for long-running analyses** — status bar now updates during analyses
  instead of showing a static "Running…" message:
  - **Leiden clustering**: animated spinner (`| / - \`) with stage messages
    (`"Preparing data…"` → `"Computing neighbors…"` → `"Running Leiden algorithm…"`)
  - **Rank genes**: animated spinner while scanpy computes
  - **Neighborhood enrichment**: live permutation counter (`"Enrichment permutations: 450/1000 (45%)"`)
  - **Ligand-receptor**: live permutation counter (`"L-R permutations: 230/1000 (23%)"`)
  - **Co-occurrence**: live interval counter (`"Co-occurrence: 23/50 iterations"`)
- New helpers in `tabs/_helpers.py`: `attach_spinner`, `attach_tqdm_progress`,
  `qt_tqdm_context`, `ProgressMailbox`

### Changed
- **Cursor hover debounce** — `_on_cursor_move` in `02_xenium_viewer.py` now fires
  `get_value()` at most every 80 ms via a `QTimer` debounce, eliminating main-thread
  stutter during rapid pan/zoom.
- **Cluster coloring off main thread** — cluster filter (numpy masking) and
  `build_direct_label_colormap()` (100K dict allocations) now run inside the
  `@thread_worker`; the main-thread callback only sets `colormap` and calls
  `refresh()`, eliminating the UI stall on Apply.
- **Faster colormap dict** — `build_direct_label_colormap()` uses bulk `.tolist()`
  (C-level) before building the `color_dict`, eliminating 100K+ Python object
  allocations per colormap rebuild.
- **Revert vectorized cluster loop** — `get_cluster_colors()` reverted to a simple
  Python `for` loop; result is cached so the vectorization complexity wasn't
  visible to users.

### Fixed
- **"Filter by selected clusters" in Transcript Density** — the checkbox was
  silently ignored when `label_to_cluster` was not yet populated (e.g. before
  clicking "Apply Cell Coloring" or after switching to gene coloring). The fix
  pre-computes filter data on the main thread in `on_compute_density()`: if a
  clustering is selected it builds `label_to_cluster` on-the-fly when needed;
  if no clustering is available it shows "No clustering applied — filter
  skipped" in the status label and returns early instead of silently falling
  through.

### Added
- **Transcript cache completeness check** (`00_preprocess_transcripts.py`) — writes a
  `.complete` sentinel file to `transcript_cache/` on successful completion, storing the
  parquet mtime and size. Re-running the script detects an up-to-date cache and exits
  immediately. If the parquet has changed, stale feather files are wiped before
  reprocessing, preventing silent data duplication.

- **Interactive minimap** (`utils/minimap_widget.py`) — floating overlay in the
  top-right corner of the napari canvas showing the DAPI thumbnail with a white
  viewport rectangle. Updates live as the camera pans/zooms; clicking navigates
  the camera. Repositions itself when the window is resized. Instantiated
  automatically after dataset load in `_do_full_init()`.

- **Transcript density bins** — new "Transcript Density" section in the
  Transcripts tab. Aggregates per-gene transcript counts into spatial bins
  (user-controlled µm bin size via slider) and displays a `transcript_density`
  napari Image layer (hot colormap). Supports optional cluster-based filtering
  using centroids from the Cell Coloring tab. Computation runs in a background
  thread worker; contrast limits are set to the 99th percentile of non-zero bins.
  New `transcript_bins_layer` field added to `ViewerContext`.

- **Transcript density normalisation** — "Normalise by cells per bin" checkbox
  in the Transcript Density section. When ticked, each bin value is divided by
  the number of cells in that bin (transcripts-per-cell), correcting for local
  cell-density variation. Bins with zero cells remain zero. When combined with
  "Filter by selected clusters", the cell count per bin uses only selected-cluster
  cell centroids.


- **Empty-viewer startup** — cancelling the startup folder dialog (or omitting the
  CLI path argument) now opens a bare napari window instead of calling `sys.exit(0)`.
  A minimal File menu (Open Dataset... `Ctrl+O` + Preprocess Dataset...) is attached
  directly to the bare viewer via the new `_attach_bare_file_menu()` helper. When the
  user picks a dataset from File → Open Dataset, `_do_full_init()` builds the full UI
  (layers, dock widget, session restore) and replaces the bare menu. Works identically
  for the first load and subsequent dataset switches.

- **`_do_full_init()` helper** — extracted from `main()`: loads a dataset, creates
  `UMAPViewer`, populates napari layers, builds `ViewerContext`, rebuilds the control
  panel dock widget, and restores the session. Called both at initial startup (when
  `data_path` is given) and from `_on_open_dataset()` for every subsequent switch.
  `main()` is now ~100 lines shorter and `_on_open_dataset()` no longer duplicates the
  load/build sequence.

### Changed
- **Zarr skip logic in `PreprocessWorker`** — `PreprocessWorker.run()` now performs a
  staleness pre-check before calling `load_sdata()`: if `sdata_cached.zarr` already
  exists and its mtime ≥ `experiment.xenium` mtime, the zarr step is skipped and the
  progress dialog shows "Zarr cache already up to date, skipping: …" instead of the
  misleading "Creating zarr cache…" message. Uses the same mtime comparison as
  `load_sdata()`.

- **`_on_open_dataset()` refactored** — now uses `nonlocal ctx` and delegates the
  load/build sequence to `_do_full_init()`. Cleanup path (session save, UMAP close,
  generation counter increment, ref nulling, gc) is guarded by `if ctx is not None`
  so it also handles the first load from an empty viewer correctly.

- **`_on_preprocess_dataset()` hardened** — all `ctx` accesses are guarded by
  `if ctx is not None`, so preprocessing works from an empty viewer. The confirmation
  dialog's "session will be cleared" line is omitted when no dataset is loaded.

- **`_parse_args()`** — dialog cancel now returns `(None, no_cache)` instead of
  calling `sys.exit(0)`. `experiment.xenium` validation only runs when a path was
  actually selected.

### Added (continued from 2026-03-16)
- **File > Preprocess Dataset...** — new menu item that runs both preprocessing
  steps (zarr cache creation via `01_load_sdata.py` and per-gene transcript
  feather splitting via `00_preprocess_transcripts.py`) in a background thread
  so the UI stays responsive. Auto-detects whether the selected folder is a
  single Xenium dataset or a parent folder containing several datasets, and
  processes all found datasets sequentially. A modal progress dialog with a
  live progress bar is shown during preprocessing. The current viewer data is
  cleared before starting to avoid a RAM peak of old dataset + I/O overhead.
  Implemented via `PreprocessWorker(QThread)`, `_find_xenium_datasets()`, and
  `_make_progress_dialog()` in `02_xenium_viewer.py`; `create_file_menu()` in
  `tabs/_helpers.py` updated to accept the new `on_preprocess_dataset` callback.

### Changed
- **Memory: free old dataset before loading new one** — in `_on_open_dataset()`,
  all heavy `ctx` fields (`sdata`, `adata`, `clusterings`, `color_manager`,
  `transcript_loader`, layer references, etc.) are now explicitly set to `None`
  and `gc.collect()` is called after clearing napari layers (step 8b) and before
  `_load_dataset()` is called (step 9). This prevents peak RSS from reaching
  old-dataset + new-dataset simultaneously during a dataset switch.

## [Unreleased] — 2026-03-14

### Added
- **File > Open Dataset** — new menu entry (and `Ctrl+O` shortcut) lets the user
  switch to a different Xenium output directory without restarting the application.
  The current session is saved automatically before the swap; the new dataset's
  session is restored if a zarr cache exists. All napari layers are replaced and
  the full control-panel dock widget is rebuilt in-place.

  Implementation details:
  - `utils/viewer_context.py` — added `dataset_generation: int = 0` field.
  - `tabs/_helpers.py` — new `create_file_menu(ctx, on_open_dataset)` function;
    `_build_control_panel()` now accepts an optional `on_open_dataset` callback.
  - `02_xenium_viewer.py` — extracted `_load_dataset()`, `_populate_viewer()`,
    `_snapshot_layers()`, `_make_initial_state()`, `_make_initial_he_state()`, and
    `_make_initial_arms_state()` as standalone helpers; introduced `_app` container
    in `main()`; implemented `_on_open_dataset()` closure with re-entry guard,
    session save/restore, and full viewer reload sequence.
  - Stale thread-worker guard added to six tab modules (`tab_transcripts.py`,
    `tab_cell_coloring.py`, `tab_clustering.py`, `tab_roi.py`,
    `tab_he_registration.py`, `tab_arms.py`): each worker captures
    `ctx.dataset_generation` at dispatch time and silently drops its result if
    the generation has changed when the worker returns.

## [Unreleased] — 2026-03-10

### Fixed
- **Co-occurrence plot colors now match labels layer** — replaced seaborn `tab20`
  fallback in `make_co_occurrence_plot()` with `CLUSTER_PALETTE` from
  `utils/coloring.py`, ensuring line colors match the cell labels layer colors.

### Changed
- **Custom nhood enrichment heatmap now matches squidpy native style** — rewrote
  `make_nhood_enrichment_plot()` in `spatial_analysis.py` to use `imshow` +
  `make_axes_locatable` instead of `sns.heatmap`. Uses `viridis` colormap with
  symmetric normalization for z-score, colored cluster category bars on left + top
  edges (using `CLUSTER_PALETTE`), cluster labels on the left bar, right-side
  vertical colorbar with `%0.2f` ticks, no gridlines. Added `annotate` parameter
  (default False) with luminance-based white/black text. Figsize scales with cluster
  count matching squidpy's `(2*n//3, 2*n//3)` formula.

- **Decomposed `02_xenium_viewer.py` into per-tab modules** — the 4295-line monolith
  has been split into 14 new files. `_build_control_panel()` shrinks from ~3850 lines
  to ~80 lines (thin orchestrator). Each of the 11 tabs is now a standalone module in
  `scripts/tabs/` with a single `build_tab(ctx)` entry point. Cross-tab state is managed
  through a `ViewerContext` dataclass (`scripts/utils/viewer_context.py`) instead of 18
  closure-captured parameters. Shared helpers (status proxy, code recording, clustering
  refresh, label editor, plot save, etc.) live in `scripts/tabs/_helpers.py`. No
  behavioral changes — all features, session persistence, and code recording work
  identically.

  New files:
  - `scripts/utils/viewer_context.py` — ViewerContext dataclass
  - `scripts/tabs/__init__.py` — re-exports all tab builders
  - `scripts/tabs/_helpers.py` — make_tab, StatusProxy, create_shared_helpers, create_preferences_menu
  - `scripts/tabs/tab_clustering.py` — Tab 0: Leiden, import/export
  - `scripts/tabs/tab_cell_coloring.py` — Tab 1: gene/cluster coloring, filter checkboxes
  - `scripts/tabs/tab_transcripts.py` — Tab 2: multi-gene transcript overlay
  - `scripts/tabs/tab_umap.py` — Tab 3: UMAP viewer
  - `scripts/tabs/tab_roi.py` — Tab 4: ROI analysis + ROI DEG
  - `scripts/tabs/tab_he_registration.py` — Tab 5: H&E registration
  - `scripts/tabs/tab_gene_analysis.py` — Tab 6: rank genes, dotplot, volcanos
  - `scripts/tabs/tab_ligrec.py` — Tab 7: ligand-receptor
  - `scripts/tabs/tab_nhood.py` — Tab 8: neighborhood enrichment
  - `scripts/tabs/tab_co_occurrence.py` — Tab 9: co-occurrence
  - `scripts/tabs/tab_arms.py` — Tab 10: ARMS overlay

### Added
- **ARMS tab code recording** — All ARMS operations (H&E load, landmark registration,
  GeoJSON/CSV load, tile DEG, DEG export, volcano plots) now emit `_record_code()`
  entries so they appear in the reproducible `code.py` journal.

## [2026-03-09]

### Fixed
- **ARMS metadata not persisting across sessions** — `save_session()` overwrote
  real-time-saved ARMS attrs (filename, affine, shape, paths) with empty snapshot
  values, causing ARMS registration to be lost on reload. Fixed by: (1) preserving
  existing ARMS attrs before the overwrite and using them as fallback, (2) saving
  all essential ARMS metadata in `_save_arms_affine_to_sdata()` (not just
  affine/flips), (3) calling `_save_arms_affine_to_sdata()` after ARMS restore to
  repair any previously corrupted attrs.
- **NameError on close** — `_on_viewer_closing` referenced `_arms_state` which is
  defined later in `main()`; wrapped in `try/except NameError` so closing the
  viewer before the ARMS tab initializes no longer crashes.

### Added
- **ARMS Tile DEG analysis** — "Run ARMS Tile DEG" button in ARMS Overlay tab
  performs differential expression between cells grouped by ARMS tile cluster ID.
  Tile polygons are transformed from GeoJSON space to Xenium pixel space via the
  registration affine before point-in-polygon tests. Clusters with <10 cells are
  excluded. Results display in a text widget and can be exported to CSV. Session
  persistence saves/restores DEG results across viewer restarts.
- **ARMS pairwise volcano plots** — "Generate ARMS Volcano Plots..." button
  (enabled after ARMS Tile DEG completes) runs pairwise DEG for every pair of
  ARMS tile clusters and saves volcano PNGs to a user-selected directory.
  Progress updates shown in status bar. `compute_arms_tile_deg()` now returns
  the normalized subset adata as a 3rd element for downstream pairwise analysis.
- **Pairwise volcano plots** — "Generate All Volcano Plots..." button in Gene
  Analysis tab runs DEG for every pairwise cluster comparison and saves volcano
  PNGs (3-color scatter with threshold lines and top gene labels) to a
  user-selected directory. Progress updates shown in status bar.
- **ARMS Overlay tab (Tab 10)** — load a larger ARMS H&E image, align it to
  Xenium coordinates via manual landmark registration (same workflow as the H&E
  Registration tab), then load GeoJSON tile boundaries + CSV cluster IDs to
  display 288 tile polygons colored by genomic cluster (1–8). The affine
  transform is applied to both the H&E image and the polygon shapes layer, so
  re-registering landmarks automatically re-aligns everything.
- **Auto-save reproducible code on exit** — code recording is now on by default
  and automatically saves to `{data_dir}/code.py` when the viewer closes (if
  recording is enabled and journal is non-empty).
- **Persist custom clusterings across sessions** — Leiden and imported clusterings
  are saved as parquet files in `viewer_session/clusterings/` inside the zarr
  cache, and restored on next launch.
- **ARMS overlay session persistence** — ARMS H&E image, registration affine,
  flip state, landmarks, and GeoJSON/CSV file paths are saved to the zarr cache
  and restored on next launch. The H&E image is stored as
  `sdata.images["arms_he_image"]` with a multiscale pyramid; the affine and
  landmarks are saved in `viewer_session/arms/`. GeoJSON/CSV tiles are re-parsed
  from the original file paths on restore (with a warning if files have moved).

### Changed
- `record_code` default changed from `False` to `True` (always-on recording).

## [Unreleased] — 2026-03-06

### Added
- **Import clustering from CSV/TSV** — "Import Clustering..." button in the
  Clustering tab. Reads a file with `cell_id` + `group` columns (auto-detects
  tab vs comma separator). Supports string-valued group names (e.g. imported
  cell type annotations). Imported clusterings appear in all downstream dropdowns.
- **Export clustering to CSV/TSV** — "Export Clustering..." button exports the
  currently selected clustering with cell_id and group columns. Custom cluster
  labels are applied to the exported group names when available.
- **Per-clustering label storage** — cluster labels are now stored per-clustering
  (nested dict `{clustering_name: {cluster_id: label}}`), so labels don't collide
  across different clusterings. Session save/restore updated with backward
  compatibility for the old flat format.
- **Multi-column label editor** — the "Edit Cluster Labels..." dialog now uses a
  multi-column grid layout (up to 3 columns, ~10 rows per column) with a scroll
  area, preventing the dialog from being too tall for clusterings with many groups.
  Handles both integer and string cluster IDs.
- **Label propagation to all plots** — cluster labels now appear in:
  - Neighborhood enrichment heatmap (axis tick labels)
  - Co-occurrence plots (subplot titles + legend entries)
  - Ligand-receptor dotplot (column axis labels)
  - Rank genes panel plot (group titles)
  - Mouse hover status bar (shows label instead of raw cluster ID)
  - Cluster filter checkboxes (show labels when available)
  - Dotplot (already worked, now uses per-clustering lookup)

### Changed
- **Cluster label editor shared between tabs** — Gene Analysis and Clustering
  tabs now share the same `_build_label_editor_dialog()` helper, eliminating
  code duplication. Both editors use the per-clustering label storage.
- **String cluster ID support** — `add_clustering_to_obs()` in gene_analysis.py
  and `_get_cluster_ids_per_obs()` in the viewer now handle string-valued
  cluster IDs (from imported clusterings) via factorization fallback.

## [Previous] — 2026-03-05

### Added
- **Cluster label editor in Clustering tab** — "Edit Cluster Labels..." button
  in the Clustering tab opens a dialog to rename clusters (manual cell type
  annotation) using the Cell Coloring tab's selected clustering.

### Changed
- **Co-occurrence plot colors match napari** — line colors in co-occurrence plots
  now use the same palette as napari cell coloring (CLUSTER_PALETTE) instead of
  seaborn tab20. Falls back to tab20 for clusters without a stored color.
- **Clustering sync across tabs** — selecting a clustering in the Cell Coloring
  tab now auto-sets the same clustering in Gene Analysis, Ligand-Receptor,
  Nhood Enrichment, and Co-occurrence tabs (one-directional sync).

## [Previous] — 2026-03-04

### Added
- **Leiden clustering tab** — new "Clustering" tab (Tab 0) with configurable
  `n_neighbors` (5–50), `n_pcs` (10–50), and `resolution` (0.1–5.0) parameters.
  Runs `sc.pp.neighbors` + `sc.tl.leiden` on a worker thread, stores results as
  `leiden_r{resolution}` in the clusterings dict, and refreshes all downstream
  ComboBoxes (Cell Coloring, Gene Analysis, L-R, Nhood, Co-occurrence).
  Reproducible code recording supported.

## [Previous] — 2026-03-03

### Added
- **Interaction database filtering for L-R analysis** — 4 checkboxes (OmniPath,
  LigRecExtra, PathwayExtra, KinaseExtra) to select which OmniPath interaction
  datasets are queried, plus a "CellPhoneDB only" toggle to restrict to canonical
  ligand-receptor pairs. Selections are passed via `interactions_params` to
  `sq.gr.ligrec()`. Unchecking PathwayExtra + KinaseExtra removes intracellular
  proteins (e.g. TP53) that leak through as false "ligands".

- **Plot format preference (PNG / SVG)** — new Preferences → Plot format menu in the
  napari menu bar with exclusive PNG/SVG radio actions. All 4 plot save handlers
  (dotplot, L-R, nhood enrichment, co-occurrence) now respect the chosen format.
  Default is PNG. `matplotlib.savefig()` infers format from file extension.

- **Native squidpy heatmap** for Neighborhood Enrichment — `on_show_nhood_plot()`
  now uses `sq.pl.nhood_enrichment()` when no cluster filter is active, falling
  back to the custom `make_nhood_enrichment_plot()` when clusters are filtered or
  on session restore. `_adata_norm` and `_cluster_key` are now stored on the result
  dict (matching co-occurrence).

- **"Filter targets" checkbox** in Co-occurrence tab — when checked with a cluster
  filter active, restricts target lines (not just query subplots) to the selected
  clusters. Three-path plot logic: (1) filter targets ON → custom plot with both
  query and target restricted, (2) filter targets OFF → squidpy native plot with
  cluster-filtered subplots, (3) session restore fallback → custom plot.
  `make_co_occurrence_plot()` gains a `target_clusters` parameter.

- **Co-occurrence tab** (Tab 9) — squidpy-based spatial co-occurrence analysis
  (`sq.gr.co_occurrence`). Computes how cluster types co-occur spatially across
  increasing distance radii. Configurable: clustering, distance bins (10-100).
  Displays summary with cluster count and distance range. Line-plot visualization
  showing co-occurrence probability vs distance for selected clusters — one subplot
  per query cluster with colored lines for each target cluster and baseline at y=1.
  Uses Cell Coloring tab's cluster filter to select which clusters get subplots.
  Export long-form CSV, save plot as PNG. Session persistence via zarr. New functions
  in `spatial_analysis.py`: `run_co_occurrence()`, `make_co_occurrence_plot()`.

- **Neighborhood Enrichment tab** (Tab 8) — squidpy-based neighborhood enrichment
  analysis (`sq.gr.nhood_enrichment`). Computes which cluster types are spatially
  enriched or depleted in each other's neighborhoods via permutation testing.
  Configurable: clustering, permutation count (100-1000), neighbor count (3-20).
  Displays summary with top enriched/depleted cluster pairs. Heatmap visualization
  (z-score or count mode) via seaborn with diverging colormap, supports filtering
  by the Cell Coloring tab's cluster selection. Export z-score matrix as CSV, save
  plot as PNG. Session persistence via zarr. New functions in `spatial_analysis.py`:
  `run_nhood_enrichment()`, `make_nhood_enrichment_plot()`.

### Changed
- Renamed `_get_lr_cluster_filter()` to `_get_cluster_filter()` — shared helper for
  cluster filtering in both L-R plot and nhood enrichment plot.

## [Unreleased] — 2026-03-02

### Added
- **Session persistence** — viewer state (ROIs, H&E registration, analysis results,
  cluster labels) is automatically saved to `sdata_cached.zarr/` when the viewer closes.
  The H&E image is stored as a spatialdata multiscale image element (`images/he_image`)
  with its affine transformation, so the next launch restores the H&E overlay with
  registration already applied — no need to re-load or re-register. ROI polygons,
  cluster labels, and analysis results (rank genes, ROI DEG, L-R) are persisted in
  `viewer_session/` and restored on startup. Affine transformations are saved in
  real-time (on each registration/flip change). Skipped when `--no-cache` is used.
  New module `scripts/utils/session.py` with `save_session()` and `load_session()`.

### Previously added
- **Gene Analysis tab** (Tab 6) — rank marker genes per cluster using Wilcoxon, t-test, or
  logreg methods via scanpy's `rank_genes_groups`. Features: configurable top-N genes,
  dotplot visualization with optional dendrogram, editable cluster labels (dialog to rename
  cluster IDs for publication-ready plots), rank genes panel plot, results preview in
  monospace text area, and full CSV export. New module `scripts/utils/gene_analysis.py` with
  `get_normalized_adata()` (cached log-normalization + PCA), `add_clustering_to_obs()`,
  `run_rank_genes()`, `make_rank_genes_dotplot()`, `make_rank_genes_plot()`.
- **ROI differential expression** — extended ROI Analysis tab (Tab 4) with genome-wide
  DEG between drawn ROI polygons. Draw >= 2 ROIs, select method (Wilcoxon/t-test), optional
  cluster filtering, then compute DE across all genes. Uses `compute_roi_deg()` which
  normalizes the subset, runs `rank_genes_groups` with `reference='rest'`, and returns
  full results table. Export to CSV supported.
- **Ligand-Receptor tab** (Tab 7) — squidpy-based spatial ligand-receptor interaction
  analysis. Computes spatial neighbor graph, runs permutation-based L-R testing
  (`sq.gr.ligrec`), displays summary of significant interactions, and generates
  interaction dotplot. Configurable: clustering, permutation count (100-1000), neighbor
  count (3-20), p-value threshold. Export means and p-values as CSV, save plot as PNG.
  Warns when the 480-gene Xenium panel lacks common L-R pairs. New module
  `scripts/utils/spatial_analysis.py` with `compute_spatial_neighbors()`, `run_ligrec()`,
  `make_ligrec_plot()`.

## [Unreleased] — 2026-03-01

### Added
- **H&E image registration** — new "H&E Registration" tab (5th tab). Load an H&E
  OME-TIFF/SVS image as a multiscale overlay, place paired landmark points on both
  Xenium and H&E images, then compute a similarity transform (rotation + uniform scale +
  translation) to align the H&E onto Xenium coordinates. Features: opacity slider,
  per-landmark residual display, save/load landmarks as JSON, save affine matrix.
  Flip checkboxes (vertical/horizontal) let you mirror the H&E before placing
  landmarks, useful when the H&E was scanned in a different orientation.
  New module `scripts/utils/registration.py` with `load_he_pyramid()`,
  `compute_landmark_affine()`, `save_landmarks()`, `load_landmarks()`.
- **Coarse tissue-outline alignment** — "Coarse Align" button in H&E Registration tab
  automatically snaps H&E roughly into position before manual landmark placement.
  Extracts tissue masks from morphology_focus (max-projection across channels) and H&E
  (HSV saturation channel), then computes a similarity transform using OpenCV image
  moments with multi-rotation hypothesis testing (0/90/180/270 deg, best IoU wins).
  The coarse affine is a visual stepping-stone; computing landmark registration replaces
  it entirely. New functions in `registration.py`: `extract_tissue_mask()`,
  `extract_tissue_mask_fluorescence()`, `extract_tissue_mask_he()`,
  `compute_coarse_affine()`.
- **ROI polygon analysis** — new "ROI Analysis" tab (4th tab). Draw polygons on the
  napari Shapes layer ("ROIs"), click "Calculate Expression" to get per-region stats
  (cell count, mean, median, std, min, max) for the selected gene. "Export CSV" saves
  per-cell data (region_id, cell_id, x/y centroid in microns, expression). Uses shapely
  `contains_xy` for fast point-in-polygon queries on precomputed cell centroids.
  Includes pairwise Welch's t-tests between regions with Benjamini-Hochberg
  correction when >1 comparison.

### Changed
- **Multi-cluster filter** — "Filter by cluster" now supports selecting multiple clusters
  simultaneously via checkboxes (3-column grid with Select All / Deselect All). Works
  in both Gene Expression and Cluster coloring modes. Replaces the single-cluster
  ComboBox dropdown.
- **White background toggle** — new checkbox at the top of the Cell Coloring tab sets the
  napari canvas background to white, useful for H&E overlay visibility or light-colored
  cluster inspection.
- **Cluster filter works in both modes** — "Filter by cluster" checkbox and cluster ID
  selector are now available in Cluster coloring mode (not just Gene Expression). When
  enabled, only cells belonging to the selected clusters are shown; all others are
  transparent.
- **UMAP window no longer auto-opens** — `color_by_gene()` and `color_by_cluster()`
  store colors/metadata but only update the viewer if already open. Click "Show UMAP
  Window" to open it manually; colors and title are applied on first open.
- **Min/max colormap scaling** — gene expression coloring now normalizes non-zero cells
  using vmin/vmax of expressing cells (instead of 0-to-max). This spreads the full
  viridis range across expressing cells so low-expressors are visually distinguishable.
  Zero-expression cells remain transparent.

### Previously added
- **Multi-gene transcript overlay** — add up to 10 genes simultaneously with distinct
  colors (yellow, cyan, magenta, orange, green, sky blue, red, violet, pink, brown).
  New gene list widget with Add/Remove/Clear All buttons and color legend.
  `TranscriptLoader.get_multi_gene_points()` merges per-gene points with palette colors.
- **Cluster filter for gene expression** — checkbox "Filter by cluster" in Cell Coloring
  section. When enabled, only cells in the selected cluster are colored; rest are
  transparent. `CellColorManager.get_gene_colors_filtered()` method.
- **Generic dataset loading** — CLI argument `python 02_xenium_viewer.py /path/to/data`
  or Qt file dialog if no argument given. Reads `pixel_size` from `experiment.xenium`.
  All scripts accept `--data-dir` / positional argument. Transcript cache stored alongside
  data (`data_dir/transcript_cache/`) instead of inside `scripts/`.
- **Zarr cache** — `01_load_sdata.py` writes `sdata_cached.zarr` on first load;
  subsequent launches use the cache (~60-70% faster). Staleness detection via
  `experiment.xenium` mtime. `--no-cache` flag to force rebuild.
- **Deferred UMAP scatter** — UMAP widget shows placeholder text until first
  "Apply Cell Coloring" click, avoiding 318K-point scatter at startup.
- **Timing instrumentation** — wall-clock times printed for each startup phase.

### Changed
- `load_sdata()` now skips `cells_boundaries` and `nucleus_boundaries` (were hidden
  anyway; 318K polygon shapes freeze napari). Saves ~30-40% of xenium() load time.
- `load_umap()` and `load_clusterings()` accept a `path` parameter instead of using
  module-level globals.
- Transcript cache default location changed from `scripts/transcript_cache/` to
  `data_dir/transcript_cache/`.

## [0.1.0] — 2026-02-28

### Added
- `scripts/00_preprocess_transcripts.py` — one-time script to split `transcripts.parquet`
  (1.3 GB, 128M rows) into per-gene feather files in `scripts/transcript_cache/`.
  Enables <100 ms per-gene transcript loading in the viewer (vs 4–5 s from parquet).
- `scripts/01_load_sdata.py` — SpatialData loader using `spatialdata_io.xenium()` with
  5-level software image pyramid for the morphology_focus TIFFs (no internal pyramid).
  Also loads UMAP coordinates and all cluster assignments from `analysis/`.
- `scripts/utils/coloring.py` — `CellColorManager`: builds napari `DirectLabelColormap`
  from gene expression (viridis/magma/plasma/RdBu/YlOrRd) or cluster assignments.
  Transparent alpha for zero-expression cells. LRU cache on repeated gene queries.
- `scripts/utils/transcript_index.py` — `TranscriptLoader`: loads per-gene transcript
  x/y coordinates from feather cache with parquet fallback.
- `scripts/utils/umap_widget.py` — `UMAPWindow`: linked matplotlib UMAP scatter
  (318K cells, rasterized=True for performance). Syncs colors with napari Labels layer.
  Handles 91-cell UMAP/adata mismatch via `.reindex()`.
- `scripts/02_xenium_viewer.py` — main entry point. Opens napari with:
  - 4-channel morphology_focus image (DAPI, ATP1A1/CD45/E-Cad, 18S, AlphaSMA/Vim)
  - Cell and nucleus labels layers (raster masks for fast coloring via DirectLabelColormap)
  - Transcript points layer (populated on demand from feather cache)
  - Magicgui control panel: gene/cluster coloring, colormap, transcript toggle, QV slider
  - Linked matplotlib UMAP window
- `.gitignore` — excludes parquet, ome.tif, zarr.zip, h5, and transcript_cache/
