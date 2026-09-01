# Changelog

Dates are ISO. Versions follow [Semantic Versioning](https://semver.org). New work goes
under an `## [Unreleased]` heading above the newest release; the dated per-session
entries under **Development log** are the closed pre-1.0.0 record.

## [Unreleased]

### Added
- **PALMS has a DOI.** Zenodo archived the `v1.0.1` release and minted
  **10.5281/zenodo.22218654** (concept, always the newest version) alongside
  10.5281/zenodo.22218655 for that release specifically. `CITATION.cff` carries the
  concept DOI, which is what GitHub's "Cite this repository" button reads, and the README
  gains a badge and says which of the two to cite — a reader following the concept DOI
  reaches whatever superseded it, which is right for citing the software and wrong for
  citing a result. `docs/Home.md` carries the same "Citing PALMS" section, since the wiki
  and the Read the Docs site are the front page for everyone who never opens the repo, and
  neither renders `CITATION.cff` or the README's badge.

## [1.0.1] — 2026-09-01

### Added
- **`tests/test_version_consistency.py`.** The version is stated in four places —
  `pyproject.toml`, `palms.__version__`, `CITATION.cff` and the newest `CHANGELOG.md`
  heading — and nothing derives any of them from another. A release that bumps three of
  the four ships a package whose metadata disagrees with itself and a DOI whose version
  field is wrong, which a minted DOI makes permanent. Found while cutting this release,
  with all four in hand.

- **A custom segmentation now reaches the recorded notebook.** Tools → Segmentation
  swaps the Xenium cells for a custom set, rebinding `ctx.adata` and clearing every
  derived result — and it recorded nothing at all. The preamble went on saying
  `adata = sdata["table"].copy()`, so every node recorded after the swap claimed to be
  about the Xenium cells. That is the worst shape a provenance defect can take: replaying
  such a notebook does not fail. It runs the whole analysis against a different set of
  cells and reports success, so `allow_errors=False` — and `verify_notebook.py` with it —
  sees nothing wrong. The numbers come out; they are about something else.

  The swap is recorded by **upserting `preamble`**, not by adding a node of its own,
  because the preamble is what says where `adata` comes from. Two things follow, both
  wanted: the notebook's load cell reads `tables["custom_table"]`, and every result
  already in the graph is flagged stale — they were computed on cells that are no longer
  bound, which the GUI already acts on by dropping them from its own state. The custom
  table is not in the 10x output (it comes from `extract_seurat_segmentation.R` +
  `build_custom_segmentation.py`), so the recorded cell reads it from the store the tab
  caches it into — directly on a crop export, whose `sdata` *is* that store, and via one
  `sd.read_zarr` on a raw dataset.

  `app.py` seeds `segmentation_source` from the session **before** the launch re-emit of
  the preamble. Letting the tab's own restore handler correct it afterwards would upsert
  a changed preamble on every launch and mark the whole notebook ⚠ for nothing — the same
  defect a manual dataset rename used to cause. `tests/test_segmentation_recording.py`
  covers both directions, both store layouts, the staleness, and the launch round trip,
  and guards at the source that neither swap path stops recording.

- **`scripts/prepare_demo_dataset.sh`** — stages a publishable copy of a dataset for
  `capture_screenshots.py`. The screenshots are published and two panels print the
  dataset's absolute path into the widget, so the path in the copy is the path on the
  wiki. Three steps, and the third is the one that gets forgotten: copy (dropping the
  legacy `xenium_viewer.log` the Cache tab lists by name), repoint the absolute paths
  recorded in the provenance graph, and replace the ARMS scan's filename — which reaches
  the napari layer list and is therefore in *every* full-window screenshot, and which is
  a slide identifier. It refuses to remove anything at the destination that is not
  already a Xenium dataset, and refuses a destination inside the source in either
  direction — that second case would `rm -rf` a subdirectory of the source, and `cp`
  only catches it afterwards. The replacement step goes through `zarr_safe`'s
  `safe_group_update` like every other write to a store, and deliberately does not echo
  the previous filename: a pasted log is a published log.

- **`docs/user-authored-analyses.md` is tracked.** A 751-line design note, written
  2026-08-26 and untracked since: can a user add a *new* analysis and a GUI element for
  it, with the app linting it and generating its tests. The answer is yes and most of the
  machinery exists — but not the part it looks like: `~/.config/palms/templates/` is a
  *shadowing* mechanism, not a plugin one, so a template id that this package did not
  ship cannot exist. Untracked, that finding lived in one working copy.

  It joins the design notes already in `docs/` (`user-configurable-templates-todo.md`,
  `reproducible_notebook_plan.md`, the two migration notes). Its lower-case filename is
  what keeps it out of the published site and the wiki — `mkdocs.yml`'s `exclude_docs`
  and `scripts/push_to_wiki.sh` both take Title-Case as "page" and lower-case as
  "internal note".

- **`.github/workflows/release.yml`** — a `v*` tag builds the sdist and wheel and
  publishes them with PyPI **Trusted Publishing** (OIDC), so no API token is stored in
  the repo or anywhere else. Build and publish are separate jobs, and only the publish
  job holds `id-token: write`, so build code from the tree cannot reach the credential.
  Two checks run before the upload, because both failures otherwise surface only after a
  tag is public: the metadata must contain no direct references, and `pyproject`'s
  version must equal the tag.

### Changed
- **The citation carries an institutional address.** `CITATION.cff` named a personal
  gmail account; Zenodo reads that file for the DOI record, and a DOI is permanent.

- **The repository moved to `rao-genomics-lab/palms`.** GitHub redirects the old
  `sraorao/palms` URLs, so nothing was broken — but a redirect is not what a citation
  should carry, and `CITATION.cff` is what Zenodo reads when it mints the DOI. All 17
  in-repo references now name the new owner: `CITATION.cff`, `pyproject.toml`'s three
  `[project.urls]`, `mkdocs.yml`, the README badge and clone lines, `docs/Installation.md`,
  and `scripts/push_to_wiki.sh` (the wiki moved with the repo — verified by resolving
  `palms.wiki.git` under the new owner, since that was an open question at the rename).
  The four `sraorao/insituCNV-copykat` links are deliberately unchanged: the fork did not
  move.

  **`release.yml` says what the move costs.** PyPI Trusted Publishing registers a
  *repository*, so a transfer invalidates it: the OIDC claim carries the new owner and
  PyPI rejects the upload at the end of the run, **after the tag is already public**. The
  publisher has to be re-added under the new name before the next tag. 1.0.0 published
  from the old one and is unaffected.

- **Two reference tabs are photographed with something in them.** Tools → Dataset was
  captured before **Scan Dataset** had been pressed, so the page documenting the
  inventory tree showed an empty box; Tools → Templates was captured with no template
  selected, so the contract, the default-vs-yours panes and the live preview were all
  blank. The capture now primes both. `--only <substring>` re-takes matching shots
  without redoing the other fifty-one, since a full run drives real analyses and takes
  half an hour.

- **The screenshot capture stages the canvas before it shoots, and shoots wider.** A
  dataset carrying a real session restores registration landmarks, an ARMS scan and patch
  overlays; they sit outside the Xenium extent, so the camera framed their union and left
  the tissue a few pixels across. The capture now hides everything outside an allow-list —
  an allow-list, because a session can restore arbitrarily named overlays and an
  unrecognised one must not end up in a published image — and sets the camera from
  `cell_labels`' own extent, since napari's `fit_to_view` measures
  `layers._extent_world_augmented` and ignores `visible`. The window is 1800×1000 rather
  than 1400×900: at the old size the Controls dock and the layer panel left the canvas
  about 400px across.

- **The CNV dependency comes from PyPI instead of a git URL, which is what makes PALMS
  itself publishable.** PyPI rejects any distribution whose metadata carries a direct
  URL reference, and `palms`'s `cnv` and `full` extras carried
  `insitucnv @ git+https://…@v0.2.0`. The upload fails on that line alone, after the tag
  exists. Upstream's own `insitucnv` on PyPI cannot stand in: 0.1.0 pins `anndata<0.12`
  and `pandas<3`, which is exactly what the fork relaxed and why it exists.

  So the fork is published under its own name, `insitucnv-copykat` — `insitucnv` on PyPI
  is upstream's and stays theirs — and `environment.yml`, `environment-copykat.yml` and
  the `cnv` extra all require `insitucnv-copykat>=0.3`. It still **imports** as
  `insitucnv`, so no code changed; installing it beside upstream's distribution would
  have the two fight over that import name.

  `[tool.hatch.metadata] allow-direct-references` is removed with the last direct
  reference: it existed only to let hatchling emit one, and without it a future direct
  reference fails at build time rather than at upload time.

### Fixed
- **Every "Edit on GitHub" link on the docs site 404'd.** mkdocs defaults `edit_uri` to
  `edit/master/docs/` and the default branch is `main`. Set explicitly, on the first site
  this repo has ever published.

- **The docs site had no page at its root, so Read the Docs refused to serve it.** The
  first RTD build after the repo went public failed with "Index file is not present in
  HTML output directory" — after `mkdocs build --strict` had passed locally and in CI on
  every push for weeks. `docs/` is GitHub Wiki source, where the landing page must be
  `Home.md`; mkdocs's is `index.md`, and mkdocs builds a rootless site without complaint.
  `--strict` does not change that, which is why only the publisher caught it.

  Renaming works in neither direction: `Home.md` is what `scripts/push_to_wiki.sh`
  publishes as the wiki home, and a lower-case `docs/index.md` would be dropped by
  `exclude_docs` (`[a-z]*.md`) — silently, since that pattern is how internal notes are
  kept off the site. So the mapping is a build-time hook, like the wiki link rewrite
  beside it, and the source keeps one convention. `mkdocs_hooks.on_files` retargets
  `Home.md` to the site root; the page *moves* rather than being duplicated, so there is
  one URL for it, the nav resolves, and `sitemap.xml` carries the root.

  CI gained the check `--strict` cannot make: `test -f /tmp/site/index.html` after the
  build. Verified against mkdocs 1.6.1 — before, `site/` had `Home/index.html` and no
  root page; after, no dangling `Home/` href survives anywhere in the built HTML.

- **The pre-publication audit's own entry reprinted the path it removed.** The
  2026-08-29 entry quoted the deleted `capture_screenshots.py` constant verbatim —
  absolute path, collaborator's name, slide ID — so making the repo public would have
  published exactly what the audit existed to prevent, in the record of its own removal.
  It now describes the leak the way the 2026-08-31 screenshot entry does, which is the
  form to keep: a changelog is a published file, and a fixed leak quoted in full is
  still a leak. Two references to a planning document outside this repo went with it;
  they named a phase a reader here cannot resolve.

- **The docs never told anyone to `pip install palms`.** PALMS 1.0.0 has been on PyPI
  since 2026-08-30 and the fresh-machine install was verified on Linux and macOS, but
  `README.md`, `Home.md` and `Installation.md` all opened with `git clone` +
  `./scripts/install.sh`, and the extras table gave only the checkout form
  (`pip install -e ".[cnv]"`). A reviewer who installs the published package was reading
  instructions for a different route. All three now lead with the wheel and keep the
  conda install as what it actually is: the way to develop PALMS, run the suite, or reach
  the CopyKAT backend. `Home.md` also listed conda/mamba as a *requirement*.

- **Preferences → Plot format was documented as the wrong menu with the wrong default.**
  `Interface-Overview.md` said "**PNG** or **SVG** … Defaults to SVG". The menu offers
  PNG + PDF, PNG, PDF and SVG, and defaults to **PNG + PDF** — a PNG to look at, a PDF to
  hand to a journal. Stale since the plots rework.

- **Two CNV controls were quoted in the docs with a spelling the app does not use** —
  "Neighbors (expression graph)" and "Smoothing neighbors" against the UI's
  "Neighbours". A reader searching the page for the label in front of them missed it.
  `README.md` also still named the CNV dependency `insitucnv`, which is upstream's
  distribution; PALMS requires `insitucnv-copykat`.

- **`push_to_wiki.sh` could publish but never unpublish.** The copy loop only ever wrote,
  so a page or screenshot deleted from `docs/` stayed on the wiki indefinitely — the five
  screenshots whose tutorial steps were dropped were still being served, and an internal
  note that leaked before the naming convention existed could never be taken down. It now
  prunes anything the wiki holds that `docs/` no longer has, or that is not a wiki page by
  the convention.

- **Every tutorial illustrated all of its steps with a photograph of its first one.**
  30 `tutorial-*.png` were 9 distinct images: all five H&E-registration steps were one
  picture, all six ARMS steps another, all five ROI steps a third. `TUTORIAL_SHOTS`
  named a tab per step and consecutive steps of a tutorial name the same tab, so with
  nothing acting on the app between grabs the capture wrote byte-identical files. It was
  invisible while every one of them was wallpaper.

  The capture now **drives the app into each step's state** before grabbing — it runs
  Leiden, ranks genes, colours by gene and by cluster, loads transcripts, opens the linked
  UMAP window, calculates ROI expression and ROI DEG, coarse-aligns the H&E, runs the ARMS
  tile DEG, draws annotation polygons and assigns types to them — through the real widgets,
  found by magicgui's `native._magic_widget` back-reference rather than by Qt layout
  position. All 52 images are now distinct, checked by pixel hash.

  Five files are gone rather than duplicated: the steps whose only action is a file dialog
  (Save Landmarks…, Export GeoJSON…, Save Volcano Plot…) cannot be driven, and their
  tutorial pages carry one fewer image instead of a repeat of the previous one.

- **Every full-window docs screenshot was a picture of the desktop wallpaper.**
  `scripts/capture_screenshots.py` grabbed the window with
  `QScreen.grabWindow(qt_window.winId())`. Under Qt6 that argument form is unsupported
  on this platform and returns a fragment of the *root* window, so `interface-overview`
  and all 30 `tutorial-*.png` — one on nearly every docs page and on the wiki — showed
  wallpaper and GL debris with no window, no Controls dock and no canvas. The 26
  `tab-*.png` were unaffected: they already went through `dock.widget().grab()`.
  Now `qt_window.grab()`, a widget render, which does not depend on a compositor, works
  over remote X, and excludes the WM title bar — so `viewer.title` can no longer leak a
  dataset folder name into a published image.

  A widget grab alone was still wrong, though: it reads the vispy canvas's last painted
  framebuffer, and the part of that framebuffer covering an overlay drawn before the
  camera moved is never repainted. Shots came back with a block of an earlier draw frozen
  in the corner, pixel-identical across frames while the rest of the canvas tracked the
  layers correctly; `repaint()` on the canvas widget does not clear it. The canvas is now
  rendered on demand through vispy (`viewer.screenshot(canvas_only=True)`) and painted
  over the canvas widget's rectangle in the grab.

- **The published screenshots showed a local path, a collaborator's name and a real
  slide ID.** `docs/screenshots/tab-dataset.png` and `tab-cache.png` printed a dataset
  directory into the panel, and that directory was a working dataset. The capture already
  took the dataset as an argument; the images are recaptured against one whose path is
  publishable. `tab-crop-dataset.png` also still read "opened directly with xenium-viewer"
  and predated the "Controls" dock rename — both were stale pictures of current, correct
  source.

- **A dead import wearing an invalid `# noqa`.** `tab_he_registration._restore_session`
  opened with `from palms.tabs._helpers import StatusProxy as _SP  # noqa: avoid
  circular`. `_SP` is never used, `StatusProxy` is already imported at module level and
  used there, and "avoid circular" is prose where ruff expects rule codes — so the
  directive was invalid, suppressed nothing, and printed a warning on every lint run
  including CI. Removed rather than corrected: there is nothing left to suppress.
- **A crop export rewrote recorded paths into files it does not copy.**
  `rewrite_graph_paths` repoints *every* absolute path in the carried provenance
  graph at the export. That is right for the preamble's `data_path` and wrong for
  anything under a directory the export does not carry. The case that bites is 10x's
  own clustering CSVs: a parent's `clustering:<key>` node legitimately reads
  `<parent>/analysis/clustering/<key>/clusters.csv`, the rewrite turns it into
  `<export>/analysis/…`, and an export carries `experiment.xenium`,
  `transcripts.parquet` and the zarr — never `analysis/`. Measured on
  `demo_data/crop_6`: with the preamble fixed, the replay got three cells in and died
  with `FileNotFoundError`, and all three `clustering:*` nodes were of that shape.
  A relaunch could not repair it, because `_record_clustering` returns early when the
  node exists — deliberately, so a loader never overwrites a producer's code.

  The export now checks each rewritten path against **what the staging directory
  actually holds**, rather than against a list of what an export is believed to hold,
  and substitutes the reload-from-the-store cell for a `clustering:<key>` node whose
  file is not there. Copying `analysis/` instead would be wrong: a crop is a cell
  subset, so the parent's CSV reindexes to NaN for every cell the crop does not
  contain. A dangling path with no such substitution is **reported in the export's
  own notes** rather than silently left, since the notebook will fail on it.

  **A reloaded clustering is not a recomputed one**, and the cell says so where the
  reader is. A crop export's notebook can reproduce the parent's labels but cannot
  re-derive them, which matters to anything quoting an ARI off a crop export.
  `utils/clustering_code.py` now holds that cell's text once, because
  `_record_clustering` and the export must not drift.

- **A Crop Dataset export declares itself cache-only, and nothing read the
  declaration.** `crop_export` stamps `cache_only: True` into the export's cache
  manifest, and `write_manifest`'s own docstring says it is there "for readers that
  would otherwise have to infer it from absent files". Five readers inferred it
  anyway — both `loader` call sites, the dataset inventory's description of the raw
  section, the Force Rebuild guard and the recorded preamble all recomputed
  `not has_raw_xenium_source(path)`, which is the inference the stamp exists to
  replace.

  It agreed by luck: an export writes none of the raw markers, so absence and
  declaration said the same thing. The moment an export also writes a raw-shaped
  file — the raw-format export this is a prerequisite for — the inference flips, the
  loader reopens rebuild paths on a dataset whose cache is the *only* copy of the
  data, and the recorded preamble silently reverts to `xenium(data_path)`, which
  reads the raw half and drops every derived layer the crop carried. That notebook
  runs and produces less, which is worse than the first-cell crash fixed in the
  previous release.

  `loader.is_cache_only()` is now the single question every reader asks, with
  `has_raw_xenium_source` left as the definition of the *inference* it falls back to.
  **Only a `True` stamp overrides**: a `false` one is read but not trusted to unset
  cache-only, because being wrong in that direction sends the loader down a rebuild
  path on data that cannot be rebuilt — the stamp may add certainty, never remove
  protection. A manifest that is missing, unparseable or silent makes no claim and
  never raises on the load path.

  `loader.cache_only_reason()` exists so the two *messages* stay true: Force
  Rebuild's refusal enumerates the absent raw files, which is a sentence contradicted
  by a declared export that has them.

- **The libGLX warning fired on an install that cannot have the problem, and told the
  user to run a command they did not have.** `gl_check` warned whenever the environment
  lacked the unversioned `libGLX.so` and the host had one. Its own docstring says "both
  halves of the collision must be present", but the abort needs *two glvnd builds in one
  process* — and the half it checked was the wrong one.

  Found by the fresh-machine install test, on the first `pip install palms[cnv]` ever
  run: that environment has **neither** `libGLX.so` nor `libGLX.so.0`, because PyQt6
  brings Qt inside its wheel and nothing pulls conda's `libglx`. PyOpenGL and Qt load the
  one host copy, nothing can disagree — and the viewer greeted a working app with a
  warning about an abort that could not happen. Measured beside the conda dev env, which
  has both names and was correctly silent.

  The check now also requires the environment's own `libGLX.so.0`. The message says which
  condition actually holds rather than naming a conda package unconditionally, and the
  remedy is derived from the running prefix: `mamba install -n <that env> libglx-devel`,
  with the `environment-linux.yml` route offered second and only where a conda env is
  running. It previously said `-n palms` whatever the environment was called — the one
  that reported this was named `test` — and pointed at `./scripts/install.sh`, which a
  user who installed from PyPI has no checkout to run.

## [1.0.0] — 2026-08-29

The first release. Work started on 2026-02-28 as a handful of scripts for looking at
Xenium output on Linux, where 10x's own Explorer has no build; `0.1.0` was a number in
`pyproject.toml` and was never published anywhere. Everything under **Development log**
below is part of this release — which is why, until now, no entry in this file carried a
version.

PALMS (Provenance-Aware Linking of Multimodal Spatial-omics) is a napari viewer for
Xenium 3.x output that brings transcripts, cells, histology and genomic overlays into one
coordinate space, and records every action the user takes as code that replays.

### Added

- **Reproducible-by-construction recording.** Every analysis is a node in a provenance
  DAG with a stable id, its code, and its dependencies; the notebook is *derived* from
  that graph by topological sort, so it respects dependencies regardless of the order
  actions were taken, across sessions. `ctx.run_step()` renders a template **once** and
  hands the same string both to `exec` and to the graph, so the code that ran and the
  code that is recorded are equal by construction rather than by discipline. Outputs are
  a flat `analysis.py` written live and an `analysis_notebook.ipynb` that replays from
  the raw Xenium output.
- **The claim is measured, not asserted.** `tests/test_notebook_replay.py` exports the
  graph as a real `.ipynb`, executes it in a clean kernel with `allow_errors=False`, and
  requires ARI exactly 1.0 on every clustering and identical top-N ranked genes;
  `scripts/verify_notebook.py` does the same against a real dataset and emits a JSON
  report with per-cell timing, per-clustering ARI and package versions.
- **Analysis tabs**, 26 of them in five groups: cell colouring and clustering (Leiden,
  Novae, imported labels), gene analysis and ranked genes, marker-gene and correlation
  plots, ROI differential expression, squidpy spatial statistics (neighbourhood
  enrichment, co-occurrence, ligand–receptor), CNV inference via **inferCNV** in-process
  and **CopyKAT** in a detached second conda env, UMAP in a linked napari window, and a
  Plots dock every figure passes through.
- **Multimodal overlays.** Landmark-based similarity-affine registration of H&E and ARMS
  images onto the transcriptomic frame, persisted with the session.
- **User-configurable analysis templates.** Every template ships as a `.tmpl` of plain
  scverse source with a declared contract, and is overridable per user **per block**, so
  untouched blocks keep tracking the shipped version and still receive upstream fixes.
  Tools → Templates shows the shipped text beside yours with a live preview of the exact
  string that would be executed; a notebook built from customised templates says so.
- **Data management that does not lose data.** A zarr cache with content-hash staleness
  (60–70% faster launches), crash-safe element writes that stage-journal-rename and keep
  the previous copy, a read-only cache verifier that works on a store too broken to open,
  a per-item dataset inventory with guarded deletion, `palms-rename-dataset` for moving a
  dataset without invalidating its provenance, and crop export.
- **Cross-platform install and docs.** Linux, macOS and WSL2; PyQt6 + napari ≥ 0.8;
  console scripts `palms`, `palms-preprocess`, `palms-build-cache`,
  `palms-rename-dataset`, `palms-fetch-references`, `palms-build-custom-segmentation`;
  documentation on Read the Docs; CI running the full suite plus a lint gate on every
  push.
- **1,392 tests** across 60 files, including source guards that fail if a fixed defect is
  reintroduced.

### Known limitations

- **Native Windows is not supported** (WSL2 is).
- The napari GUI proper and the spatial-analysis tabs have no automated coverage; that
  testing is manual.
- The annotation-neighbourhood and annotation-distance tabs are not yet recorded as
  steps — their synthetic virtual cells are sampled from a napari shapes layer the
  notebook cannot reach. Their figures are recorded as declared viewer-state notes.
- CopyKAT needs a second conda environment (`palms_copykat`, python 3.11 + R 4.3),
  because its R stack cannot coexist with the main environment's python 3.12.

---

## Development log

Every entry below predates the 1.0.0 tag and is part of it. They were written one per
working session, newest first, and are kept in full rather than summarised: several are
the only record of *why* a fix took the shape it did.

### 2026-08-29 (b) — dependencies declared and pinned

#### Changed
- **`insitucnv` is pinned to a tag instead of tracking `master`.** All three install
  sites — `environment.yml`, `environment-copykat.yml` and the `cnv` extra — carried a
  bare `git+https://…/insituCNV-copykat.git`, which resolves to whatever the fork's
  default branch happens to be that day. Two installs a week apart could run different
  CNV code while reporting the same everything else, which is not a dependency a
  reproducibility claim can rest on. The fork had no tags at all; `v0.2.0` now marks
  `172e753` and all three sites pin `@v0.2.0`. (The tag is ahead of the fork's own
  `version = "0.1.0"`, inherited from upstream and never bumped — the git ref is what
  pins the code, so `pip list` still says 0.1.0.)

  The three must move together: the viewer's env and `palms_copykat` share no
  site-packages, so a pin that lands in one and not the other has the CopyKAT worker
  running different code from the viewer that launched it.
  `tests/test_declared_dependencies.py` now fails if the three disagree, or if any of
  them loses its `@<ref>`.

#### Fixed
- **Seven modules the app imports were declared nowhere.** `novae` (Spatial → Domains)
  appeared in no dependency list at all; `platformdirs`, `pygments`, `geopandas`,
  `seaborn`, `tqdm` and `numcodecs` were reached only as transitive dependencies of
  napari, spatialdata, scanpy and zarr. Nothing was broken here, because every box that
  has ever run PALMS also had them — but the failure mode is bad: a `pip install palms`
  that resolves a napari or scanpy release which has dropped one of them fails at an
  import line naming neither PALMS nor the package that used to supply it.

  All six non-optional ones are now in `dependencies`, each with a comment saying which
  feature needs it. `novae` becomes its own extra (`pip install "palms[novae]"`, also
  folded into `full`) rather than a core dependency: it pulls torch, torch-geometric and
  lightning, and requires python ≥ 3.11. The Domains tab's "not installed" message now
  names the extra rather than `pip install novae`.

#### Added
- **`tests/test_declared_dependencies.py`** — a source guard in the idiom of
  `test_persistence_safety.py`. It walks the AST of every module under `src/palms`,
  maps each imported top-level name back to the *distribution* that provides it
  (`importlib.metadata.packages_distributions`, because `cv2` is `opencv-python` and
  `sklearn` is `scikit-learn`), and fails on anything missing from `pyproject.toml`'s
  dependencies or extras. Four of the seven above were found by writing it, not by
  reading the code. Exemptions are listed with reasons and a second test deletes one
  the moment nothing imports it any more.

### 2026-08-29 — reproducibility measurement

#### Fixed
- **inferCNV ran without the gene-mapping guard, and reproduced the exact defect that
  guard was written to prevent.** `418cf07` added `cnv_analysis.check_gene_mapping`,
  which refuses a panel where under 5% of genes map to genomic coordinates. Migrating
  inferCNV to `ctx.run_step` moved execution into the `genes.cnv_infercnv` template,
  and `run_cnv_pipeline` — the guard's only caller — was left serving the CopyKAT
  worker alone. So the inferCNV path called
  `prepare_cnv_input(add_gene_positions=True, drop_unmapped_genes=True)` and carried
  straight on; `n_genes_mapped` was computed *after* the run, only for the summary.

  The failure is silent and late. On the dataset that prompted this, a **mouse** panel
  against insitucnv's default **human** reference (infercnvpy Maynard 2020) mapped 8 of
  5,006 genes into 5 genomic windows. The run completed, clustered into 27 clusters and
  persisted a result; the first symptom came several steps later when the heatmap
  crashed. A notebook replay died at `cnv:infercnv` with `ArpackError -9`, which names
  nothing useful.

  The check now lives in the **template**, so it travels with the recorded code and a
  replayed notebook refuses for the same reason the GUI does — rather than dying
  opaquely somewhere downstream. The threshold is passed in from
  `MIN_MAPPED_GENE_FRACTION` rather than spelled again, so the template and
  `check_gene_mapping` cannot drift, and the CopyKAT worker still reaches the helper
  directly.

  **The decision is the mapped fraction, not the species.** A mouse panel given a mouse
  annotation table should run, and a human panel against a broken reference should not;
  species only sharpens the message. `tests/test_cnv_step.py` now runs the *rendered
  template through a real executor* and asserts it refuses and never reaches
  `run_infercnv` — the previous tests called `check_gene_mapping` directly and stayed
  green throughout the regression, which is how it went unnoticed.

#### Changed
- **`verify_notebook.py` treats 10x's own clusterings as inputs, not results.**
  `compare_clusterings` marked any persisted clustering the replay did not produce
  as `not_in_replay`, and `passed = not failures` counted those — so a session
  could reproduce every step it recorded and still report `passed: false`. On the
  `crop_6` demo dataset all four recorded clusterings replayed at ARI exactly
  1.000000 with identical labels, and the verdict was still `false`, purely
  because `graphclust` and nine `kmeans_*` columns came in with the dataset and
  nothing in the session produced them.

  Those are **inputs**: 10x computed them, the loader read them, no recorded step
  makes them, and a notebook replaying from raw output is not expected to. They
  now carry `status: "input"`, are excluded from `passed`, and are still listed
  in the report and named in the summary — "not compared" and "compared and
  agreed" are different statements, so neither is allowed to hide the other.

  **A viewer-derived clustering with no step behind it stays a failure**, which is
  the distinction that had to be made carefully: that is the defect
  `tests/test_clustering_recording.py` guards against, and the two are
  indistinguishable in `obs` — both are `clustering_<key>` columns. So the native
  ones are identified by the new `loader.is_native_clustering`, not by the absence
  of a node. Disk decides where it can: membership in `analysis/clustering/` is
  definitive. Only when a dataset has no `analysis/` folder at all does the name
  decide — a Crop Dataset export drops that folder while keeping the columns, and
  refusing to answer there would make every crop export look as though it had lost
  ten analyses.

### 2026-08-29 — stale cleanup

#### Fixed
- **A step whose recorded code does not vary with its settings could never be
  un-flagged as stale.** `ProvGraph.upsert` returns early when the code, deps and
  kind are unchanged, and that early return sat *above* `existing.stale = False` —
  so re-running the step, which is exactly what the `⚠ stale — input changed;
  re-run in the viewer` badge tells you to do, left the flag set. Permanently:
  nothing else clears it, and it round-trips through JSON, so it survived every
  restart.

  Found on the `crop_6` demo dataset, where `clustering:cnv_leiden_res0.2` showed
  as stale in the DAG immediately after the inferCNV run that produced it. That
  node publishes the CNV labels onto the table with a fixed three-line snippet
  carrying none of the run's parameters, while its parent `cnv:infercnv` carries
  all of them. Re-running inferCNV with a different reference therefore changed
  the parent (staling the child) and re-recorded the child byte-identically —
  which did nothing. Reproduced in four `upsert` calls.

  Only the node's own flag is cleared; descendants are deliberately untouched,
  since nothing about their input changed (`test_rerun_identical_is_noop` still
  holds). Safe because every caller that reaches the branch has just executed the
  step — `StepExecutor.run` upserts only after a successful `exec`, and the
  "ensure a node exists" paths never get there: `record_clustering` returns early
  when the id is present, and `ensure_normalized` / `ensure_spatial_neighbors`
  short-circuit before running their step.

  This mattered more once staleness became actionable: **"Select Stale Results"
  would have offered to delete a result that was in fact current.**

#### Added
- **Two ways to act on staleness, which until now was display-only.** A node is
  flagged `stale` when an ancestor was re-recorded with different code after that
  node last ran, so the artifact on disk was produced by an input that no longer
  exists. The Notebook tab drew a ⚠ badge and `dag_view` outlined the node orange;
  nothing else looked at the flag.

  - **Tools → Dataset → "Select Stale Results…"** ticks the rows that hold a stale
    step's results and hands them to the same confirm-and-delete path every other
    deletion goes through. It is a *selector*, not a second executor: `plan_deletion`
    still vets the batch, `assert_node_deletable` still has the last word, and
    `_remove_tree` remains the only function that touches the filesystem. The
    provenance graph is left alone, so the notebook still replays and recreates
    whatever was cleared — which is what `_plan_warnings` has always told a user who
    deleted a clustering column by hand.
  - **Tools → Notebook → "Drop Stale Nodes…"** prunes the steps themselves, in
    reverse topological order, after backing the graph up to
    `viewer_cache/prov_graph.backup_<time>.json` (a name `store_inventory` already
    protects from deletion). Before this there was **no way to remove a node through
    the GUI at all**: the per-cell "Delete" removed only the widget, and
    `_sync_from_graph` put the cell straight back.

  New pure module `utils/stale_results.py` is the **id bridge** between the two
  halves, and the only place the correspondence is written down — provenance ids are
  namespaced by the artifact they produce (`clustering:<key>`, `rank_genes:<key>`,
  `cnv:<backend>`) while inventory keys name a row on disk
  (`obs:table/clustering_<key>`, `sidecar:adata_cnv_cache_<backend>.h5ad`), and
  nothing on a `ProvNode` said which is which. Qt-free and filesystem-free like
  `store_inventory`, so it is tested on its own.

  Four things the implementation had to get right, each a real failure mode:

  - **Two monotonic guards assume nodes are never removed.**
    `app._load_prov_graph_items` rejects the sidecar whenever it holds fewer nodes
    than the session attr, and `session._build_session_attrs` refuses to write a
    graph smaller than the stored one — both stating "nothing in the GUI removes
    nodes". A prune that wrote only the sidecar was therefore undone at the next
    exit *and* the next launch, with a printed line about a "partial graph" as the
    only clue. Measured on a real store: pruning 2 of 4 nodes and writing only the
    sidecar loaded all 4 back. `tab_notebook._persist_pruned_graph` writes both
    copies in the same action, so the sizes agree and neither guard fires; neither
    guard needed weakening, and a *lone* shrink still means what it always meant.
  - **`uns['nhood_enrichment']`, `uns['co_occurrence']` and `uns['ligrec']` are
    single unkeyed slots** while the nodes that write them are keyed per clustering.
    A stale `nhood:A` and a fresh `nhood:B` name the same bytes, so clearing on the
    stale one alone would destroy a perfectly current result. Such a slot is only
    cleared when every member of its family is stale, and the report names the fresh
    step that spared it. `cnv:copykat_propagated` is excluded from that vote because
    it stores nothing of its own.
  - **The stale set is not dependency-closed.** `upsert` clears `stale` on the node
    it re-records while flagging that node's descendants, so a step recorded *after*
    a stale one is fresh and depends on it — and `ProvGraph.remove` refuses any node
    another still names. `plan_prune` excludes anything a fresh node can reach
    (transitively, not one hop) and emits the rest leaves-first, so every removal is
    legal when it happens.
  - **A moved dataset restales everything for nothing.** `app.py` re-emits the
    preamble for the current `data_path` on every launch, so the first launch after a
    manual `mv` flags every descendant stale — and a one-click clear would then
    delete every analysis in the dataset because a directory changed name. The
    confirm dialog says so and points at `palms-rename-dataset --repair`.

  Out of scope and reported as such rather than silently skipped: figures under
  `<data_path>/plots/`, which is outside every `deletable_roots()` directory by
  design, and CSV exports, which went to a path chosen in a save dialog and recorded
  nowhere. A stale step whose result is simply absent is distinguished from one this
  module has no rule for — reporting both as "nothing to remove" reads as a hole in
  the mapping table when it is not one.

  `ProvGraph.remove` also gains its first test coverage.

### 2026-08-29 — crop-export notebooks

#### Fixed
- **A Crop Dataset export's recorded notebook could never replay.** The preamble
  recorder always emitted `from spatialdata_io import xenium` /
  `sdata = xenium(data_path)`. A crop export has no raw 10x output — the zarr store
  the crop wrote *is* the data — so that call cannot work, and every notebook
  recorded on one failed on its **first** cell with
  `FileNotFoundError: …/cells.zarr.zip`.

  Nothing caught it because recording and persisting the graph both succeed; only
  *executing* the notebook fails, and the only thing that executes it is
  `scripts/verify_notebook.py`. Measured on the demo dataset shipped in the release
  bundle (`demo_data/crop_6`, an 18-node graph): the replay stopped at `preamble`
  after one cell.

  `_record_preamble` now branches on `loader.has_raw_xenium_source()` — the single
  existing definition of "can this be read from raw output" — and records
  `sd.read_zarr(data_path / "sdata_cached.zarr")` when there is none. The store is
  derived from `data_path` rather than recorded as its own absolute path, so
  `palms-rename-dataset`'s single `data_path = Path(r"…")` rewrite still moves the
  notebook with its dataset; a recorded absolute store path would have been rewritten
  in the first line and left stale in the second, which looks repaired and is not.

  This matters beyond one dataset: the crop export is how the project makes a small
  shareable dataset, so it is exactly the case where someone else replays the
  notebook. Note that re-launching the viewer on an existing crop export re-emits the
  corrected preamble, and the upsert flags every descendant stale — on `crop_6`, all
  16 of them. That is the graph telling the truth: those results were produced by a
  preamble that no longer exists.

### 2026-08-29 — pre-publication audit

#### Fixed
- **`scripts/capture_screenshots.py` hardcoded a collaborator's dataset.** The
  module constant was an absolute path into a working dataset directory, carrying a
  collaborator's name and a real slide ID in a tracked file in a repo about to go
  public — so it is described here rather than quoted. It is now `argv[1]` or
  `$PALMS_SCREENSHOT_DATASET`, resolved and validated *before* the napari import,
  so a missing argument fails in a second instead of after a ten-second Qt load.

  The docstring records why it must not become a constant again: two of the panels
  it captures — Tools → Dataset and Tools → Cache — print the dataset path into the
  widget, so whatever is passed ends up legible in `docs/` and on the wiki. That
  already happened, and `docs/screenshots/tab-dataset.png` and `tab-cache.png` show
  the path today. Recapturing them is tracked separately; an edit here cannot fix a
  PNG.

- **The control dock was still called "Xenium Controls".** The dock title, the View
  menu's "Show Xenium Controls" and the docs table row. Unlike `experiment.xenium`
  or "Xenium Explorer", that is the app naming *its own* panel — the exact use the
  rename existed to remove. It is now **"Controls"**: no product name at all, since
  the dock is unambiguous inside its own window. The five `setObjectName("xenium_*")`
  identifiers follow it to `palms_*`; nothing reads them back (there is no
  `objectName()` call in the tree), so they were pure identity strings.
  `Ctrl+Shift+X` is deliberately unchanged — a moved keyboard shortcut is a worse
  surprise than an arbitrary one.

- **`install_copykat` told users to run a command that has never existed.** Both its
  docstring ("Exposed as the `xenium-install-copykat` console script") and the error
  raised when the R install fails named an entry point absent from
  `[project.scripts]` under either name. Not wiring it up is correct, and now stated:
  the module needs rpy2 and R, which only the `palms_copykat` env has, so a console
  script installed by the main env would be a command that could never work. The
  failure message now points at re-running the analysis, which is what actually
  retries the install.

- **Two smaller leaks in `docs/`.** `pyqt6-migration.md` pointed at a directory
  outside the repo, which tells a reader nothing they can act on; it now says the
  upstream reports are drafted but not filed. `readthedocs-setup.md` pointed at a
  phase of a planning document that lives outside this repo — the only such reference
  in it — and reads the same without it.

#### Removed
- **`more_datasets.tsv`.** Ten candidate prostate/atlas datasets at the repo root,
  tracked since `64fc950` and read by no code, with ChatGPT citation artifacts
  (`:contentReference[oaicite:N]{index=N}`) inside the URLs. Moved to the private
  planning repo rather than deleted.

#### Verified
- The rest of the pre-publication audit came back clean: no secrets or credentials in
  tracked files; `.gitignore` does cover `manuscript/`, `data/` and `report.json`;
  `.claude/` is untracked; `LICENSE` is the MIT one `pyproject.toml` declares; both
  open GitHub issues read acceptably in public; and the only absolute-path hits in the
  whole tree were the two fixed above.

### 2026-08-27 — rename to PALMS

#### Changed
- **The project is now PALMS — Provenance-Aware Linking of Multimodal
  Spatial-omics.** The old name borrowed 10x Genomics' instrument trademark for
  the identity of an independent tool, and undersold the scope: this is not a
  transcriptomics viewer, it is where transcriptomics, histology and genomic
  overlays are brought into one coordinate space with every step recorded as
  replayable code. It sits beside ARMS (Adaptive Resolution Multiscale Spatial
  DNAseq) in the same ecosystem. Renamed now because the citable release would
  otherwise have minted the old name into a DOI permanently.

  Import package `xenium_viewer` → `palms`; distribution `xenium-viewer` →
  `palms`; the six console scripts to `palms`, `palms-preprocess`,
  `palms-build-cache`, `palms-rename-dataset`, `palms-fetch-references`,
  `palms-build-custom-segmentation`; conda envs to `palms` and `palms_copykat`;
  the six `XENIUM_*` environment variables to `PALMS_*`.

  **No shim.** Nothing was on PyPI and the version is 0.1.0, so no external
  consumer can break — and a compatibility package would have kept the
  proprietary name in the tree, which is the thing being removed.

- **References to the Xenium *data format* are untouched, deliberately.**
  `experiment.xenium`, `spatialdata_io.xenium()`, "Xenium 3.x output" and the
  comparison to Xenium Explorer all describe what the tool reads, not what it
  is called. A blind substitution would have broken the loader; every
  replacement was anchored on an identity token instead, and the format
  references were diffed before and after to prove they did not move.

#### Fixed
- **Three names that reach onto users' disks keep working.** A rename that
  silently orphans existing work is a data-loss bug, not a cosmetic change:
  - Template overrides written to `~/.config/xenium-viewer/templates/` still
    resolve. `search_path()` reads the pre-rename directory below the current
    one, and only while it exists, so it retires itself.
  - `.tmpl` files whose banner reads `# xenium-viewer template` still parse —
    the banner is prose, not a header field. Now pinned by a test, because a
    stricter parser would deactivate every old override at once.
  - `xenium_viewer.log` in an existing dataset directory is still listed and
    deletable in Tools → Dataset. Dropping the key would not have hidden the
    file, it would have relabelled it "original 10x output, never modified by
    the viewer" — false, and un-deletable.

### 2026-08-26 — doc counts

#### Fixed
- **Four counts in the docs were wrong, three of them badly.** The test suite is
  **1318 tests across 56 files**, not the "~320" `CLAUDE.md` had carried for
  months — off by a factor of four. There are **26 tabs** in 5 groups, not the
  11 in `CLAUDE.md` or the 21 in `README.md`; `app.py`'s `addTab` calls are the
  authoritative count and `src/xenium_viewer/tabs/*.py` agrees with them.
  `README.md`'s per-group listing was also short: the Tools group showed 3 of its
  7 tabs, omitting Crop Dataset, Dataset, Cache and Templates. `app.py` is ~1800
  lines, not ~1300.

  Counted rather than recalled: `pytest --collect-only -q` (1314 passed + 4
  skipped confirms it), `grep '\.addTab(' src/xenium_viewer/app.py`, `wc -l`.
  `CLAUDE.md` now says how to re-measure the test count instead of stating a
  number that will drift again.

#### Removed
- **`tabs/__init__.py` no longer re-exports `build_*_tab`.** It aliased
  `build_tab` from 11 of the 26 tab modules — the 11 that existed when it was
  written — and nothing has ever imported them: `app.py` imports each module
  directly, *inside* the function that builds the panel, precisely so the
  napari-heavy tab modules load only when a viewer is being constructed. The
  eager re-exports quietly undid that for anything touching the package, which
  includes three test modules doing `from xenium_viewer.tabs import tab_cache`.
  Completing the list to 26 would have made that worse, so the dead surface is
  gone and the docstring says what to import instead.

### 2026-08-26 — plots

#### Fixed
- **Every control label in the Xenium Controls panel was invisible.** `make_tab`
  added `w.native` — the *bare* Qt control. magicgui keeps a widget's caption in a
  wrapper only a `Container` builds, so all 72 captions were discarded at layout
  time. The Clustering tab rendered as six anonymous controls reading `15`, `40`,
  `1.00`, `igraph`, `2`, `2000`, with nothing to say which was which; `CheckBox`
  and `PushButton` looked fine only because Qt paints their text on the control.

  `_helpers.labelled()` restores the caption with a one-widget `Container` —
  public API, and measured identical to the private `_LabeledWidget` (same 142px
  minimum, same 22px height). It returns the bare control for a `ButtonWidget`
  (mirroring what magicgui's own `Container` does, so a checkbox does not show its
  text twice) and for a widget with no label, which magicgui reports as `''` — it
  never derives one, so nothing sprouts a caption nobody wrote.

  `tab_annot_nhood` and `tab_annot_distance` hand-rolled the `QVBoxLayout` +
  `QScrollArea` that `make_tab` already is, which is exactly why their 10 labels
  would have stayed invisible after a `make_tab`-only fix; both now call the
  helper. No other call site changed.

  Found by screenshotting the running viewer over the new `--mcp` bridge. The
  1315-test suite passed throughout, because none of it renders a tab — so
  `tests/test_control_panel_size.py` gained three guards: the label survives
  `make_tab`, a `CheckBox` is not captioned twice, and an unlabelled widget stays
  unlabelled.

#### Changed
- **Control captions read as English rather than scanpy parameter names**, now
  that they render at all: `n_neighbors` → **Neighbours**, `n_pcs` → **Principal
  components**, `flavor` → **Clustering backend**, `n_top_genes` → **Highly
  variable genes**. Each carries a tooltip naming the template parameter its value
  lands in, so the correspondence with the exported notebook — which the captions
  used to *be* — is kept rather than lost.

  Captions that repeated a section heading or the tab title are trimmed
  (`CellTypist Model` → **Model** under the "CellTypist Annotation" heading;
  `UMAP pt size` → **Point size** on the UMAP tab), the three CNV controls whose
  captions all contained "resolution"/"backend" are now distinct (**Cluster
  resolution**, **Inference backend**, **Heatmap backend**), and `N neighbors` /
  `N neighbours` are spelled the one way. Docs updated to match; the two generated
  reference pages still name the parameters, since they describe code.

#### Added
- **`--mcp`: a dev-only MCP bridge onto the running viewer** (`src/xenium_viewer/dev_mcp.py`,
  `mcp` extra). Starts `napari-mcp`'s `NapariBridgeServer` over the viewer that already has
  the dataset loaded, so an assistant can take **full-window** screenshots — docks included,
  which is where this app's UI defects actually live — and run `execute_code` on the Qt main
  thread. `_push_to_console` now also publishes `ctx` as `viewer._xv_ctx`; that is the one
  place a *fresh* ctx reaches an interactive surface, so it survives a dataset switch, where
  hanging it anywhere else would leave a stale object behind.

  Off by default, out of `environment.yml` and out of `full`: the bridge is an
  unauthenticated localhost server that executes arbitrary Python in the viewer process.
  `napari_mcp`'s own `init_viewer` is the wrong entry point here — it would create a second,
  empty viewer.

  `fastmcp` is pinned `<3` in the extra. napari-mcp declares only `>=2.10.3` but reaches
  into `FastMCP._tool_manager._tools`, which fastmcp 3 renamed, so a fresh install resolves
  3.x and the bridge dies at start with `'FastMCP' object has no attribute '_tool_manager'`.
  Measured, not guessed.

  Known limits, both inherent: every bridge call marshals onto the Qt main thread, so one of
  this app's modal dialogs blocks the bridge until a human clicks it, and long main-thread
  work (cache build, pyramid load) hits `NapariBridgeServer`'s 300 s timeout.

- **A Plots dock** (issue #35). Every figure the viewer produces — dotplots, UMAPs,
  neighbourhood-enrichment heatmaps, co-occurrence curves, L-R dotplots, CNV heatmaps,
  the provenance DAG — now appears in one gallery at the bottom of the window, newest
  first, each card carrying a thumbnail, the paths it was written to, and Open / Save
  as… / Remove. Hidden until the first plot; **View → Show Plots** (`Ctrl+Shift+P`)
  toggles it. Capped at 20 figures, because the viewer produces one per action and a
  rank-genes heatmap is not small.

  **Open** draws into the dock's own `FigureCanvasQTAgg` rather than calling
  `plt.show`, which is what makes display independent of the process-wide matplotlib
  backend — see the `matplotlib.use` fix below.

- **UMAP coloured by gene expression** (issue #34). The UMAP tab gained a multi-gene
  picker (up to 15) with a colormap and a column count: one gene renders as a single
  panel with its colour bar, several as a grid. "Save UMAP Plot…" became **Plot UMAP by
  cluster** — same template, no file dialog, and it no longer builds a throwaway AnnData
  by hand. The gene list is session-persisted.

  New `umap.plot` template. Its `embed.xenium` block reads
  `analysis/umap/gene_expression_2_components/projection.csv` — **the same file the
  viewer's UMAP window reads** — so a replayed notebook reproduces the figure that was on
  screen. The old `plot:umap` node recorded `sc.pp.neighbors` + `sc.tl.umap` instead,
  which produces an equally valid but *different* layout; `docs/Tab-UMAP.md` had to warn
  about the discrepancy. `embed.recompute` remains for a Crop Dataset export, which has
  no `analysis/` folder, and says so in a comment. Verified against a real dataset:
  `reindex(obs_names)` matches 318,617 of 318,708 cells — the 91 Xenium's UMAP omits,
  which become NaN rows that `sc.pl.umap` skips.

#### Changed
- **One save policy: `<dataset>/plots/`, PNG and PDF.** Preferences → Plot format now
  offers **PNG + PDF** (default), PNG, PDF and SVG, and every plot is written in each
  chosen format. Previously six sites used `ctx.auto_save_plot` (one format), CNV
  hard-coded PNG+PDF, UMAP and marker genes asked via `QFileDialog`, the three volcano
  batches hard-coded PNG into a directory of the user's choosing, and the rank-genes
  panel plot was **never saved at all** — the one figure in the app a user could not keep.

  Names are now keyed by what the figure is about (`dotplot_leiden_r1.0`,
  `nhood_enrichment_graphclust`, `umap_EPCAM_KRT5`). They were not, so a second
  neighbourhood plot silently overwrote the first.

  Volcano batches still ask for a directory — an N×N run is the one output people
  routinely want elsewhere — but default to `<dataset>/plots/volcano_<key>/` and honour
  the format setting. They stay out of the gallery on purpose: fifty figures would bury
  everything else.

- **Recorded `savefig` paths name the file that was actually written.** Every
  hand-written `plot:*` node recorded a bare relative guess (`"dotplot.svg"`,
  `"umap.png"`, `"nhood_enrichment.svg"`) that matched nothing on disk; the volcano nodes
  recorded `Path("<basename of the chosen directory>")`. They now emit the real path
  relative to the dataset, preceded by the `mkdir` that a replayed notebook needs —
  `savefig` does not create its parent, so those cells would have failed on a clean
  checkout anyway.

- `genes.marker_plot` and `genes.correlation` take `paths: list` in place of
  `path: str`, and marker_plot's `save.plain`/`save.dpi` pair collapsed into one `save`
  block (20 assemblies → 10). A user override of those blocks is flagged **stale** by the
  existing upgrade machinery rather than silently doing the wrong thing.

- The two annotation plots (nhood, distance) recorded nothing at all and now record a
  `NOTE`. Not a `TERMINAL`: both are computed over virtual cells sampled from shapes
  drawn in the viewer, which a notebook cannot reach, so a cell calling
  `plt.gcf().savefig(...)` with no preceding plot would replay as a silent no-op writing
  an empty figure. `NOTE` says "not code, and not meant to be" and is counted apart from
  the comment-only punch list.

- `render_dag` no longer takes a `path=` — it returns the Figure and lets the caller
  write it, the convention `make_cnv_heatmap` already documented.

#### Fixed
- **`matplotlib.use('Agg')` leaked out of two workers and disabled plotting for the rest
  of the session.** `tab_umap` and `tab_marker_genes` set it globally, from a background
  thread, and never restored it — so once a user saved a UMAP or ran any marker plot,
  every later `plt.show(block=False)` in the process became a silent no-op and figures
  from *unrelated tabs* simply stopped appearing, with no error anywhere. This is a large
  part of what issue #35 described as plots "showing up in completely different
  interfaces". Both calls are gone, and nothing displays through pyplot any more.

- **Two plots blocked the viewer.** `tab_annot_nhood` and `tab_annot_distance` called
  bare `plt.show()`, which spins its own event loop: the window had to be closed before
  the viewer responded again.

- The CNV tab's "Show Heatmap" button now shows the heatmap. It built the figure, wrote
  two files and closed it without ever putting it on screen.

- **The rank-genes dotplot crashed the moment it reached the Plots dock.**
  `sc.pl.rank_genes_groups_dotplot(return_fig=True)` returns a scanpy `DotPlot`, not a
  `Figure` — and `make_rank_genes_dotplot` was *annotated* as returning a `Figure`, which
  is how it came to be handed straight to Qt. `DotPlot` has `savefig`, so the old
  save-only path never noticed; drawing one raises `AttributeError: 'DotPlot' object has
  no attribute 'set_canvas'`. New `fig_render.to_figure()` resolves it, once, in
  `ctx.show_plot`. Two wrinkles it has to handle: a `BasePlot`'s `.fig` is `None` until
  the plot is built (so `get_axes()`, which builds it and is idempotent), and
  `BasePlot.savefig` writes `plt.gcf()` rather than `self.fig` — resolving to the concrete
  Figure removes the chance of saving whichever figure happened to be current. Anything
  else raises a named `TypeError` instead of failing inside Qt. The misleading annotation
  is gone.

- **Closing the Plots dock destroyed it, and the View-menu toggle then did nothing.**
  napari's dock title bar has an "×", and it does not hide the dock: it calls
  `destroyOnClose` → `Window.remove_dock_widget`, which reparents the gallery to `None`
  and `deleteLater()`s the dock. So one click left a dangling C++ pointer, and
  **View → Show Plots** called `setVisible` on it, raised `RuntimeError` into a bare
  `except: pass`, and looked broken. That is both halves of the report — the dock
  "disappearing" and the toggle not working.

  The dock is now disposable and the gallery is what persists: `ensure_plots_dock`
  re-creates a dock around the surviving panel whenever the old one is gone, so the
  figures are never lost. The "×" is shimmed to *hide* rather than destroy (closing a
  gallery should not throw away what is in it). And `reveal_plots_dock` raises the window
  and, if a floating dock has been dragged somewhere no connected screen reaches, docks it
  back into the main window rather than re-showing it out of sight.

- **Two clusters could not share a display name.** `.cat.rename_categories()` raises
  *ValueError: Categorical categories must be unique*, so naming clusters 0 and 2 both
  "Tumour" — an ordinary thing to want, meaning "these are the same cell type" — failed
  the whole step. Both templates that relabel (`umap.plot`, `genes.marker_plot`) now use
  `.map()` with a mapping and `dict.fromkeys` to dedupe, which **merges** the clusters,
  which is what was meant; cluster order is preserved in the legend. `categories` is a
  `dict` (original id → display name) rather than a positional list, so it also cannot
  silently mis-align if the category order shifts. `genes.marker_plot` carried the same
  defect before this branch and is fixed with it.

`tests/test_plot_consistency.py` is the source guard that keeps this from unravelling one
tab at a time: it parses every `tabs/tab_*.py` and fails on any `plt.show`, any
`matplotlib.use`, and any `savefig` outside `plot_output.save_figure` — plus a check that
`<data_path>/plots` has exactly one definition, where five modules used to build it by
hand. New `tests/test_plot_output.py` (pure logic) and `tests/test_plots_panel.py`
(offscreen Qt, including the dock-minimum bound from the August resize fix, which a
second dock inherits).

### 2026-08-26

#### Fixed
- **The Xenium Controls dock could not be resized.** Dragging the separator between the
  panel and the canvas did nothing: the control panel reported a minimum size of
  **536 × 534 px** (up to ~617 wide once the Templates tab's "Take new default for
  changed blocks" button appeared), so there was nothing left to give. Now **131 × 175**,
  and the dock itself reports 127 × 190 instead of ~540 × 534.

  This is the March fix at `CHANGELOG.md:2889` (*Console not resizable when opened*)
  recurring on the width axis, and for the same reason. The panel is a `QTabWidget` of
  `QTabWidget`s, and **a stacked widget's minimum is the maximum over all its pages,
  hidden ones included** — so one unwrapped page sets the floor for the whole dock, even
  when the user never opens it. `make_tab()` wraps page content in a `QScrollArea`
  precisely to stop that (a scroll area reports a fixed ~68 px minimum whatever it
  holds); two pages added after March bypassed it:

  - **Notebook** (528 × 107) — its six-button toolbar sits *outside* the tab's own scroll
    area, and a row of buttons cannot shrink below its labels. The toolbar is now a
    `toolbar_row()`: still pinned above the cells, but scrolling sideways when narrow, so
    it costs 68 px instead of 528. Wrapping the whole tab instead would have worked, but
    would have scrolled the toolbar out of view.
  - **Templates** (389 × 468) — `build_tab` returned a bare `QWidget` of nested
    `QSplitter`s with no scroll area anywhere. Now returned through `scrollable()`, with
    its four action buttons in a `toolbar_row()` so a narrow dock does not hide
    "Save && Activate" behind a scrollbar.

  New `_helpers.scrollable()` (which `make_tab` now uses, so there is one implementation)
  and `_helpers.toolbar_row()`. `tests/test_control_panel_size.py` measures both pages and
  the assembled panel, including the Templates tab with every button shown — the worst
  case is the state a user is in when they visit the tab to resolve an upgrade.

  Ruled out by measurement, and recorded so nobody re-derives it: napari 0.8's dock widget
  (`setMinimumWidth(50)`, an 8 px separator, `Movable` set — an empty panel gives a 121 px
  dock), the central canvas (`minimumSizeHint` 8 × 4), the PyQt6 enum rewrite, and every
  `make_tab`-wrapped page (a 580 px pixmap inside one still reports 68 px).

### 2026-08-25 — installable on macOS

#### Fixed
- **`conda env create -f environment.yml` failed to solve on macOS.** The PyQt6
  migration added `libglx-devel` to fix a Qt6 startup abort on remote X displays, but
  that package is part of conda-forge's linux-only glvnd stack — there is no
  `osx-arm64` or `osx-64` build — so the solve failed on a Mac before it reached the
  Qt stack at all. **PyQt6 itself was never the problem**: `pyqt6` publishes 6.8.1 and
  6.11.0 builds for `osx-arm64` on py310–py314, and every other dependency
  (napari 0.8.0, spatialdata 0.8.0, squidpy 1.8.2, scanpy 1.12.3, zarr 3.1.5 within
  the `<3.2` pin, `py-opencv` on a `qt6_` build) resolves there too. So this was
  packaging, not a reason to walk the Qt backend back.

  conda env files have **no platform selectors** — `# [linux]` is a conda-*build*
  `meta.yaml` feature that `conda env create` ignores — so the fix is a split:
  `environment.yml` is now cross-platform, and `environment-linux.yml` is a small
  *overlay* carrying `libglx-devel` alone. An overlay rather than a second full copy,
  so the shared dependency list stays in exactly one file and the two cannot drift.

  `scripts/install.sh` does both steps, branching on `uname -s`, so the install is one
  command on every platform. **WSL2 is covered for free** — `uname -s` is `Linux`
  there and conda installs `linux-64` packages, so it takes the Linux branch.

  Skipping the overlay on Linux is no longer silent. New `utils/gl_check.py`, called
  from `app.py` *before* `import napari` (the import that aborts), checks whether the
  env is missing the unversioned `libGLX.so` while a host copy exists, and prints the
  `env update` command. It deliberately only reports: preloading conda's `libGLX.so.0`
  with `RTLD_GLOBAL` does not help, because PyOpenGL still `dlopen`s the host's
  unversioned copy as a separate mapping. Both halves of the collision must be present
  before it says anything — a warning that fires on every working machine is a warning
  people learn to scroll past. Note the check is for the *unversioned* name on both
  sides: `ctypes.util.find_library('GLX')` returns `libGLX.so.0`, which every box with
  working graphics has, so it cannot answer this question.

- **A live CopyKAT worker was reported as `interrupted` on macOS.** `tab_cnv`'s
  liveness probe reads `/proc/<pid>/cmdline` with a documented `os.kill(pid, 0)`
  fallback for non-Linux — but on macOS `open("/proc/...")` raises `FileNotFoundError`,
  an `OSError` subclass caught *first* by the "no such process" branch, so the fallback
  was unreachable. The platform check now comes before `/proc` is touched.

- **The whole test suite ran offscreen on macOS.** `tests/conftest.py` forced
  `QT_QPA_PLATFORM=offscreen` when `DISPLAY` was unset, and `DISPLAY` is an X11
  variable that is absent on a Mac with a perfectly good screen. Same flaw in
  `test_units.py::_has_real_display()`, where it silently skipped every
  `requires_display` test. Both now special-case `darwin`.

#### Changed
- The ICE/X11 startup patch in `app.py` is gated on Linux. It was setting
  `SESSION_MANAGER=''` unconditionally, including on macOS, which has no session
  manager.
- CI gains a **macOS (arm64) leg** that builds `environment.yml` *without* the overlay
  and asserts the resolved Qt backend, then runs the suite under `cocoa`. A
  Linux-only matrix cannot see this class of breakage — the first report of it was a
  user whose MacBook could not create the env. The Linux leg now applies
  `environment-linux.yml` after creating the env, so it keeps testing the stack a real
  Linux user actually gets.
- README, `docs/Installation.md`, `docs/Home.md` and the `pyproject.toml` classifiers
  now say Linux/macOS/WSL2 rather than Linux. **The CopyKAT second environment stays
  Linux-only** and is documented as such: `r-dlm` is published only on the Anaconda `r`
  channel, which has no `osx-arm64` builds. inferCNV runs in the main env and is
  unaffected everywhere.

#### Known limitations on macOS
- Per-element memory reporting degrades to "unknown": `utils/mem_probe.py` reads
  `/proc/self/statm` and `/proc/self/status`, and resolves `malloc_trim` from libc.
  All three already fail soft, so this is a missing number rather than an error.

### 2026-08-25 — quiet the units warning

#### Fixed
- **Loading a dataset no longer prints `Inconsistent units across layers` once per
  layer.** Roughly twenty lines on a real load, and spurious every time: a new layer
  arrives carrying napari's default pixel units while every existing layer is already
  in micrometres, and napari's own canvas handler triggers a draw inside that window —
  `on_draw -> _update_scenegraph -> add_layer_visual_mapping -> _update_world_units` —
  which reads `viewer.layers.extent.units`, finds the disagreement, and warns. By the
  next draw the layer is stamped and everything agrees; measured, `extent.units`
  settles to micrometres and the scale bar renders correctly. The message describes a
  state that has already stopped being true.

  Connection order does not fix it (`position="first"` on `inserted` changes nothing —
  the draw is not ordered by it), and setting units at the ~20 `add_*` call sites is
  the design `utils/units.py` exists to avoid, because one missed site leaves that
  layer *misplaced* rather than merely mislabelled.

  So the message is dropped for the span of one insertion — `inserting` to the end of
  the `inserted` handler — and only that message. A genuine mismatch outside that
  window still reaches the user; there is a test for each direction.

  **The one real failure this could have masked is now reported better than the
  message it replaces.** A layer that refuses the µm scale really does leave the world
  inconsistent, so `apply_to_layer` returns whether it took, and
  `reporting.report_layer_scaling_failure` names the layer and says what it means —
  which napari's message never did. Once per layer per session.

  The patch is lifted when the window closes rather than left installed behind a depth
  counter. The first version kept a permanent wrapper around whatever `show_warning`
  was bound at first use, which made it order-dependent with anything else that
  rebinds the name — it surfaced as a test that passed alone and failed in a suite.

  Note for anyone adding tests here: a napari `Viewer` needs a GL context and CI's
  runner has none. It is worse than a clean failure — the `Viewer` *constructs*, then
  vispy returns `None` from `glGetParameter` on the first draw and the run ends in a
  segfault (exit 139), so there is nothing to catch. The two canvas-level tests are
  therefore gated on a real `DISPLAY` and never build a Viewer in CI at all. Because a
  skipped test guards nothing, the same property is covered canvas-free by driving the
  wiring directly — and that one was verified to fail with the fix neutered, under
  `QT_QPA_PLATFORM=offscreen`, which is exactly how CI runs.

### 2026-08-25 — a cropped dataset opens like a dataset

#### Fixed
- **A store with no `viewer_session` restored nothing at all.** `restore_fn` — the call
  that fans out to every tab's `restore_session`, and the only thing that adds the H&E,
  ARMS, external-image and patch layers at startup — sat inside `if session is not None`,
  under `if not no_cache and zarr_path.exists()`. The `elif` meant to catch a store with
  elements but no session was attached to the *outer* `if`, so it could never run while
  the zarr existed. A freshly exported crop therefore opened with no H&E, no ARMS, no
  ROIs and no cluster names, while every `load_*_from_sdata` result computed immediately
  above it was discarded. Restore is now unconditional and the session is an overlay on
  top of what the elements already say. Not crop-specific: it hit any cache built by
  `xenium-build-cache`, any store whose session node was deleted in Tools → Dataset, and
  any recovered cache.
- **An overlay can now place itself** (`utils/registration_seed.py`). Placement used to
  come only from session state, so a store without one added the image at the origin.
  The element's own transform is read back when the session offers no affine — and
  strictly all-or-nothing: the element stores `fine @ flip` while the session stores
  `fine` and re-derives the flip, so taking the affine from one and the flip from the
  other applies the flip twice. When the element is the source the flips are forced
  `False` and the shape is taken from the element.
- **Restoring could overwrite a registration on disk.** `tab_arms._on_arms_restored`
  calls `_save_arms_affine_to_sdata()` on every restore; with no session the affine is
  `None`, the composed transform is the identity, and that identity was written over the
  real one — losing a registration during a launch that only meant to read. Both overlay
  tabs now refuse to write an identity over a stored transform. This had to land *before*
  the restore fix above, which would otherwise have destroyed the one overlay a real
  export got right.
- **The crop sliced overlays out of the wrong region.** `crop_overlays` decided an
  overlay's frame from its stored transformation, but a registration is only written
  there when you re-register or flip — a real dataset's `he_image` declared `identity`
  while its affine lived in the session. The Xenium crop box was therefore mapped
  straight into H&E pixels: measured on an export, `he_image` came out as
  `(3, 2254, 16371)` — rows 11436–13690 × cols 6721–23092 of a `13690 × 23092` slide, a
  strip off the bottom-right corner unrelated to the cropped tissue. The external image
  went the same way. Frames are now resolved from the live viewer
  (`utils/crop_state.py`), falling back to the element, then to identity.
- **`arms_tiles` and `patch_*` overlays were read as Xenium data.** Both hold an
  overlay's pixels while declaring identity on disk, because they are drawn with a
  *linked* layer's affine that is never written to the element. All three `patch_*`
  elements in a real export were dropped as "nothing inside the crop", and 37 of 288
  `arms_tiles` survived by numeric coincidence at meaningless coordinates. They are now
  clipped in their companion image's frame.
- **An overlay whose frame cannot be established is no longer sliced.** An identity
  frame is credible only when the raster is on the morphology grid; otherwise the export
  skips it and names both extents in the summary and the provenance note. Slicing on a
  guess is the irreversible step — a transform can be corrected later, absent pixels
  cannot.
- **Image-space landmarks were offset by their companion's slice origin.** The exported
  raster is re-origined so its pixel `(0, 0)` is the source's `(r0, c0)`, but landmarks
  were passed through in *unsliced* pixels and are drawn with the image's affine. The two
  agreed in the source store only because nothing was sliced; in an export they differed
  by exactly `(r0, c0)` — 148 px on the regression fixture.
- **A cropped export now carries a `viewer_session` and its own provenance graph**
  (`utils/crop_session.py`), with paths rewritten to the export and a `crop_export:<name>`
  node built by the same function the source dataset's graph uses. Registration affines
  are deliberately *not* copied into it — the element owns placement, and a second copy
  drifts as soon as anyone re-registers inside the export. `arms_geojson_path` /
  `arms_csv_path` are nulled, because the ARMS restore falls back to those absolute paths
  and would pull the whole slide's tiles into a crop.
- `_carry_over_clusterings` no longer writes a `cluster_labels_*` column consisting
  entirely of empty strings when a clustering's names have all been cleared.

### 2026-08-24 — crop carries the overlays

#### Added
- **Crop & Export now carries every registered overlay and drawn region** (issue #27).
  A crop used to keep five core elements and silently drop the rest, so exporting a
  core threw away the registered H&E, the ARMS tiles and the ROIs — the work that
  made the dataset worth cropping. `utils/crop_overlays.py` brings them along, each
  cropped to the same region, and a checkbox on the tab turns it off for a bare
  core-only export.

  **What decides how an element is carried is the element itself.** An overlay with
  its own affine holds coordinates in its own pixel space, so its geometry is left
  alone and the crop is composed into its transform; an element with no
  transformation is already in Xenium pixels, so its geometry moves into the crop
  frame. That is derived rather than listed, because the loader's own comments
  record that a fixed name list has already gone stale once, and `ext_*` / `patch_*`
  elements are named per file precisely so they cannot be enumerated up front.

  > **Superseded on 2026-08-25** — see the entry above. The element turned out not to
  > be a trustworthy authority: a registration reaches the element only when you
  > re-register or flip, and `arms_tiles` / `patch_*` never get one. Frames are now
  > resolved from the live viewer, with the element as the fallback. The paragraphs
  > below still describe the geometry correctly; only the *source* of the frame changed.

  Three things that are wrong in ways that raise no error, each pinned by a test:

  - **H&E-space landmarks carry no transformation on disk.**
    `save_overlay_affine_to_sdata` is only ever called with an image or patch
    element's name, and `_save_he_affine_to_sdata` writes only to
    `images["he_image"]` — so `he_he_landmarks`, `arms_he_landmarks` and `*_image_lm`
    arrive looking exactly like Xenium-space geometry. Read by the rule above they
    would be translated by the crop origin, corrupting the registration they exist
    to reconstruct. They are recognised by name and passed through verbatim.
    (Since 2026-08-25 they also move by their companion raster's slice origin — the
    exported image is re-origined, and they are drawn with *its* affine.)
  - **Landmarks are never filtered.** Dropping a landmark that falls outside the crop
    re-fits the registration: the affine is a least-squares fit over the whole set,
    so a subset yields a different transform than the one the export ships with.
  - **`cell_circles` is stored in microns**, with a 1/pixel_size scale, and patch
    overlays in their source image's pixels. Carrying either whole would ship circles
    and patches covering cells the exported table no longer contains, so the crop
    region is mapped down into each element's own frame and the clip happens there.

  Rasters stay lazy end to end — the level-0 dask array is sliced, never computed —
  for the same reason the core crop does; an overlay is no smaller than the
  morphology image it was registered against.

  `crop_and_export` now returns `(path, overlay_notes)`. An overlay that cannot be
  carried is skipped and named in the summary dialog and in the recorded provenance
  note, rather than failing an export the user did get most of.

#### Notes
Tested against synthetic elements (`tests/test_crop_overlays.py`, 21 tests) rather
than a dataset: a rotated, scaled, registered overlay is exactly what no local
dataset has, and the failure mode is geometric — a misplaced overlay still exports
cleanly, still opens, and still looks like a registered H&E. The tests assert where
a known pixel lands, and the two subtlest guards (the image-space landmark rule and
the (y, x) → (x, y) transposition when mapping the region into a rotated frame) were
each confirmed to fail when mutated.

### 2026-08-24 — PyQt6

#### Fixed (post-review)
- **napari aborted at startup under Qt6 on a remote X display** — `Could not initialize
  GLX`, `SIGABRT`, no Python traceback. This nearly got the migration shelved as a
  Qt6-versus-remote-X incompatibility. It is neither: it is a library-loading conflict, and
  the display is incidental.

  conda's `libglx` provides only `libGLX.so.0`; the unversioned `libGLX.so` lives in
  `libglx-devel`. PyOpenGL's loader tries the unversioned name **first**, with `RTLD_GLOBAL`,
  so without that package it misses the env and loads the *host's*
  `/usr/lib/.../libGLX.so`. Two glvnd builds then share one process, Qt's GLX plugin
  resolves `glX*` across both, `glXGetVisualFromFBConfig()` returns NULL for a config the
  other library allocated, and Qt calls `qFatal`.

  Reduced to a reproduction where only the import order differs — `import PyQt6.QtGui` then
  `import OpenGL.GL` aborts; the reverse works — with no napari and no vispy involved.
  napari triggers it because it imports Qt and then PyOpenGL; bare vispy never imports
  PyOpenGL at all, which is why vispy rendered fine on the same display.

  `environment.yml` now depends on **`libglx-devel`**. `/proc/self/maps` goes from two
  `libGLX.so.0.0.0` mappings to one, and the viewer launches, restores its session and
  renders on the display that previously dumped core.

  **PyQt5 maps the same two copies and merely tolerates them**, so this was latent long
  before the migration and will not appear as a regression in a PyQt5 run. No environment
  variable fixes it — roughly fifteen were tried — because none creates the missing name
  inside the env. `tests/test_qt_backend.py` guards the dependency.

- **CI came up on PyQt5 despite `pyqt6` in `environment.yml`** — `napari 0.8.0 |
  qtpy PyQt5 5.15.15`, with PyQt6 6.8.1 also installed. Two things combined:

  **qtpy resolves in the order PyQt5, PySide2, PyQt6, PySide6**, so any environment
  that merely *contains* PyQt5 runs on it, whatever the environment file asked for.

  **conda-forge's `matplotlib` metapackage bundles a Qt binding, and which one
  depends on the version resolved** — 3.9.1 depends on `pyqt >=5.10` (Qt5), while
  3.9.3+ depend on `pyside6`. CI landed on 3.9.1 and got Qt5; a local solve landed
  on 3.10.9 and got PySide6. Same file, different answer, which is why this
  reproduced in CI and not on a developer machine.

  `environment.yml` now asks for **`matplotlib-base`**, which never depends on a Qt
  binding. Nothing here needed the metapackage: scanpy, squidpy and
  matplotlib-scalebar all depend on `matplotlib-base`, and this app supplies its own
  binding. The solved environment now contains PyQt6 alone — no PyQt5, no PySide6,
  no Qt5 stack at all, which also removes ~70 MB of duplicate Qt.

- **The Qt backend is now stated rather than inherited.** `xenium_viewer/__init__.py`
  sets `QT_API=pyqt6` via `setdefault`, and only when PyQt6 is importable — so an
  explicit `QT_API` still wins and an environment with only PySide6 is left alone.
  Fixing `environment.yml` fixes a *fresh* env; this covers an existing one that a
  user upgrades in place, where whatever binding is already installed would
  otherwise keep winning silently.

  `tests/test_qt_backend.py` pins all of it, including a source guard that fails if
  a bare `matplotlib` dependency comes back. Verified in both environments:
  PyQt5 1047 passed / 25 skipped, PyQt6 1048 passed / 24 skipped (the differing
  skip is the backend-specific test in each).

#### Changed
- **Migrated the Qt backend from PyQt5 to PyQt6** (issue #15). napari deprecated the
  PyQt5 backend for removal in autumn 2026 and warns about it at every startup, along
  with a second warning that system theme detection needs Qt6.

  `environment.yml` now asks for `pyqt6`, `pyproject.toml` for `PyQt6`, and all 98
  unscoped Qt enum sites plus 8 `.exec_()` calls are in their Qt6 form. The
  replacement map was derived from the live bindings rather than written by hand — for
  each `Klass.MEMBER` in the tree, look up which nested enum of `Klass` defines MEMBER
  — and every scoped form was checked to evaluate to the same integer under PyQt5
  before anything was rewritten, so the edits were a no-op on the old backend.

- **`napari` is now pinned `>=0.8`.** With `pyqt6` swapped in, the solver quietly
  resolved napari to **0.7.0** — not a Qt conflict, but the solver easing the unrelated
  `zarr>=3.0,<3.2` pin by walking napari back a minor version. `napari=0.8.0` + `pyqt6`
  solves cleanly on its own. A silent downgrade of the thing this application *is* was
  worth closing off.

- **CI asserts the resolved backend** before running the suite — that qtpy reports
  PyQt6 and that napari is still >=0.8. Both regressions are silent, both happened
  once during this migration, and a green suite on the wrong backend is exactly what
  the migration was for.

#### Notes
`docs/pyqt6-migration.md` was rewritten from a plan into a record, because the plan was
wrong in both directions. It listed 8 unscoped enums (there were 98) and treated those
edits as the work — but qtpy's `enums_compat.promote_enums()` and its `exec_` aliases
mean the *pre-migration* code would have run unchanged under PyQt6. The edits were
hygiene; the dependency solve was the real task, and the plan did not mention it.

Two packaging facts cost the most time and are recorded so they cannot cost it again:
conda-forge ships PyQt6 as **`pyqt6`**, not as a 6.x of `pyqt` (which stops at 5.15, so
`pyqt=6` fails as though PyQt6 were unavailable); and **`QT_API=pyside6 pytest` never
smoke-tested strict enums** as previously claimed — PySide6 runs in forgiveness mode and
accepts every unscoped form.

Verified: full suite green on PyQt5 (1045 passed, 24 skipped) before the pin flip, and
again on a real PyQt6 6.8.1 / Qt 6.8.1 / napari 0.8.0 environment after it.

### 2026-08-24 — scale bar in micrometres

#### Added
- **The canvas scale bar now reads in micrometres**, switching to millimetres as you
  zoom out (issue #25). It is on from the moment a dataset loads, and the conversion
  comes from `pixel_size` in `experiment.xenium` — nothing is assumed.

  napari 0.8 has no `scale_bar.unit`; units live on the layers (`layer.units`,
  pint-backed) and the *magnitude* has to live in `layer.scale`. The napari-0.5-era
  shortcut of a scaled unit string is worse than unavailable — `layer.units =
  "0.2125 um"` is **accepted and silently discards the magnitude**, leaving one pixel
  labelled as one micrometre: a wrong scale bar with no error anywhere. A test pins
  that behaviour so nobody reintroduces it.

  So every layer is given `scale = (pixel_size, pixel_size)`, which makes napari's
  world coordinates micrometres. Applied through the viewer's `layers.inserted` event
  rather than at each `add_*` call site: there are more than twenty of those across
  eight modules, and one missed would leave that layer in pixels — *misplaced*
  relative to everything else, not merely mislabelled.

#### Fixed
- **Registered overlays would have silently moved.** napari composes
  `world = affine(scale(data))` — the affine is applied *after* the scale, so its
  translation is in world units, while every affine this codebase stores is in Xenium
  pixels (registration fits it from landmark pixels; `adata_persistence` writes it to
  the zarr in that frame; the crop export composes with it there). Switching the world
  to micrometres without touching them would shift every H&E, ARMS section, external
  image and patch overlay by ~4.7×, with nothing raised and nothing logged.

  `utils/units.py` is that boundary and the only place the two frames meet: stored
  affines stay in pixels, napari layer affines are in world units. The conversion is
  a similarity conjugation, which reduces to scaling the translation column — a
  rotation is the same rotation in any unit, an offset of 100 px is 21.25 µm.

  Copying an affine between layers (`utils/affine_linking.py`, `img_lm.affine =
  lyr.affine`) is unit-agnostic and deliberately left alone; converting there would
  be a double-scaling bug.

- **The minimap read the camera as if world units were pixels.** Its viewport
  rectangle and click-to-navigate both convert between `camera.center` and the
  morphology shape, so both would have been off by ~4.7× — landing inside the tissue,
  plausibly, which is what would have made it easy to miss. It now works in world
  units throughout, and defaults to `pixel_size=1.0` so an unscaled viewer is
  unchanged.

#### Notes
Layer **data** is still in Xenium pixels everywhere, so nothing that reads
`layer.data` — the crop export, the ROI tab, the ARMS tile ingest — is affected.

`tests/test_units.py` (16 tests) asserts *placement*, not labelling: a test that only
checked "the scale bar says µm" would pass with every overlay in the wrong place. It
also pins the upstream premise (`napari` applies affine after scale), since the whole
conversion rests on a behaviour we do not control. The tests build bare
`napari.layers.Image` objects rather than a `Viewer`, so they need no OpenGL context
and run in CI.

### 2026-08-24 — readthedocs

#### Added
- **The docs site builds** (issue #16). `mkdocs.yml` had shipped a complete 44-page
  nav for months and `site_url` already pointed at `xenium-viewer.readthedocs.io` —
  but the site had **never been built**, and would not have worked if it had.

  `docs/` is authored as GitHub Wiki source, where links are extensionless
  (`[Clustering](Tab-Clustering)`). The wiki resolves those; mkdocs reads them as
  paths to files that do not exist. Measured: **158 dead cross-links** — every
  cross-reference on the rendered site.

  `mkdocs_hooks.py` rewrites them at build time, so `docs/` keeps *one* link
  convention rather than gaining a second one that the wiki publisher does not check.
  The hook lives at the repo root, not in `docs/`, because the publish script selects
  wiki pages by filename convention and a stray module there has leaked to the public
  wiki before.

- `.readthedocs.yaml` and `requirements-docs.txt`. The build needs no conda env, no
  scanpy and no Qt — `Analysis-Templates.md` and `API-Reference.md` are generated and
  checked in — and that separation is deliberate: a docs build that needs the
  application is a docs build that breaks whenever the application does.

- **CI builds the site with `--strict`** on every push, in the *lint* job.
  `mkdocs.yml` also raises link validation to `warn`, without which a dead link is an
  INFO line and the build still "passes" — which is exactly how this shipped broken.

- `exclude_docs` keeps internal planning notes and the wiki's `_Sidebar.md` out of
  the published site, using the same lower-case/Title-Case convention
  `scripts/push_to_wiki.sh` already uses.

#### Notes
**Connecting the Read the Docs project needs the repository to be public** — RTD
Community does not build private repos. Everything above is done and verified; the
connection is the one step that cannot be. `docs/readthedocs-setup.md` records the
remaining steps. It unblocks when the repository is made public, rather than needing
anything bought.

`tests/test_docs_links.py` gained 11 tests covering the rewrite rule, including that
every bare link in `docs/` rewrites onto a page that actually exists — derived from
the files present, not from a list.

### 2026-08-19 — H&E image memory

#### Fixed
- **Loading a second H&E over one restored from the cache displayed the new image but
  never persisted it** — the next launch brought back the old one, and any landmarks
  placed in between were saved against an image the cache did not hold. The same defect
  hit **ARMS H&E**, **external images**, and **custom segmentation labels**.

  `safe_write_element` opens with `_assert_not_dask_backed`, which refuses to rename an
  element's directory into `.xv_trash/` while something still reads it. That guard is
  right — a dask graph resolves lazily *by path*, so a surviving reader would go on to
  read the *new* element's bytes rather than fail cleanly. But it inspects only
  `sdata._gen_spatial_element_values()` — **not napari layers, not open file handles** —
  so a caller that has already torn down every reader still tripped it, purely because
  `sdata.images["he_image"]` was the element it was about to replace. The guard runs
  before staging, so the write was a clean no-op: the symptom was silence.

  `safe_write_element(..., replace_backed=True)` is how a caller says it has done the
  teardown. It drops the in-memory binding (`_unbind_backed`, built on the existing
  `_drop_in_memory`) and **restores it if the write fails before it commits** — not
  after, since a failure in `consolidate()` happens once the new data is already live.
  It still refuses when the element being written *is* the store-backed one, so
  `test_refuses_to_move_a_dask_backed_element` holds with the flag set: "I have torn down
  every reader" cannot be true of the value being written.

  Passed at the three sites where the old napari layer is provably gone by the time the
  write runs — `_save_he_to_sdata`, `_save_arms_he_to_sdata`, and
  `save_external_image_to_sdata`. The last of those had **hand-rolled the same unbind**
  (`del ctx.sdata[element_name]` under a bare `except: pass`, with a comment admitting it
  did not know whether it worked); it now gets rollback and stops swallowing.

- **Loading a new custom segmentation over a cached one had the identical bug.**
  `save_custom_seg_to_sdata` takes `replace_backed` and only the safe caller passes it:
  `_on_done` writes a *different* segmentation after `_apply_custom_segmentation` has
  removed the old labels layer. **Tools → Update SpatialData deliberately does not** — it
  writes `ctx.cell_labels_layer.data` back while that layer is on screen, so the guard is
  correct to refuse there. Pinned by a test that parses both call sites and requires
  exactly one to opt in.

#### Fixed (memory)
- **Loading an H&E built RAM fast enough to kill the session.** Measured on a
  16384x12288 RGB TIFF: load + write to the cache peaked at **2.61 GB**, and the
  session restore of the same image at **1.71 GB** — for a 604 MB image. On a real
  slide that scales to tens of GB, and the failure mode is a killed terminal with no
  traceback (`journalctl -u systemd-oomd`), not an error.

  **This was never a regression in the H&E code.** `utils/registration.py` — the only
  code that reads an H&E from disk — had not changed since `564179c` (2026-05-07), and
  that commit only swapped `da.from_zarr(store)` for `zarr.open` + `da.from_array`
  to survive zarr 3; dask's `from_zarr` *is* `from_array(z, z.chunks, name=…)`, so it
  is inert for memory. `_save_he_to_sdata`'s `np.asarray(pyramid[0])` had been there
  verbatim since the March decomposition. What changed is that the 2026-08-17 memory
  work — the tiled read in `raster_io`, `loader._reopen_written_cache`, and
  `sdata_write`'s per-element caps — all stop at the loader/morphology boundary, and
  every H&E entry point sits outside it. Morphology got quiet; H&E did not.

  Two sites, both eager, both fixed:

  - **The cache write** (`tabs/tab_he_registration.py`, `tabs/tab_arms.py`) did
    `np.asarray(pyramid[0])` — the full-resolution slide as dense numpy, on the **Qt
    main thread** — then a second full copy via a numpy `.astype(np.uint8)`, then
    handed `Image2DModel.parse` a dense base whose four `scale_factors` levels
    spatialdata computes in a single `da.compute` through float64. Now
    `registration.parse_rgb_image_for_store()`, shared by both tabs and built on the
    same `da.map_blocks` detach that `adata_persistence.save_external_image_to_sdata`
    already used. **2.61 GB → 0.91 GB**, and the written element is **byte-identical**
    to what the eager version produced (asserted per level, both tiled and untiled).
  - **The session restore** (`_load_he_from_sdata`) called `.compute()` on *every*
    level including `scale0`, then transposed each into a second array — on every
    launch of a dataset with an `he_image`. It now hands napari dask, which is what
    the line-for-line identical ARMS code has done since `9cad210` (2026-03-11); the
    H&E copy beside it was missed for five months. **1.71 GB → 0.49 GB.**

- **`app._warn_if_pyramid_is_not_stored` now covers every image element**, not just
  `morphology_focus`. `he_image`, `arms_he_image` and the `ext_*` images reach napari
  the same way and cost the same when their levels are a computation; scoping the
  warning to one element is why an H&E could kill a session without the viewer having
  said anything first. The "one task per chunk means it was read, not built" test is
  now `raster_io.level_is_computed()`, shared rather than inline.

- **The H&E/ARMS/external-image tabs log what they loaded**
  (`registration.describe_pyramid`): level count, base shape and dtype, chunking,
  whether the levels are the file's own or a chain, and — the fact that actually
  predicts the memory — whether the base is **tiled**. "The viewer died loading an
  H&E" carried no information about the file that killed it.

#### Investigated, no change
- **The `da.coarsen` chain `load_he_pyramid` synthesises for a TIFF with no internal
  pyramid is *not* the problem, and materialising each level as it is built makes
  things worse.** This looked like the morphology bug wearing a different hat — napari
  draws the smallest level first (`_data_level = len(data) - 1`), so touching it walks
  the chain — and it was implemented before being measured. It is not the same:
  morphology's chain stood on **whole-page** dask chunks (5.93 GB each, per
  `utils/raster_io.py`), whereas a coarsen over a *tile-chunked* base streams. Measured
  peak, load + write, same 604 MB image:

  | | tiled | untiled (one strip) |
  |---|---|---|
  | before | 2.61 GB | 2.55 GB |
  | lazy write only (shipped) | **0.91 GB** | 2.69 GB |
  | + materialise each built level | 2.14 GB | 2.88 GB |

  So the materialisation was reverted and the chain left lazy. For an **untiled** H&E
  nothing here helps: the floor is tifffile decoding the whole plane (1.34 GB of the
  2.67 GB peak is the decode alone, before any pyramid exists), and decoding once and
  re-exposing it as tiles was measured too — 2.67 GB, no better. That is the same
  conclusion `--no-cache` morphology reached, so it warns rather than pretending
  otherwise.

### 2026-08-18 — dataset rename

#### Added
- **`xenium-rename-dataset`** — a console script that renames or moves a dataset
  directory and rewrites the absolute paths recorded inside it
  (`src/xenium_viewer/scripts/rename_dataset.py`).

  Nothing stores a dataset *name*: not the AnnData table, not the SpatialData store,
  not the cache manifest — the folder name reaches only the napari window title.
  Cache freshness is a content hash of `experiment.xenium`, so moving a dataset never
  triggers a rebuild, and every cache/sidecar location is re-derived from `data_path`
  at launch. Verified directly: every `zarr.json` / `.zattrs` / `.zgroup` /
  `.zmetadata` under a real `sdata_cached.zarr` was grepped for an absolute path and
  none contains one.

  What a bare `mv` *does* break is the paths recorded in the provenance graph — the
  `preamble` node's `data_path = Path(r"…")` and each `clustering:<key>` node's
  `read_csv` of `analysis/clustering/<key>/clusters.csv` — plus the session's
  `he_path` / `arms_*_path` attrs. And it does not stay quiet: `app.py` re-emits the
  preamble for the current `data_path` on every launch, so `ProvGraph.upsert` flags
  **every transitive descendant stale** and the first launch after a manual rename
  marks the whole notebook ⚠ despite nothing being recomputed. Repairing the graph
  before that launch is the point.

  ```
  xenium-rename-dataset /data/old_name new_name      # rename in place
  xenium-rename-dataset /data/old_name /other/place  # move
  xenium-rename-dataset /data/new_name --repair      # fix one moved by hand
  xenium-rename-dataset ... --dry-run                # report, write nothing
  ```

  `--repair` infers the previous path from the recorded preamble, so the common case
  needs no extra argument. Substitution is **path-prefix only**, which is what leaves
  an H&E or GeoJSON living outside the dataset directory untouched — that file did not
  move. `analysis.py` and `analysis_notebook.ipynb` are regenerated from the repaired
  graph (the same recipe `app.py` uses on session restore) rather than string-patched,
  and only when they already exist.

  Store writes go through `safe_group_update` / `atomic_json` only — a source guard
  parses the module and fails on a bare `zarr.open`, `shutil.move`, or any `unlink`.
  `is_dataset_dir` reuses `loader.has_raw_xenium_source` (issue #17) rather than
  re-testing `cells.zarr.zip`, so a Crop Dataset export — whose zarr *is* the data —
  is recognised by the one definition rather than a second copy of it.

  It refuses rather than half-acting: an existing destination, a running CopyKAT job
  (`plots/copykat_RUNNING.txt`), unreadable root metadata, a cross-device move (no
  `shutil.move` fallback — a half-copied multi-gigabyte zarr is worse than a clear
  message), or a directory that is not a dataset. Interrupted safe writes are
  recovered before the move.

#### Fixed
- **The rename tool's derived-file writes now swap an inode instead of truncating
  one.** Found by probing the real tool against a `cp -al` snapshot of a live dataset:
  `analysis.py` and the notebook were the only writes not going through
  `atomic_json`, so regenerating them wrote *through* the hardlink and silently edited
  the snapshot too. Both now use a temp file plus `os.replace`.
  `test_derived_outputs_are_replaced_not_written_in_place` states the property —
  a hardlink taken before the rename still holds the old bytes after it — and fails
  against the previous version.
### 2026-08-17 — squidpy 1.9 readiness

#### Changed
- **`sq.gr.spatial_neighbors` → `sq.gr.spatial_neighbors_knn`** (issue #19), before
  squidpy 1.9 removes the old name. Two call sites that had to change together —
  `utils/spatial_analysis.py::compute_spatial_neighbors` (what the GUI runs) and
  `step_templates/builtin/spatial_neighbors.tmpl` (what the exported notebook replays) —
  plus the `coord_type='generic'` prose in `scripts/generate_docs.py` and the regenerated
  `docs/Analysis-Templates.md`, which carries the template body verbatim.
  `compute_spatial_neighbors` keeps its signature, so `tab_annot_nhood` is untouched.

  **No result changes and nothing to recompute.** Re-measured on the installed squidpy
  1.8.2: `spatial_neighbors_knn(n_neighs=6)` and
  `spatial_neighbors(coord_type='generic', n_neighs=6)` write byte-identical
  `spatial_connectivities` and `spatial_distances` and the same
  `uns['spatial_neighbors']` — the new function still records `coord_type: 'generic'`
  there, so even code reading that key is unaffected. `coord_type='generic'` with
  `n_neighs=k` *is* k-nearest-neighbours.

  A notebook **already exported to disk** keeps the old call, because a provenance node
  stores the code verbatim from when it ran. This self-heals rather than needing a
  migration: the next session that opens Neighbourhood Enrichment or Ligand-Receptor
  re-runs `ensure_spatial_neighbors`, `upsert` revises the node, and its descendants are
  flagged stale — which is what a changed upstream step is supposed to do.

#### Corrected
- **Raised the squidpy floor to `>=1.8.2`** in `environment.yml` and `pyproject.toml`.
  The first cut of this change left it at `>=1.8` on the strength of
  `spatial_neighbors`'s own `.. deprecated:: 1.7.0` directive, read as "the replacements
  exist from 1.7.0". They do not. Checked against the wheels,
  `spatial_neighbors_knn` is absent from **1.7.0, 1.8.0 and 1.8.1** and first appears in
  **1.8.2** — so the old pin permitted two versions on which every Nhood and
  Ligand-Receptor run would raise `AttributeError`. The lesson is narrow and worth
  keeping: a `deprecated::` directive dates the *deprecation*, not the *replacement*.
  `test_the_squidpy_floor_is_a_version_that_has_spatial_neighbors_knn` now checks both
  the installed module and the declared floors, and fails on the old pin.

- The 2026-08-17 tracking entry below, and `docs/squidpy-spatial-neighbors-migration.md`,
  both said **three** tabs broke together. Co-occurrence was never a dependent:
  `tab_co_occurrence` declares only `deps=[clustering:…]` and `spatial.cooccur.tmpl`
  calls `sq.gr.co_occurrence`, which computes its own radii and needs no `obsp` graph.
  Two tabs, not three. Also: the `cnv` extra contains no squidpy, and
  `tests/test_notebook_replay.py` *does* cover this path — it replays a
  `spatial_neighbors` node in a clean kernel, which is the gate the doc said did not
  exist.
### 2026-08-17 — cache-only datasets

#### Fixed
- **A dataset with no raw Xenium output is never told to rebuild its cache**
  (issue #17). A Crop Dataset export ships `experiment.xenium`,
  `sdata_cached.zarr/`, `transcripts.parquet` and `transcript_cache/` — and
  nothing `spatialdata_io.xenium` can read. It also shipped **no
  `.xv_manifest.json`**, so freshness fell back to comparing `experiment.xenium`'s
  mtime against the cache directory's, which a `cp -r`, an unzip or a cloud-sync
  restore reorders. Three separate branches then attempted a rebuild that cannot
  happen, and **two of them renamed the only copy of the data aside first**:
  - `tab_cache._on_rebuild` — the reported case. It checked only that the cache
    existed and that there was disk space, moved the store to
    `sdata_cached_prev_<stamp>.zarr`, and told the user to restart. The next
    launch died in h5py on a `cell_feature_matrix.h5` they never had, and the
    dataset stayed unopenable until someone renamed the folder back by hand.
  - the silent mtime-staleness branch (`_stale_preference` returning `None`),
    which left the cache in place but raised on every launch, permanently.
  - `_ask_corrupt_cache`, whose only non-quit answers both move the cache aside.

  The fix makes "this has no raw source" a fact rather than an inference.
  `loader.has_raw_xenium_source()` is now the single definition — it replaces an
  inline `cells.zarr.zip` check that existed but guarded only the `--no-cache`
  override, which is why every rebuild path walked past it. It is deliberately
  conservative: `True` unless *none* of `cells.zarr.zip`,
  `cell_feature_matrix.h5`, `cell_feature_matrix/`, `morphology_focus/` is
  present, because partially-missing raw output is broken raw output and should
  keep raising the error that says so.

  For such a dataset the loader skips the staleness question entirely (the cache
  *is* the source of record), refuses with a new `NoRawSourceError` naming the
  recovery routes instead of offering a rebuild when the cache will not open, and
  carries a precondition immediately before `spatialdata_io.xenium` so any branch
  added later is safe by construction. Nothing is ever moved — asserted with a
  `shutil.move` spy rather than a source grep, since what matters is that no
  rename happens.
- **Crop Dataset exports now stamp their manifest at export time**, with
  `cache_only: true` and `derived_from`. That alone closes the mtime trigger for
  new exports; the predicate covers the ones already on disk, and stamps them on
  their next successful load.
- **Force Rebuild is disabled, with the reason as its tooltip**, on a cache-only
  dataset, and re-guarded inside the callback — `_set_busy(False)` re-enables
  every mutating button when a worker finishes, and this is the one that must
  stay off. The precondition lives in a pure `_rebuild_blocked_reason()` so it is
  testable; `_on_rebuild` itself opens a modal. The Cache tab header and
  `xenium-build-cache --check` both report the dataset as cache-only, and
  `--check` reports freshness as *not applicable* rather than stale — condemning
  the cache would be advising a rebuild that cannot happen.
- **Tools → Dataset stops calling a crop export's own files "original 10x output
  — the viewer never modifies it".** The viewer wrote them. They stay
  non-deletable (they are the only copy), but that wording was the same wrong
  mental model that let Force Rebuild loose on these datasets. A test docstring
  in `test_loader_policy.py` carried the mirror-image error — "a Crop Dataset
  export has no `experiment.xenium`" — and is corrected too.

### 2026-08-17 — headless cache build

#### Added
- **`xenium-build-cache` console script** — builds `sdata_cached.zarr/` for a dataset
  without starting the GUI, so the cold read (tens of minutes, tens of GB) can happen
  over ssh or overnight. The capability already existed as `loader.main()` but was never
  registered in `[project.scripts]`, never documented, and its `--help` still advertised
  `scripts/01_load_sdata.py`, deleted in the March 2026 refactor. `xenium-preprocess`
  only ever built the *transcript* cache, which is what made the zarr one look like
  GUI-only work.

  Flags: `--on-stale {ask,keep,rebuild,restore}`, `--no-pyramid`, `--n-jobs`, `--check`.
- **`--check`** reports a cache's freshness, integrity and user-generated contents and
  exits without building anything, non-zero if it is missing, stale or does not verify.
  It reuses `cache_repair.verify` and `_detect_user_data`, both filesystem-level, so it
  works on a store too broken for zarr to open — which is when it gets run.

#### Changed
- **`load_sdata(on_stale=)`** answers the stale-cache question in advance instead of
  prompting. This closes a real gap rather than relaxing a safety rule: with no dialog
  available `_ask_rebuild_preference` returns `'keep'` and `_ask_corrupt_cache` raises,
  deliberately — but that left a terminal user with no way to say yes at all, short of
  moving the cache aside by hand. The default (`None`) is unchanged, and the GUI passes
  it. `'keep'` still does not satisfy `_ask_corrupt_cache`: a cache that will not open
  cannot be kept, and pretending otherwise would hide the breakage.
- The stale-cache decision moved into `loader._stale_preference`, a near-pure function.
  It checks `on_stale` **before** the branch chain, because two of those branches (the
  pre-manifest rebuild, and a certain-stale cache holding nothing) return without ever
  consulting a dialog — a preset threaded only through the ask-helper would have been
  silently ignored there. Extracting it also made every branch testable without data.

### 2026-08-17 — documentation

#### Tracked, not fixed
- **`sq.gr.spatial_neighbors` is deprecated in squidpy 1.8.2 and removed in 1.9.0** —
  `docs/squidpy-spatial-neighbors-migration.md`, flagged at the top of `CLAUDE.md` so it is
  seen before the next dependency bump. Found while executing the API-Reference snippets.
  It matters more than a leaf deprecation because the neighbour graph is a **dependency**:
  Neighbourhood Enrichment, Co-occurrence and Ligand-Receptor break together, and so does
  every **already-exported notebook**, whose recorded `spatial_neighbors` cell calls the
  removed function — which is precisely the promise the provenance graph exists to keep.
  `squidpy>=1.8` is pinned with no upper bound in `environment.yml` and `pyproject.toml`,
  so an ordinary `conda env update` is enough to trigger it.

  The migration is nonetheless **low-risk, and measured rather than assumed**: on 500
  random points under squidpy 1.8.2, `spatial_neighbors_knn(n_neighs=6)` writes the same
  `obsp`/`uns` keys and a **byte-identical** `spatial_connectivities` and
  `spatial_distances` to `spatial_neighbors(coord_type='generic', n_neighs=6)`, and
  `nhood_enrichment` runs unchanged against it. `coord_type='generic'` with `n_neighs=k`
  *is* k-nearest-neighbours. So it is a rename across two call sites, not a change of
  results — no saved analysis is invalidated. The doc carries the equivalence check to
  re-run against whatever squidpy version is current at the time.

#### Added
- **[Analysis Templates](docs/Analysis-Templates.md)** — a catalogue of all fourteen analysis
  templates: what each computes, which tab runs it, its full contract, and its **complete
  default source, block by block**. The templates *are* the analysis, and until now none of
  them could be read without launching the viewer and clicking through the Templates tab.
- **[API Reference](docs/API-Reference.md)** — the ~50 functions that are usable outside the
  GUI, with live signatures, a runnable snippet each, and the distinction that matters:
  `run_rank_genes(adata_norm, groupby, …)` takes a plain `AnnData` and works in any notebook,
  while `crop_and_export(ctx, …)` needs a `ViewerContext` and is documented as *not* API. No
  page in `docs/` previously contained the string `from xenium_viewer`.
- `scripts/generate_docs.py` generates both pages — the template bodies from the registry,
  the signatures from the live objects — with the prose held in the script so regenerating
  never destroys writing. **`tests/test_generated_docs.py` fails if a checked-in page differs
  from what the generator produces**, if a shipped template has no prose, or if a documented
  function no longer exists. Code in the docs cannot go stale without CI noticing.

#### Fixed
- **The four missing screenshots** (CNV, Dataset, Cache, Templates) are captured; those pages
  had carried unfilled placeholders since the tabs were added.
- **Three documented calls were wrong, and executing the snippets is what found them** —
  which is the argument for the snippets being executable at all. `load_clusterings` takes a
  *dataset path* and returns 10x's own `analysis/clustering/` results, not an `AnnData` and
  not the clusterings you computed; `TranscriptLoader.cached_genes` is a property, not a
  method; and the rank-genes example grouped by a `'leiden'` column that no dataset actually
  has, since the viewer names its clusterings `clustering_*`.
- **A generated page had this machine's paths in it.** `TranscriptLoader`'s `cache_dir` and
  `parquet_path` defaults resolve against the working directory at import time, so the
  extracted signature contained an absolute path that is wrong for every reader. Shortened to
  `.../name`, with a test that fails on any absolute path in a generated page.
- `Installation.md` gains a **Memory** section and a "the viewer vanished and took the
  terminal with it" entry naming `journalctl -u systemd-oomd` — the failure looks nothing
  like an out-of-memory error, because oomd fires on memory *pressure* and kills the whole
  cgroup scope. It also no longer recommends `--no-cache` for read-only filesystems without
  saying what it costs on a full slide.
- **Two docs tests that could not fail**, both found while adding a third:
  `test_internal_links_resolve` split the anchor off a `Page#anchor` link and checked only
  the page, so a link to a heading that never existed passed — `test_link_anchors_resolve`
  now checks the heading, across 30 such links. And the new screenshot-index guard matched
  nothing on its first version (0 of 26 entries) because it compared page names that still
  had `.md` on them; it now asserts a count as well as the property, since a guard that
  silently checks nothing is worse than no guard.
- `Tab-Crop-Dataset.md`, `Tab-Cache.md`, `Tutorial-Getting-Started.md` and
  `Tutorial-Clustering.md` describe current behaviour: bounded crop memory, per-element
  memory reporting during a rebuild, per-clustering ranked genes, and the two clustering
  preprocessing checkboxes the tutorial never mentioned.

### 2026-08-17

#### Performance
- **Building `sdata_cached.zarr` needed ~24 GB before it could write a single byte**, and
  the cause was the *chunking of the read*, not anything in the write. `spatialdata_io`
  reads `morphology_focus` through `dask_image.imread`, which yields **one dask chunk per
  full channel page** — measured on a 57887×51217 4-channel slide, `chunksize
  (1, 57887, 51217)` = **5.93 GB per chunk**. `Image2DModel.parse(chunks=(1, 4096, 4096))`
  then rechunks that into ~195 tiles per channel, but dask cannot produce *any* tile
  without decoding a whole page, and must hold it until every consumer is done — including
  the chained pyramid levels above it, whose top two levels are single chunks depending on
  the entire level below. All four channel pages therefore land in memory at once.

  The file was never the problem: it is tiled 1024×1024 and readable a tile at a time
  through `tifffile.imread(..., aszarr=True)` (opens in 0.85 s; a 1024² sub-region reads
  from the NAS in 0.36 s). `utils/raster_io.py` reads morphology_focus that way and the
  loader re-parses it, leaving the pyramid computation to spatialdata exactly as before —
  **so every written byte is unchanged**; only the chunking of the read differs.

  **Measured on the dataset that prompted this, both runs under an `RLIMIT_AS` cap so an
  over-budget run fails instead of taking the machine down:**

  | | peak RSS | outcome |
  |---|---|---|
  | before | **41.2 GB** | died at a 60 GB cap: `Unable to allocate 5.52 GiB for an array with shape (1, 57887, 51217)` |
  | after | **8.5 / 9.2 GB** (VmHWM 8.97 / 9.26 GB) | full 14 GB store written under a 24 GB cap, in 473 s / 491 s |

  (Two independent complete runs; the second was a re-measurement on the final code.)

  And the reported symptom — memory climbing element after element — is gone; RSS is now
  flat across the whole write, with the peak set once during the image and never revisited:

  ```
  [mem] wrote morphology_focus: rss=4.74GB peak=8.81GB
  [mem] wrote nucleus_labels:   rss=4.93GB peak=9.26GB
  [mem] wrote cell_labels:      rss=4.97GB peak=9.26GB
  [mem] wrote transcripts:      rss=5.79GB peak=9.26GB
  [mem] wrote table:            rss=5.82GB peak=9.26GB
  ```

  The resulting store was diffed against a cache built by the old code on the same dataset:
  identical level counts, shapes, dtypes, chunk grids and sampled data for **all six image
  pyramid levels and both label pyramids**, identical table `X`/`obs`/`var_names`, identical
  transcript row count (218,457,525) and partition count. The one difference is the *order*
  of the transcripts columns (same set, same data) — `PointsModel.parse` iterates
  `set(data.columns)` and documents its column order as "not guaranteed", so it varies run
  to run upstream of anything here, and nothing in the viewer reads those columns by
  position.

  The TIFF's own pyramid levels 1+ are deliberately **not** reused, even though they exist,
  are tiled, and match `scale_factors=[2, 2, 2, 2, 2]` shape for shape. 10x downsamples with
  a different filter, and measured against the real slide the difference is not a rounding
  error: mean |diff| of **8.63** against `coarsen().mean()` on a level whose mean value is
  57.6 (≈15% of every pixel drawn when zoomed out), 5.98 against decimation, 14.74 against
  max-pooling. Reading only level 0 gets the same 24 GB back with no change to the output at
  all. The reader **declines** — leaving spatialdata_io's element untouched — if there is no
  readable OME-TIFF, or if its full-resolution level disagrees with what spatialdata_io
  produced in shape or dtype, so older Xenium layouts behave exactly as before.

- **dask's concurrency is now capped per element type rather than per write.** The 2-worker
  cap introduced for the crop path applied to whole writes, which is both too much and too
  little:
  - **Images do not need it.** Measured on the real slide, capping them moved the peak only
    8.71 GB → 8.31 GB while making the raster write several times slower. Their chunking was
    the problem, not their concurrency.
  - **Labels do**, for a reason images do not share. Both pyramids are lazy, but labels go
    through `spatialdata.models.pyramids_utils.to_multiscale`, which does an explicit
    `prev.astype(float)` and **two rechunks per level** (either side of
    `ome_zarr.dask_utils.resize`'s `map_blocks`); a rechunk is a many-to-many gather and the
    float64 intermediates are 134 MB a chunk. Uncapped, it was the *label* write — not the
    30× larger image write — that died: `Unable to allocate 177. MiB for an array with
    shape (3620, 6404) and data type float64`. Labels are only ~277 MB on disk, so the cap
    costs a couple of minutes.

  `sdata_write.ELEMENT_WORKERS` declares the exceptions; `crop_export` and the cache build
  share one constant so the two write paths cannot drift apart.
- The morphology re-read is **not** gated on `load_sdata(build_pyramid=...)`, because the
  whole-page read happens either way: `spatialdata_io` fills in its own default
  `scale_factors` whenever the caller passes none, so `build_pyramid=False` never actually
  meant "single scale". For the same reason the replacement's `scale_factors` are derived
  from the element being replaced rather than from a constant.

#### Added
- `utils/sdata_write.write_sdata()` — writes a store **one element at a time** instead of
  through `SpatialData.write()`'s internal loop, so there is a point at which progress can
  be reported, dask concurrency chosen per element type, and memory released. Uses only
  public spatialdata API (an empty `SpatialData` writes the store shell and performs the
  path-safety validation; `write_element` does the rest), and a test asserts it produces a
  store byte-identical to `SpatialData.write()`.
- `utils/mem_probe.py` — RSS/peak-RSS probes, a `malloc_trim` wrapper and a reader for
  napari's global dask cache. The cache build now logs a memory line per element, to the
  terminal and the dataset log, because "how much is this using" is a question people ask
  *while* a tens-of-minutes build is happening.

#### Fixed
- **Loading a dataset for the first time killed the terminal while the GUI was adding
  `morphology_focus`** — build the cache and display it in one session and it died;
  restart, load the same dataset, and it was fine. Both halves have the same cause:
  `load_sdata` wrote the cache and then returned the **in-memory** object it had built,
  not the store it had just written. Setting `sdata.path` pointed later *writes* at the
  cache but rebound no array, so every layer was still backed by the build-time dask
  graph, in which the pyramid levels above `scale0` exist only as a chained
  `coarsen().mean()`.

  That matters because **napari displays the smallest level first**
  (`_scalar_field.py`: `_data_level = len(data) - 1`), and so does the thumbnail, and so
  does the coarse-alignment thumbnail in `_populate_viewer`. On a full slide `scale5` is
  20 MiB of pixels standing on **13,444 tasks over all 780 `scale0` chunks**, promoted to
  float64 on the way — three layers over. Touching the smallest level materialised the
  largest one.

  `load_sdata` now re-opens the cache it just wrote (`_reopen_written_cache`) and returns
  that, so the first session holds exactly what every later session holds. Measured on
  `Run4/…PCA030_PCA031…`, layer-building only, under a 40 GB `RLIMIT_AS` cap:

  | | peak RSS | outcome |
  |---|---|---|
  | before (build-time object, uncapped) | **23.1 GB** | died in 51 s: `MemoryError: Unable to allocate 32.0 MiB for an array with shape (1, 2048, 2048) and data type float64` — a 32 MiB coarsen intermediate failing is what address-space exhaustion looks like |
  | after (re-opened from cache) | **1.70 GB** | 9 s, complete |

  The re-open is lazy and costs a few seconds. It also settles `sdata.path` without a
  separate assignment, and fixes the same latent problem for `points/transcripts` and the
  table, which were likewise left pointing at the raw output after a build.

  **What was actually killing the terminal is `systemd-oomd`**, and it is worth knowing
  because the numbers look wrong otherwise: it fires on sustained memory *pressure*
  (50% PSI for 20 s), not on absolute usage, which is why RSS was never seen near the
  125 GB the box has — and it kills the **cgroup scope**, i.e. the terminal and every tab
  in it, rather than the offending process. `journalctl -u systemd-oomd` names the scope
  and the pressure figure outright, and is the first thing to check the next time a
  session disappears without a traceback. Left as it is: it is the reason a runaway
  allocation costs a terminal instead of the desktop.
- **`--no-cache` now says what it is about to do**, because the fix above cannot help it.
  With no cache to read, the lower pyramid levels have to be computed, and that
  computation does not stream: each level is rechunked either side of the coarsen, which
  is a many-to-many gather, so a whole level of float64 must be live at once — ~24 GB for
  `scale1` on a full slide. Capping dask's concurrency was tried first and **does not
  work**, which is worth recording because it sounds like it should: same dataset, same
  40 GB cap, **23.1 GB uncapped vs 25.2 GB at 4 workers, both dead**. The memory is a
  property of the graph's shape, not of how many workers walk it. So `_populate_viewer`
  detects a pyramid whose levels are computed rather than stored and warns, naming
  `journalctl -u systemd-oomd` — turning a session that vanishes without a traceback into
  one that said why first. `--no-cache` on a full slide still needs that memory; it always
  did.
- **Freed memory now returns to the OS between elements.** glibc keeps freed blocks in
  per-thread arenas, so RSS ratcheted to the high-water mark and stayed there — which is
  what made the write loop look like it leaked when nothing was actually retained.
  `write_sdata` calls `gc.collect()` + `malloc_trim(0)` after each element.

#### Investigated, no change
- **napari's global opportunistic dask cache is not implicated.** `Layer.__init__` calls
  `configure_dask(data, cache=True)` for any dask-backed layer, which sizes
  `dask.cache.Cache` at 25% of RAM — **33.77 GB** here — and `dask.callbacks` registration
  is process-global, so a slicing window overlapping a long write does register it. But
  cachey scores entries by cost *per byte*, and large chunks score near zero: measured
  retention during a registered window was **0.2 MB out of 32 MB computed**. Bounding it
  would have been a plausible-sounding fix for a problem it does not cause.
- `FLUSH_THRESHOLD` in `preprocess.py` stays at 10M (see 2026-08-16).
- **Crop & Export still uses `sdata.write()` under a blanket 2-worker cap**, rather than
  being moved onto `write_sdata`. It would gain per-element progress and let its image write
  run at full width, but that path was verified at real scale only days ago and is bounded
  as it stands; changing it is a separate, separately-verified piece of work.

### 2026-08-16

#### Performance
- **The transcript feather cache rebuilt ~13× more data than it needed to** — measured, not
  estimated: of the ~69 minutes a verified real-scale crop took, ~53 minutes (77%) was
  `preprocess.py`. Feather (Arrow IPC) has no append, so `_flush_buffers` **read back and
  rewrote every gene's entire accumulated file on every flush cycle**; with ~13 cycles that
  is `sum(1..13)/13 ≈ 7×` the final size in writes plus ~6× in reads, and the same
  multiplier on the `pd.concat` cost — ~66,000 read+write operations where 5,101 writes
  would do. (The measurement was against local disk, so this was genuine rewrite cost, not
  NAS latency.) Replaced with **one open `pa.ipc.new_file` writer per gene**, so each row is
  written exactly once (bounded by the process's real fd limit — see the two fd bugs
  below, which is where the first attempt at this went wrong). Two
  supporting changes: all batches are written against **one fixed schema** (decoding
  dictionary columns to their value type — `feature_name` comes back from `to_pandas()` as a
  dictionary whose categories vary batch to batch, which the writer rejects on the first
  mismatch, the same class of bug hit in `crop_export.py`), and `iter_batches(columns=…)` now
  reads only the 6 columns used instead of decoding all 13 and discarding most.

  **Measured: full 218.5M-row dataset in 6.2 min / 1.40 GB peak RSS**, versus a 52.6-min
  baseline for 127.9M kept rows — ~13× faster normalised per row. Correctness was pinned
  first: old vs new over a real 3M-row slice forced through ~6 flush cycles produce
  **identical output** — same 5,101 genes, same row counts, same values, and the same dtypes
  (`float32` preserved, not upcast). Benefits the `xenium-preprocess` CLI and `app.py`'s
  one-time preprocessing worker too, not just cropping, since all three share this function.

  **Two fd-exhaustion bugs were found and fixed while validating this**, both of which had
  to be reproduced under the *viewer's* real limits rather than a shell's:
  - A pool of one writer per gene needs one fd per gene (5,101 for this panel), but a
    desktop-launched process inherits systemd's default soft `RLIMIT_NOFILE` of **1024**
    (hard ~1e6) — an interactive shell reports ~1e6 because conda's init already raised it,
    which is exactly how the unbounded design got waved through. `preprocess()` now raises
    its own soft limit best-effort for the duration (restoring it after; lowering a soft
    limit never invalidates open fds), **and** caps the pool at `limit - 256`. Genes beyond
    the cap fall back to one shard per flush, merged at the end into the same single
    `<gene>.feather` — so the reader is unaffected and the failure is structurally
    impossible rather than merely unlikely.
  - `RecordBatchFileWriter.close()` **does not release the underlying file descriptor**
    when pyarrow opened the sink from a path; the fd survives until garbage collection
    (measured: 200 writers, all closed, still holding 200 fds until `gc.collect()`). Sinks
    are now created as explicit `pa.OSFile` objects and closed alongside their writers.

  All three fd regimes are verified to produce output **identical to the pre-change
  implementation** (same genes, rows, values, dtypes): the viewer's real soft-1024 case,
  a worst case where the limit cannot be raised at all (4,322 genes sharded), and a forced
  8-writer budget (5,082 genes sharded) — the fallback is exercised deliberately rather
  than left as untested branch.

  Also: gene files are now written as `<gene>.feather.partial` and renamed only once every
  writer is closed. An IPC file has no footer until closed, and
  `TranscriptLoader._cached_genes` globs `*.feather` **without** consulting the `.complete`
  sentinel — so an interrupted rebuild previously left files the viewer would happily pick
  up. `FLUSH_THRESHOLD` moved to a module-level constant beside `CHUNK_SIZE`/`MIN_QV` so the
  flush path is tunable and testable (a single-flush test passes trivially and proves
  nothing). Its value is deliberately left at 10M rows: with the rewrite gone it only
  controls batch granularity and peak buffer memory, and raising it would trade back the
  memory headroom the crash fix above exists to protect.

#### Fixed
- **Crop Dataset could exhaust memory and crash the machine on a large crop** — the real
  cause, found after an initial (real but insufficient) transcript-filter fix still
  crashed on a retry: `crop_and_export()` (`utils/crop_export.py`) pulled the cropped
  morphology image and both label rasters fully into dense numpy arrays
  (`np.asarray(...compute())`) before ever touching disk. For a crop spanning roughly
  half a slide (confirmed from the crash's own leftover partial output: 27043×51144px,
  full width / half height — this dataset has two tissue sections sharing one run,
  cropped out separately), that's ~11GB for the image alone plus ~5.5GB each for the two
  label rasters, before `np.isin`/`np.unique` working-buffer overhead — on top of
  whatever the running viewer already held. Nowhere else in the codebase does this;
  the initial cache build already writes the full, larger, uncropped image the same
  dask-backed way. Fixed by keeping the image/label crops as dask arrays through
  `Image2DModel.parse`/`Labels2DModel.parse` (both accept a dask array directly) instead
  of materializing them first, so `new_sdata.write()` computes and writes them chunk by
  chunk. The transcript path needed the same treatment for the same crop shape (a
  full-width crop barely benefits from `filters=` bbox pruning, since the on-disk
  partitions are banded by x, not y) plus two more real bugs surfaced only by testing at
  the actual crash size: `PointsModel.parse`'s own internal index-monotonicity check
  forces an eager `.compute()` that (like the image/label case) can pull many large
  partitions into memory at once under dask's default ~40-way scheduler — capped to a
  small worker count for this path specifically; and `feature_name` (5000+ possible gene
  names) needs `.cat.as_known()` before any write, or partitions that happen to see
  different numbers of distinct genes get inconsistent categorical index widths and the
  parquet write fails outright with a schema mismatch (the fix mirrors spatialdata's own
  `write_points()`, which hits the identical issue). The exported `transcripts.parquet`
  also has to stay a single *file* (matching real Xenium output, which the rest of the
  app assumes) rather than the directory dask's own `.to_parquet()` would produce, so
  it's streamed through one shared `pyarrow.parquet.ParquetWriter` across partitions
  instead. **Verified against the actual dataset and crop size that crashed**: peak RSS
  9.76GB (previously crashed the machine), 287,548 cells / 145,605,025 transcripts
  written correctly, full pytest suite green.
- **Crop Dataset's transcript-filter step loaded the entire source transcripts table on
  every crop, regardless of crop size** — `ctx.sdata.points["transcripts"]` (spatialdata's
  dask dataframe over `points.parquet`, opened with no predicate) fed straight into
  `.map_partitions(...).compute()`: every partition of a potentially
  hundreds-of-millions-of-row table was loaded and mostly discarded, for every crop. A
  boolean bbox mask applied to the already-open dask dataframe does *not* get pushed down
  into the parquet reader (confirmed empirically — identical wall-clock with or without
  the mask); re-reading the same on-disk `points.parquet` directory with an explicit
  `filters=` kwarg does get real row-group pruning where the crop is narrow in x (the
  axis partitions are banded by) — e.g. ~40s → ~1-2s for a small crop — though as the fix
  above found, this alone doesn't help (and isn't relied on for correctness/memory
  safety) once the crop spans most of the source table's width.

### 2026-08-05

#### Changed — templates now use the library API instead of hand-rolled equivalents

The `.tmpl` files are unusual code: they are `exec`'d by the GUI *and* recorded verbatim
into the reproducible notebook, so a reader learns the analysis from them. Hand-rolled
numpy/shapely that re-implements a library call is therefore worse here than in ordinary
code — longer, likelier to be subtly wrong, and it teaches a non-idiom. A review of all
14 templates found the four squidpy ones clean and no hand-rolled neighbour-graph
construction anywhere; it found four places worth replacing. `CLAUDE.md` records the
principle as rule **(e)**, beside the four that already exist because each was once broken.

- **ROI membership is `spatialdata.polygon_query`.** `roi.deg` and `roi.expression` each
  carried the same five-line `shapely.contains_xy` loop, with three defects: a
  `[:, ::-1]` flip undone by the `[:, 1], [:, 0]` read on the very next line (an
  `n_obs × 2` copy for nothing, and a comment asserting a convention the code reversed);
  `buffer(0)` to repair a self-intersecting ROI; and a µm→pixel conversion applied to
  every centroid rather than declared once. Cell centroids are now a `PointsModel` whose
  scale is a **declared transformation**, and `roi.polygons` binds valid shapely
  geometries in (x, y) pixel space — the same frame `sdata.shapes['rois']` already uses,
  so the napari (y, x) flip happens once at the source.

  Two user-visible consequences:
  - **A self-intersecting ROI keeps every lobe.** `buffer(0)` silently *deleted* one —
    on a bow-tie it returned half the area. `shapely.make_valid` repairs the shape
    instead, so a polygon drawn across its own edge no longer loses cells without warning.
  - **Boundary semantics.** spatialdata tests `intersects` where `contains_xy` was
    strict, so a cell centroid landing exactly on an ROI edge is now inside. For float
    centroids this is measure-zero: on a real 63k-cell dataset with two overlapping
    L-shaped ROIs, membership was **identical** to the old loop, cell for cell.

  Overlapping ROIs keep their last-wins labelling, which the template now states rather
  than leaving implicit.

- **Rank genes are stored per clustering.** `sc.tl.rank_genes_groups(key_added=…)`
  replaces the notebook-local `rank_results` dict that worked around scanpy overwriting
  `uns['rank_genes_groups']` in place. The dict kept every ranking but put it where no
  `sc.pl.rank_genes_groups*` function could reach it — they all take `key=`.

  Results now live in **`uns['rank_genes_<clustering>']`** in the persisted table, so
  ranking a second clustering adds a result instead of replacing the first, and the
  recorded dotplot and panel-plot cells pass the matching `key=`. `rank_genes_groupby`
  still names the most recent one. **Existing caches keep working**: every read resolves
  through `gene_analysis.resolve_rank_key`, which falls back to the unkeyed
  `rank_genes_groups` that is all a pre-existing cache holds. Cache-rebuild warnings and
  Tools → Dataset deletion recognise the keyed names too — both matched a literal string
  before, so a keyed ranking would otherwise have stopped counting as user data.

- **Two smaller substitutions.** `roi.expression` reads expression through
  `sc.get.obs_df` (one call, no `.todense()` branch, no per-region `DataFrame` +
  `concat`); the inferCNV template and `cnv_analysis.run_cnv_pipeline` take the CNV row
  mean with `np.abs(X).mean(axis=1)`, which works natively on CSR instead of densifying
  the whole `n_cells × n_bins` matrix. `genes.correlation`'s fraction-of-total branch
  takes totals from `sc.pp.calculate_qc_metrics(…, inplace=False)` — **`inplace=False`
  deliberately**: `inplace=True` writes `obs['total_counts']`, one of Xenium's own
  columns, and would overwrite the raw value on the live object.

`spatialdata` is now bound as **`sd`** in the template namespace and imported by the
recorded notebook preamble.

#### Documentation
- **Comprehensive wiki review.** `docs/` doubles as the GitHub Wiki and the mkdocs
  source; 20 of 23 `Tab-*.md` pages had not been touched since 2026-06-26, and
  `mkdocs.yml` had not moved either. Audited every page against the code it documents.

  **Three shipped tabs had no page and appeared in no navigation** — Tools → Dataset,
  Cache and Templates, all added in the last fortnight. `docs/Tab-Dataset.md`,
  `docs/Tab-Cache.md` and `docs/Tab-Templates.md` now exist, and the Tools section of
  `_Sidebar.md`, `Home.md`, `Interface-Overview.md` and the `mkdocs.yml` nav lists all
  seven tabs rather than four. The nav had also been missing `Tab-CNV.md` and
  `Tab-Crop-Dataset.md` since they were written.

  **`Tab-CNV.md` predated the entire CopyKAT backend.** It described a single-backend
  tab; it now covers the backend selector (defaulting to *both*), CopyKAT max cells,
  call extrapolation, the cell-type restriction grid, the per-backend/per-resolution
  heatmap selectors, the human-only panel guard, and the detached second-environment
  run. Its output paths (`plots/cnv_heatmap_<backend>_<key>.*`) and clustering keys
  (`cnv_leiden_res<res>`, not `cnv_leiden_<res>`) were both wrong.

  **Corrections that mattered most**, each verified against the code rather than
  assumed: `Tab-Transcripts.md` claimed **Min QV** filtered at display time and that
  **Show transcripts** toggled the layer — neither is connected to anything at runtime;
  `Tab-Cell-Coloring.md` claimed colouring state was restored across launches, and it is
  not persisted at all; `Tab-UMAP.md` claimed a plot without a clustering saved
  uncoloured, when it writes nothing; `Tab-Rank-Genes.md` claimed LLM annotation needs
  an API key, when it shells out to a local CLI; `Installation.md` named a
  `sdata_cached_corrupt_*` file the code no longer produces and described staleness as
  an mtime check rather than a content hash. Six pages omitted that their figure is
  written to `plots/` (SVG by default) with no dialog. `--no-user-templates` and
  `XENIUM_VIEWER_TEMPLATE_PATH` were undocumented despite being the designated
  "results in doubt" escape hatches.

  New `docs/Tutorial-Recovering-a-Cache.md`, because recovery spans the startup dialog,
  the Cache tab and the Dataset tab, so no reference page owns it. New `## Menu Bar`
  section in `Interface-Overview.md`, which had never documented the three menus.

#### Fixed
- **Three widgets were built, filled, and never laid out.** `ga_results_text` (the Rank
  Genes top-50 preview), `lr_results_text` (the Lig-Rec interaction counts and top-20
  table) and `reg_residuals_qt` (the H&E per-landmark residuals) were each constructed,
  connected and populated on every run, then omitted from their tab's `make_tab(...)`
  call — so the work was done and discarded, and the results reached the user only as a
  one-line status message. The H&E one is the sharpest case: `Tab-HE-Registration.md`
  has always documented a "Residuals (read-only)" control, and per-landmark residuals
  are the only way to see *which* landmark is dragging a fit, where the status bar shows
  the mean alone. Found by diffing each page's control table against the `make_tab`
  argument list, which is the definitive record of what renders.

- **The screenshot script mislabelled a tab.** `capture_screenshots.py` mapped Tools
  index 2 to `tab-notebook.png`, but index 2 has been Crop Dataset since that tab was
  inserted — so every capture run overwrote the Notebook screenshot with a picture of
  Crop Dataset. Indices are positional into `app.py`'s `addTab` order; the missing CNV
  and Tools 3–6 entries are added.

- **An internal planning note was published to the public wiki.**
  `push_to_wiki.sh` excluded repo-only files by a hand-maintained list, which had to be
  extended after each leak and had already missed `user-configurable-templates-todo.md`.
  Replaced with a convention — a wiki page is `Title-Case.md` or `_Sidebar.md`, and
  internal notes are lower-case — so the next note is repo-only without anyone
  remembering. `tests/test_docs_links.py` parses the pattern out of the script, keeping
  it the single source of truth.

#### Added
- **`tests/test_docs_links.py`** — the first automated check on `docs/`. Asserts that
  internal wiki links resolve, that referenced screenshots exist, that the set of
  published pages, the mkdocs nav and the sidebar agree, and that every `addTab` label
  in `app.py` has a reference page (read with `ast`, so no Qt import). That third
  assertion alone would have caught every navigation gap found in this review. Pure
  stdlib, no dataset, runs in the existing CI job in milliseconds.

### 2026-07-31

#### Fixed
- **The verifier's own notebook never said a session was customised.**
  `scripts/verify_notebook.py` builds its notebook from
  `prov_graph.graph_to_cells` rather than the `notebook_export` one, because
  that is the version carrying `node_id` — which is what turns nbclient's
  "cell 4 failed" into a named step. The customisation banner is added by the
  other one, so it was silently dropped: the `--work-dir` notebook was the
  single artifact of a customised session that did not say so, while the
  viewer's own export and the `--out` report both did, and it is the artifact
  most likely to be kept or forwarded. Found by running the round-trip on a
  real dataset. The banner is now prepended explicitly, with a matching `None`
  in `node_ids` — `node_of()` indexes that list by absolute cell index, so
  inserting a cell without shifting it would have misattributed every timing in
  the report. Both halves are guarded by tests.

- **Marker-gene plots ignored user templates entirely.** `tab_marker_genes` built its
  `Step` with `template=builtin_assemble(...)` — the shipped-files-only path, which by
  design cannot see an override directory. So `genes.marker_plot` could be edited,
  validated and saved in Tools → Templates with no effect whatsoever, and its provenance
  nodes carried no `template_id`, leaving them invisible to the notebook's customisation
  banner and to `verify_notebook`'s `stock_templates`. Every other migrated call site
  splatted `step_template`; nothing checked that they all did. Now one does
  (`test_every_step_resolves_user_overrides`), by parsing every `Step(...)` in every tab.
  `builtin_assemble` remains correct where it is still used — in the `_*_template`
  helpers the pinning tests read, which must not see a developer's own overrides.

#### Changed
- **The Templates preview now shows what the button would actually run — for every
  template, in shape as well as in value.** Two gaps, one visible and one not. Only the
  Clustering tab supplied live parameters, so twelve of fourteen panes read
  `groupby='sample', n_genes=1`. And the preview always rendered the template's *first
  declared assembly*, so even Leiden's "real" preview kept the same code shape however
  the checkboxes were set — untick "use HVGs" and the numbers moved while the statements
  did not.

  Both follow from the same contract. A tab now registers a `Preview(blocks, params)`
  in `ctx.state["template_preview"]`, and **its own callback runs from that same
  call** — so the blocks selected and the parameters passed are one expression, not two
  that agree by discipline. Block selection has to be in there: the branch structure *is*
  what the widgets mean, which is exactly why it stays in Python and cannot be re-derived
  by the pane. Twelve templates have a provider; the two that do not are declared, with
  their reasons — `normalize` takes no params, and `spatial_neighbors` takes its `k` from
  whichever tab called it, so it gets a realistic literal from a new `# sample-params:`
  header field instead of a provider that would have to pick a slider arbitrarily.

  `note` covers the values that cannot come from a widget because they do not exist yet:
  a save-dialog path renders as the filename the dialog would propose, and the header
  says `(path chosen on save)` rather than showing a placeholder as settled.

  Guarded four ways, each of which caught a real defect when tried against it: every
  template has a provider or a declared exemption; every provider is *called* by its own
  tab rather than shadowed by a second inline dict; every provider actually answers when
  invoked against a live tab (a raising provider is otherwise indistinguishable from a
  half-built one, since `_preview` catches it and shows sample values); and every
  provider selects a block sequence the template declares.

### 2026-07-30

#### Added
- **A customised template survives an upgrade, and says when it needs a second look.**
  This is the `dpkg` conffile problem, and most of it turned out to be already solved by the
  shape of the storage rather than by logic: an override records *only* the blocks the user
  changed, so every other block resolves against whatever the current release ships. The
  "unmodified file, replace silently" case `dpkg` has to detect **cannot arise here** — there
  is nothing to detect, and no prompt to dismiss.

  What is left is the genuinely hard case: a block the user *did* change whose shipped version
  has since changed too. `overrides.json` records the hash of the **shipped** text each block
  was forked from — of the shipped text, not the user's, because the question a later release
  must answer is "has the thing they diverged from moved?". The edit still applies; silently
  reverting someone's method would be far worse. But it is badged `⚠ review`, with a two-way
  diff of theirs against the new default and a **Take new default** button that updates only
  the blocks that moved, leaving their other customisations intact. Saving again is itself the
  act of reviewing, and clears the flag. No three-way auto-merge: its conflict markers land in
  Python source, where a stray `<<<<<<<` is a syntax error rather than a visible annotation.

  A conflict that no longer *validates* — a release drops a param the forked block still
  references — is deactivated outright rather than flagged, since there is nothing to review.
  A corrupt or missing manifest costs the warning, not the override: losing bookkeeping must
  not lose the user's work.

- **Customisation is visible downstream, not just in the tab.** The exported notebook gains
  one markdown banner at the top when any step used a non-shipped template, naming the steps,
  their template ids and hashes — because a customised template renders to code that looks
  entirely ordinary, so the source alone cannot tell a reader this is not the stock pipeline.
  `hand-edited` is called out separately as the one origin whose code may not describe what
  produced the result. `scripts/verify_notebook.py` gains a `templates` section with per-step
  origin and hash plus a top-level `stock_templates` bool, and prints it: replay agreement
  proves reproducibility, not that the pipeline was the standard one.

- **Analysis templates can now be customised, per user.** An edited template is used both
  by the GUI and by the recorded notebook — which needs no special machinery, because the
  Step system already renders one string and hands it to both.

  Overrides live in `~/.config/xenium-viewer/templates/*.tmpl` (`platformdirs`, so
  `XDG_CONFIG_HOME` is respected) and are **resolved per block**. That is the load-bearing
  choice, not a tidiness one: most of a fork is text the user never touched, so those blocks
  keep tracking the shipped template and a later fix to them still reaches everyone who
  customised a *different* part. Whole-file override would freeze the entire template at the
  version it was forked from, which is how someone quietly stops receiving a correctness fix.
  Saving writes only the blocks that actually differ, so this is the default rather than
  something the user has to think about.

  **Nothing is trusted without validation** (`step_templates/validate.py`, which promotes
  `check_step`/`free_names` from test-only helpers to a production gate). Per legal assembly:
  it must render, parse, read only names the executor guarantees or the template declares,
  bind every output it claims, and leave frozen blocks alone. The check weighted most heavily
  is that **a required param the template no longer mentions is a hard stop** — that template
  would run, report success, and silently ignore a setting the user chose, which is far worse
  than a crash. Explicitly *not* a security boundary, and the module says so: a user can
  already run arbitrary Python in the Notebook tab.

  **A refused override is loud and never fatal.** It is skipped, the shipped template runs, a
  napari warning fires once per template per session, the Cache tab's health line mentions it,
  and Tools → Templates badges it `✕ not used`. A broken file must not make the viewer
  unlaunchable — the user would have no way in to fix it. Two off switches, both load-bearing
  rather than conveniences: `--no-user-templates` (the first thing to try when a result is in
  doubt) and `XENIUM_VIEWER_TEMPLATE_PATH` (which `tests/conftest.py` empties, so a
  developer's own customisations can never change what the suite asserts).

  Every step now carries `template_id` / `template_origin` / `template_hash` into the
  provenance graph, so a reader can tell a stock run from a customised one — which rendered
  source alone cannot show.

  Tools → Templates gained the editing half: **Default (read-only) beside Yours (editable)**,
  with Validate, Save & Activate, and Revert. **Save never refuses.** Refusing to write would
  send the user to an external editor and out of the feedback loop; what is gated is
  *activation*, and that needs no special mechanism — an invalid file on disk is simply
  rejected by the resolver, which falls back and says so.

  Three bugs found while building this, each by a test rather than in review: the header
  parser treated an unindented prose line as a continuation of the field above it (so a saved
  override's own explanatory comment was appended to `schema-version` and failed to parse);
  the frozen-block annotation was written onto the block-marker line, where anything after the
  name *is* the name, so editing the CNV Arrow shim silently created a new block instead of
  the protected one; and saving resolved its destination independently of reading, so the tab
  tests wrote into the real `~/.config`. Writes now derive their destination from the same
  search path reads use, so a write cannot land where the reader does not look.

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

#### Changed
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

#### Fixed
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

#### Added
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

### 2026-07-29

#### Fixed (found by replaying a real session)
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

#### Added
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

#### Fixed
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

#### Added
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

#### Changed
- **Sidecar analysis outputs moved out of the zarr store root** into
  `<data_path>/viewer_cache/`. Files in the store root make zarr's hierarchy walk emit a
  `ZarrUserWarning` each — the most likely source of the reported "several warnings",
  since `app.py` called `consolidate_metadata` without the filter spatialdata itself uses.
  It also meant a cache rebuild deleted them, including hours of CopyKAT compute. Readers
  fall back to the old location, so existing datasets keep working and nothing is migrated
  eagerly. (`utils/adata_persistence.py`, `tabs/tab_cnv.py`, `loader.py`)

#### Tests
- First coverage the zarr/persistence paths have ever had: `test_zarr_safe.py` (26,
  including interrupted-write simulation with both a recoverable `OSError` and a
  `KeyboardInterrupt` that bypasses cleanup the way a `kill -9` does),
  `test_persistence_safety.py` (9), `test_session_persistence.py` (14),
  `test_cache_repair.py` (20), `test_loader_policy.py` (16), `test_sidecar_location.py`
  (20), `test_reporting.py` (21), `test_tab_cache.py` (24). Plus source guards that fail if `delete_element_from_disk` is called outside
  `zarr_safe.py`, if `loader.py` `rmtree`s the live cache, or if a sidecar is written into
  the store root, or if a cache write path prints a warning instead of logging it, or if recovery opens a backup as a whole.
  289 tests pass.

### 2026-07-28

#### Added
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

#### Fixed
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

#### Removed
- **`spatial_analysis.run_ligrec` and `spatial_analysis.run_co_occurrence`.** One squidpy
  call each; the templates in the tabs now *are* that call. (`utils/spatial_analysis.py`)
- **`utils/leiden_worker.py`.** The spawned-subprocess Leiden existed for GUI
  responsiveness, but it was also the second expression of the algorithm that drifted from
  the recorded one. scanpy's `igraph` flavor is, per its own warning, "orders of magnitude
  faster" than `leidenalg`, which removes the motivation. (`utils/leiden_worker.py`,
  `cnv_copykat_worker.py` docstring reference)

#### Changed
- **`.gitignore`**: added `manuscript/` (preprint drafts and planning notes, kept out of
  the public repo) and `data/` (untracked local datasets).

### 2026-07-19

#### Fixed
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

#### Docs
- **Tracked TODO to migrate napari off the deprecated PyQt5 backend.** Napari warns
  that PyQt5 support is deprecated and will be removed in fall 2026. No code change
  yet (the viewer runs fine on PyQt5 until then); the codebase already routes all Qt
  access through `qtpy`, so the migration is small — only the backend pins plus a few
  Qt5-isms (8 unscoped enums, 7 `.exec_()` calls). Captured the migration checklist in
  `docs/pyqt6-migration.md` for a future session. (`docs/pyqt6-migration.md`)

#### Added
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

#### Fixed
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

#### Added
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

### 2026-07-17

#### Added
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

### 2026-07-16

#### Added
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

#### Fixed
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

### 2026-07-15 — b

#### Changed
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

### 2026-07-15

#### Added
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

### 2026-07-14 — e

#### Added
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

### 2026-07-14 — d

#### Fixed
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

### 2026-07-14 — c

#### Fixed
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

### 2026-07-14 — b

#### Fixed
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

### 2026-07-14 — a

#### Fixed
- **Crop Dataset: orphan nuclei in `nucleus_labels`** — nucleus IDs are their own
  independent numbering, unrelated to `cell_id`/`cell_labels`, so an ID-based overlap
  check (added, then removed, in the previous fix) couldn't reliably mask them and
  the fallback left `nucleus_labels` completely unmasked. Now masked correctly via
  spatial overlap with the already-masked `cell_labels` crop: a nucleus is kept only
  if it occupies at least one pixel within a kept cell's footprint.

### 2026-07-13

#### Added
- **Crop Dataset (Tools tab)** — draw one or more polygons in a new "Crop Regions"
  napari layer and export each as its own standalone, independently-openable
  xenium-viewer data directory (cropped morphology image, cell/nucleus labels,
  transcripts, and AnnData table). Images/labels are cropped to each polygon's pixel
  bounding box; cells and transcripts are filtered to the exact drawn polygon via
  true point-in-polygon tests. Output/name for each region is chosen via sequential
  folder-picker + name-prompt dialogs, run in a background thread with a progress
  dialog. New module `src/xenium_viewer/utils/crop_export.py`; new tab
  `src/xenium_viewer/tabs/tab_crop_dataset.py`.

### 2026-07-01

#### Fixed
- **Permission-denied dialog gaps** — two write paths were missing coverage from the
  read-only zarr dialog: (1) `delete_element_from_disk` in the Patches tab only printed
  to console on failure; (2) `record_code` in `_helpers.py` had no exception handling at
  all when writing `code.py`. Both now call `_maybe_show_permission_dialog` on
  `PermissionError`/`OSError`.

### 2026-06-30

#### Added
- **tab10 palette in Patches tab** — `tab10` (matplotlib's 10-colour categorical
  palette) is now available in the Patches tab palette dropdown alongside tab20,
  glasbey_dark, Set1, Set3, and ARMS.

### 2026-06-26 — d

#### Fixed
- **Leiden clustering UI freeze** — `sc.tl.leiden` held the Python GIL during graph
  construction and partition, causing the progress bar and status-bar spinner to freeze.
  The Leiden step now runs in a subprocess via `ProcessPoolExecutor` (spawn context),
  giving it its own GIL so the main process's Qt event loop stays responsive throughout.
  New module: `src/xenium_viewer/utils/leiden_worker.py`.

### 2026-06-26 — c

#### Added
- **Progress bar for long-running analyses** — an indeterminate `QProgressBar` now
  appears inside the control panel directly below the run button while an analysis is
  in progress (Leiden clustering, rank genes, L-R, neighbourhood enrichment,
  co-occurrence, ROI DEG, annotation nhood enrichment, Novae domains). The bar
  disappears automatically when the analysis finishes or errors out. The existing
  napari status-bar spinner/tqdm text is retained unchanged.

### 2026-06-26 — b

#### Added
- **Wiki screenshots** — all 52 documentation screenshot placeholders are now filled with
  actual PNGs captured from a running viewer instance. A new script
  `scripts/capture_screenshots.py` automates future recapture by programmatically
  navigating each tab and grabbing the control-panel and full-window views.
  One placeholder (`tutorial-clustering-step5.png`, the matplotlib dotplot window)
  remains a comment for manual capture.

### 2026-06-26

#### Fixed
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

### 2026-06-25

#### Fixed
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

#### Added
- **ARMS Overlay: save/load landmarks** — "Save Landmarks..." and "Load Landmarks..."
  buttons added to the ARMS Overlay tab (below "Compute Registration"), mirroring the
  H&E Registration tab. Landmarks and the computed affine are saved to a portable JSON
  file via the existing `save_landmarks` / `load_landmarks` API. The save button
  enables when ≥1 landmark pair is placed and disables on "Clear All".

#### Changed
- **H&E Registration: removed Save Affine button** — the "Save Affine..." button has
  been removed. There was no corresponding "Load Affine..." button, making it a dead end.
  The affine is already auto-persisted to the sdata zarr cache on every registration.

### 2026-05-05

#### Changed
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

### 2026-04-18

#### Changed
- **External Images tab: single composite layer** — multichannel images (e.g. 25-channel PhenoCycler)
  now display as a single RGB composite napari layer instead of one layer per channel. Per-channel
  visibility checkboxes, color buttons, and contrast range sliders are in the tab widget.
  Composite is built lazily via dask so large images render efficiently.
- **External Images tab: landmark registration** — external images can now be independently aligned
  to the Xenium image via landmark-based registration (same `compute_landmark_affine` as H&E tab).
  Includes flip V/H checkboxes. The "Apply transform from" dropdown remains as an alternative.
  Landmarks are persisted to sdata.shapes; affine to sdata transformations.

### 2026-04-17

#### Fixed
- **Overlay affine persistence** — affine transforms for patch overlays and external images are now
  saved to the SpatialData object via `set_transformation()` / `write_transformations()` (same
  pattern as H&E registration). On restore, the affine is read back from sdata directly, so
  overlays are correctly positioned even before the source layer (e.g. H&E) finishes loading.
  A deferred-linking listener also re-establishes live affine mirroring once the source layer
  appears.
- **QCheckBox crash on close** — `_snapshot_layers` now reads `entry["hidden_cluster_ids"]` (a
  plain set maintained by the tab) instead of iterating Qt checkbox widgets that may already be
  destroyed during shutdown.

#### Changed
- **ARMS palette for subclone predictions** — subclone prediction overlays now default to the
  ARMS palette (RColorBrewer Set1+Set2+Dark2, 24 colours) to match the R-based ARMS visualisation
  pipeline. Cluster-to-colour mapping is 0-normalised so 1-based genomic cluster IDs (1,2,3) map
  to Set1 red, blue, green — matching the R package exactly.

### 2026-04-16

#### Added
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

### 2026-04-13

#### Added
- **Cluster labels persisted to SpatialData** — user-assigned cluster names (from the label editor,
  reference atlas annotation, and LLM annotation) are now saved as `adata.obs["cluster_labels_<key>"]`
  columns (e.g. `cluster_labels_leiden_r1.0`) inside `sdata.tables["table"]` immediately on
  assignment. Labels are loaded back on the next launch and merged with the session-attrs fallback
  (sdata wins on conflict). The obs columns are readable in any standalone Python session.

### 2026-04-12

#### Changed
- **ROI DEG and ARMS tile DEG persistence migrated to SpatialData** — results are now saved as
  sidecar parquets inside the zarr cache (`roi_deg_cache.parquet`, `arms_tile_deg_cache.parquet`)
  immediately on computation rather than at session close. Restores automatically on relaunch.
  One-time migration copies old `viewer_session/*.parquet` files to the new location.
- **ARMS tile DEG code recording fixed** — generated code snippet in `code.py` is now fully
  executable: loads tile polygons from `sdata` via `load_arms_tiles_from_sdata`, reconstructs the
  ARMS registration affine (fine × flip), applies it to the tiles, filters by selected tile
  clusters, and optionally applies a Xenium cell cluster mask before calling `compute_arms_tile_deg`.

### 2026-04-10

#### Added
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

### 2026-03-26 — 2

#### Added
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

### 2026-04-01

#### Added
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

#### Changed
- **ARMS tiles: ColorBrewer Set1+Set2 color palette** — replaced the hardcoded 8-color custom
  palette with a concatenation of ColorBrewer Set1 (9 colors) + Set2 (8 colors) = 17 distinct
  colors. Cluster filter checkboxes and the legend now show `C{id}` labels without color names.

### 2026-03-27 — 3

#### Fixed
- **Annot. Nhood and Annot. Distance clustering dropdowns not refreshed on segmentation swap** —
  both tabs created their `clustering_widget` ComboBoxes without registering them on `ctx`, so
  `refresh_clustering_choices` skipped them. Fix: register as `ctx.annot_nhood_clustering_widget`
  and `ctx.annot_dist_clustering_widget`, add corresponding fields to `ViewerContext`, and include
  them in the refresh loop in `_helpers.py`.

### 2026-03-27 — 2

#### Fixed
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

### 2026-03-27

#### Added
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

### 2026-03-26

#### Added
- **View menu: Show Minimap toggle** — checkable "Show Minimap" action appended to napari's native
  View menu. Enabled and checked when the minimap overlay is available; disabled (grayed out) when
  there is no morphology data. State is reset on dataset reload.

### 2026-03-25

#### Changed
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

#### Fixed
- **Index alignment for adata persistence** — `adata.obs.index` uses integer strings (`'0'`, `'1'`, ...)
  while clustering series and UMAP are indexed by cell barcode (`'aaaagflk-1'`, ...). Save/load now
  maps between the two via the `cell_id` column.
- **Show/Export buttons disabled at startup for persisted analyses** — nhood enrichment,
  co-occurrence, and L-R results were loaded from `adata.uns` *after* tab session restore ran,
  so the buttons stayed disabled despite results existing. Fixed by loading `adata.uns` analysis
  results before calling `restore_fn()`.

### 2026-03-23

#### Added
- **Level slider in Novae tab** — exposes the `level` parameter (hierarchical tree level,
  default 7, range 1–15) for `assign_domains()`, giving finer control over domain granularity.
- **Console variable injection** — key variables (`adata`, `sdata`, `viewer`, `ctx`,
  `clusterings`, `color_manager`, `gene_names`, `data_path`) are now pushed into napari's
  built-in IPython console on dataset load, with a help message listing what's available.
  Variables are refreshed on dataset reload. Enables a GUI-to-code handoff workflow alongside
  the existing code recording feature.

#### Fixed
- **Console button hidden by layer widgets** — capped layer controls and layer list dock
  widgets to 200px max height so the console toggle button stays visible.
- **Console not resizable when opened** — the Xenium Controls dock (QTabWidget with 13 tabs)
  had a ~970px minimum height from stacked widget content, leaving no room for the console.
  Fixed by wrapping tab content in `QScrollArea` inside `make_tab()`, dropping the dock's
  minimum height to ~117px. Also defers the `resizeDocks` call via `QTimer.singleShot(0)` to
  ensure it runs after Qt finishes the visibility layout pass.

### 2026-03-20

#### Added
- **Spatial Domains tab (Novae)** — new tab that runs [Novae](https://mics-lab.github.io/novae/)
  zero-shot spatial domain inference. Select species (human/mouse), optionally specify N domains
  (0 = auto-detect), and click "Run Novae Domains". On completion the cell labels layer is
  automatically recolored by inferred domains, the `novae_domains` key is added to the clustering
  dropdowns, and results are persisted to the session cache so they are restored on relaunch.
  Requires `pip install novae`. The full pipeline is recorded to `code.py`.

### 2026-03-19

#### Added
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

#### Changed
- **Auto-save plots on Show** — clicking "Show" in the Neighborhood Enrichment,
  Co-occurrence, Ligand-Receptor, Gene Correlation, and Gene Analysis (dotplot) tabs now
  automatically saves the figure to `<data_dir>/plots/<stem>.<format>` and reports the
  path in the status bar. The separate "Save Plot" button has been removed from all tabs.
- **Default plot format is now SVG** — plots are saved as vector graphics out of the box;
  PNG remains available via Preferences → Plot Format.

#### Added
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

#### Fixed
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

#### Changed
- **"Edit Cluster Labels..." moved to Cell Colouring tab** — the button is now in the Cell
  Colouring tab (below "Apply Cell Coloring"), where it is more naturally discovered alongside
  the clustering selector and cluster filter controls. Removed from the Clustering tab.

#### Added
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

#### Fixed
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

#### Added
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

#### Changed
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

#### Fixed
- **"Filter by selected clusters" in Transcript Density** — the checkbox was
  silently ignored when `label_to_cluster` was not yet populated (e.g. before
  clicking "Apply Cell Coloring" or after switching to gene coloring). The fix
  pre-computes filter data on the main thread in `on_compute_density()`: if a
  clustering is selected it builds `label_to_cluster` on-the-fly when needed;
  if no clustering is available it shows "No clustering applied — filter
  skipped" in the status label and returns early instead of silently falling
  through.

#### Added
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

#### Changed
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

#### Added (continued from 2026-03-16)
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

#### Changed
- **Memory: free old dataset before loading new one** — in `_on_open_dataset()`,
  all heavy `ctx` fields (`sdata`, `adata`, `clusterings`, `color_manager`,
  `transcript_loader`, layer references, etc.) are now explicitly set to `None`
  and `gc.collect()` is called after clearing napari layers (step 8b) and before
  `_load_dataset()` is called (step 9). This prevents peak RSS from reaching
  old-dataset + new-dataset simultaneously during a dataset switch.

### 2026-03-14

#### Added
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

### 2026-03-10

#### Fixed
- **Co-occurrence plot colors now match labels layer** — replaced seaborn `tab20`
  fallback in `make_co_occurrence_plot()` with `CLUSTER_PALETTE` from
  `utils/coloring.py`, ensuring line colors match the cell labels layer colors.

#### Changed
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

#### Added
- **ARMS tab code recording** — All ARMS operations (H&E load, landmark registration,
  GeoJSON/CSV load, tile DEG, DEG export, volcano plots) now emit `_record_code()`
  entries so they appear in the reproducible `code.py` journal.

### 2026-03-09

#### Fixed
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

#### Added
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

#### Changed
- `record_code` default changed from `False` to `True` (always-on recording).

### 2026-03-06

#### Added
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

#### Changed
- **Cluster label editor shared between tabs** — Gene Analysis and Clustering
  tabs now share the same `_build_label_editor_dialog()` helper, eliminating
  code duplication. Both editors use the per-clustering label storage.
- **String cluster ID support** — `add_clustering_to_obs()` in gene_analysis.py
  and `_get_cluster_ids_per_obs()` in the viewer now handle string-valued
  cluster IDs (from imported clusterings) via factorization fallback.

### 2026-03-05

#### Added
- **Cluster label editor in Clustering tab** — "Edit Cluster Labels..." button
  in the Clustering tab opens a dialog to rename clusters (manual cell type
  annotation) using the Cell Coloring tab's selected clustering.

#### Changed
- **Co-occurrence plot colors match napari** — line colors in co-occurrence plots
  now use the same palette as napari cell coloring (CLUSTER_PALETTE) instead of
  seaborn tab20. Falls back to tab20 for clusters without a stored color.
- **Clustering sync across tabs** — selecting a clustering in the Cell Coloring
  tab now auto-sets the same clustering in Gene Analysis, Ligand-Receptor,
  Nhood Enrichment, and Co-occurrence tabs (one-directional sync).

### 2026-03-04

#### Added
- **Leiden clustering tab** — new "Clustering" tab (Tab 0) with configurable
  `n_neighbors` (5–50), `n_pcs` (10–50), and `resolution` (0.1–5.0) parameters.
  Runs `sc.pp.neighbors` + `sc.tl.leiden` on a worker thread, stores results as
  `leiden_r{resolution}` in the clusterings dict, and refreshes all downstream
  ComboBoxes (Cell Coloring, Gene Analysis, L-R, Nhood, Co-occurrence).
  Reproducible code recording supported.

### 2026-03-03

#### Added
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

#### Changed
- Renamed `_get_lr_cluster_filter()` to `_get_cluster_filter()` — shared helper for
  cluster filtering in both L-R plot and nhood enrichment plot.

### 2026-03-02

#### Added
- **Session persistence** — viewer state (ROIs, H&E registration, analysis results,
  cluster labels) is automatically saved to `sdata_cached.zarr/` when the viewer closes.
  The H&E image is stored as a spatialdata multiscale image element (`images/he_image`)
  with its affine transformation, so the next launch restores the H&E overlay with
  registration already applied — no need to re-load or re-register. ROI polygons,
  cluster labels, and analysis results (rank genes, ROI DEG, L-R) are persisted in
  `viewer_session/` and restored on startup. Affine transformations are saved in
  real-time (on each registration/flip change). Skipped when `--no-cache` is used.
  New module `scripts/utils/session.py` with `save_session()` and `load_session()`.

#### Previously added
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

### 2026-03-01

#### Added
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

#### Changed
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

#### Previously added
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

#### Changed
- `load_sdata()` now skips `cells_boundaries` and `nucleus_boundaries` (were hidden
  anyway; 318K polygon shapes freeze napari). Saves ~30-40% of xenium() load time.
- `load_umap()` and `load_clusterings()` accept a `path` parameter instead of using
  module-level globals.
- Transcript cache default location changed from `scripts/transcript_cache/` to
  `data_dir/transcript_cache/`.

### 2026-02-28 — initial prototype (numbered 0.1.0, never published)

#### Added
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
