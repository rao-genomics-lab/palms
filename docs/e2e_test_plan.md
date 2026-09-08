# End-to-end test

## Why

The project's standing position was recorded in `tests/test_units.py`: *"the napari
GUI proper has no automated coverage"*. In CI that was literally true. Two tests
built a real `napari.Viewer` and one of them was gated behind `requires_display`,
which is false under `QT_QPA_PLATFORM=offscreen` — so **no CI run ever constructed
a viewer**. Each of the 29 tabs was built only against a `SimpleNamespace` stub,
never against a real `ViewerContext`, and never alongside the other 28. And the
save-on-exit block was straight-line code *after* `napari.run()`, so the last thing
every session does could not be reached without an event loop.

`tests/e2e/` closes that, over one lifecycle:

```text
staged dataset -> create_viewer -> operate through the real widgets
              -> shutdown_viewer -> re-read the zarr -> create_viewer again
```

## What it covers, and what it does not

Covered: the load sequence, all 29 tabs built from one real context, both docks,
the layer set, a recorded analysis step driven through its own widget, the QC
filter (the one action that rebinds `ctx.adata`), the exit-time persistence path,
and the derived outputs (`analysis.py`, `analysis_notebook.ipynb`).

**Not covered: the raw `spatialdata_io.xenium` load path.** The fixture is a crop
export, so `has_raw_xenium_source()` is False and the zarr *is* the data. That
makes it a fixture for the never-stale / never-rebuild branch, and it means a
green E2E says nothing about a first launch on raw 10x output. Closing that needs
raw Xenium output no demo dataset has.

Also not covered: `PALMS_TEMPLATE_PATH` is emptied by `tests/conftest.py` for the
whole suite, so the user-template resolution path a real user takes runs nowhere
here either. `tests/test_template_overrides.py` is where that is covered.

## The fixture

`tests/data/crop_7/` — 6.7 MB, committed, so CI needs no download and no network.
`tests/data/README.md` records what it is, how it was made
(`scripts/prepare_e2e_fixture.py`) and why each of the three preparation steps
exists. Point `PALMS_E2E_DATASET` at another Xenium output directory to run the
same tests against it; the assertions that name this fixture's own numbers skip.

Two properties of the fixture are load-bearing and are enforced rather than
documented:

- **The source is never launched against.** `e2e_dataset` stages a copy into
  `tmp_path` and hashes the source before and after — a launch writes `palms.log`,
  repairs the cache, and rewrites `analysis.py` and `viewer_cache/prov_graph.json`
  beside whatever it is pointed at.
- **The staged copy is path-repaired.** The provenance graph records absolute
  paths, and `app.py` re-emits `preamble` for the current `data_path` on every
  launch — so a bare `cp` flags all 30 descendants stale and every staleness
  assertion becomes vacuous. The fixture ships recording the placeholder
  `/palms-e2e-fixture/crop_7`; staging rewrites it with the shipped
  `rename_dataset.repair`.

## The seams

`run_viewer()` was split, with no behaviour change:

| | |
|---|---|
| `app.create_viewer(...)` | everything before `napari.run()`; returns `(viewer, ctx, app_state)` |
| `app.shutdown_viewer(app_state, ...)` | everything after it |
| `app.run_viewer(...)` | the two, with `napari.run()` between |

`app_state["ctx"]` is the authoritative context — a dataset switch rebinds it, and
`shutdown_viewer` reads it from there so a switch does not persist the dataset the
user left. `app_state["snapshot_layers"]` is the `aboutToQuit` hook, published so a
caller with no event loop invokes *that function* rather than a copy: the layers
must be read while their Qt objects are alive and the store written after the loop
stops, and that ordering is the whole reason the two halves are separate.

`scripts/capture_screenshots.py` calls `create_viewer` too. It previously called
`_do_full_init` with a hand-rolled `_app` dict missing `plots_dock` and
`plots_panel` — the drift a shared seam removes, and `tests/test_rig_is_shared.py`
is the guard.

## The rig

`src/palms/testing/rig.py` — `Rig` navigates to a tab by its visible labels, finds
a magicgui widget by the label the user reads (via `native._magic_widget`), clicks
it and waits for the worker it started. It came out of
`scripts/capture_screenshots.py`, which still uses it; two copies would disagree
the first time either learned something about a widget the other did not.

`strict` is the one behavioural difference between the callers: a screenshot of a
step that is still running is still a picture of that step, so the capture script
passes `strict=False` and a `wait_idle` timeout only prints. A test takes the
default, where it raises.

## Three guards, each closing a way to lose a CI run rather than fail one

- **`no_dialogs`** makes any modal raise. `reporting._headless()` suppresses the
  modals *it* raises, but nothing suppresses a `QFileDialog` or a `QMessageBox`
  from a tab, and one of those blocks the runner until the job times out.
- **`deadline`** is `faulthandler.dump_traceback_later(900, exit=True)` — stdlib,
  so no plugin, and unlike `pytest-timeout` it prints where the process actually
  was, which is the diagnostic a Qt hang otherwise denies.
- **`Session.settle()`** runs the event loop until `WorkerBase._worker_set` is
  empty, both after a launch and before a close. Not tidiness: a worker whose
  `returned` callback lands after Qt has deleted the widgets it writes into raises
  inside a slot, and **PyQt6 turns an exception in a slot into `qFatal()`** — a
  core dump with no test name.

That last one is why writing this suite found two real defects, both on the
teardown path that `remove_dock_widget` also takes on a live dataset switch:

- `tab_qc._refresh_readout` is connected to four widgets' `changed`, and
  destroying the tree emits from one after its siblings are gone.
- `tab_he_registration`'s restore worker guards on `ctx.dataset_generation`, which
  a *close* does not increment.

## Does it catch anything?

A test that cannot be made to fail is not yet a test, so three defects were
injected and the suite re-run:

| injected | result |
|---|---|
| drop `ctx.ensure_qc_filter()` from `_do_full_init` | **caught** — *"the filter was stored but not re-applied", `assert 151 == 120`* |
| drop `save_session` from `shutdown_viewer` | **caught** — *`assert None == {'min_counts': 60, 'min_cells': 3}`* |
| drop `_persist_table` from `shutdown_viewer` | **not caught**, and correctly so |

The third is worth writing down rather than fixing. Analysis artifacts are
written *eagerly* — `save_clustering_to_adata` puts the column in the store when
the step runs — so `_persist_table` at shutdown is a backstop, not the writer,
and a Leiden run reaches disk without it. Any future test that means to pin the
shutdown table write has to produce a change nothing else persists.

## Running it

```bash
pytest tests/e2e -vv --durations=10      # no env vars: conftest.py selects offscreen
PALMS_E2E_DATASET=/path/to/dataset pytest tests/e2e
```

CI runs it in the existing Linux conda job — the fixture is in-repo, so no `-m`
selection and no workflow change beyond `--durations=10`. The suite is marked
`e2e` (registered in `pyproject.toml`) for local selection, and skipped on macOS,
which runs under cocoa against a real window server; enable it there once the
Linux leg has been green for a while.

**If a runner cannot give vispy a usable context**, the fix is a display, not a
redesign: wrap the `pytest` step in `xvfb-run -a`. Offscreen is what
`tests/conftest.py` already documents as the configuration CI needs, so try it
first.

## Follow-ups

Keep these independent, so a failure names the boundary that regressed.

- ROI draw → persist → reload; H&E landmark registration; Tools → Dataset deletion.
- Crop-export-then-open. crop_7 *is* a crop, so cropping it and opening the result
  covers `crop_export` and the cache-only load path together.
- The negative case for `app._load_prov_graph_items`: truncate
  `viewer_cache/prov_graph.json` and assert the launch warns and restores from the
  session attr.
- The un-repaired-copy case: stage without `repair` and assert every descendant is
  flagged stale — the rename hazard, currently covered only at unit level.
- A dataset *switch* (not just a relaunch), which is the path both teardown
  defects above actually reach in a live session.
