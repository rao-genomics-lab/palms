# TODO: user-configurable analysis templates — remaining work

Status: **Phases 0–4a implemented and committed** on `feat/template-registry`
(as of 2026-07-31). What follows is what is left. The plan this came from is
summarised in `CLAUDE.md` under "Code Recording"; the design paradigms are
Django template loaders + dpkg conffile upgrade policy + nf-core-style declared
contracts + dbt-style manifest stamping.

## Where things stand

| Phase | What it was | State |
|---|---|---|
| 0 | `EXECUTOR_BASE_NAMES`, kill the `$`-token hacks, Leiden declared output, mark hand-edited cells | done |
| 1 | Template text → `utils/step_templates/builtin/*.tmpl`, registry, registry-wide `check_step` gate, read-only Tools → Templates | done |
| 2 | Per-user overrides, validation gate, fail-to-builtin-loudly, provenance stamping, editable tab | done |
| 3 | `overrides.json` fork record, stale/needs-review detection, diff + take-new-default, notebook banner, `verify_notebook` templates section | done |
| 4a | `Preview(blocks, params, note)` providers across every tab; marker-plot override fix; registry-wide preview gates | done |
| 4b | Per-dataset scope, export/share, dry-run | **not started (optional)** |

821 passed, 24 skipped, ruff clean. Branch has 6 commits and **no PR yet**.

---

## 1. ~~Preview providers~~ — done

Twelve of fourteen templates now register a `Preview(blocks, params, note="")`
in `ctx.state["template_preview"]`, and each tab's own callback runs from that
same call, so the run and the preview are one expression. Providers carry the
**block selection** as well as the params — `_preview` used to render
`assemblies[0]` regardless, so the pane tracked the widgets in value but not in
shape.

Two declared exemptions in `tests/test_tab_templates.py::_NO_PROVIDER`:
`normalize` (no params) and `spatial_neighbors` (its `k` comes from whichever
tab called `ensure_spatial_neighbors`, so it uses the new `# sample-params:`
header field rather than a provider that would have to pick one slider).

Fixed along the way: **`genes.marker_plot` ignored user overrides entirely** —
it built its `Step` with `template=builtin_assemble(...)`, the shipped-files-only
path. Editing that template did nothing, and its nodes carried no `template_id`.
`test_every_step_resolves_user_overrides` now parses every `Step(...)` in every
tab.

## 2. Verification against a real dataset — done

Run 2026-07-31 against a 4104-cell / 477-gene Xenium output, by driving
`app.run_viewer` with `napari.run` replaced by an inspector that exits before
the `save_session` block. Nothing was written back: the zarr was untouched, and
the `viewer_cache/prov_graph.json` and `palms.log` the launch creates
were removed afterwards. **`analysis.py` in the dataset root was overwritten**
and could not be restored — it is derived and rewritten on every recorded step,
so any run of this kind destroys it. **Copy a dataset before probing it, or use
one with no `analysis.py`.**

- **Launch.** 7.2 s to ready. All 14 templates preview; all render code that
  compiles. 12 register a provider, and the two exemptions correctly show
  sample values — `spatial_neighbors` with `n_neighs=6` from `sample-params`,
  not the synthesised `1`.
- **The preview is live in value.** `genes.rank_genes` reads
  `groupby='graphclust', method='wilcoxon', n_genes=25` from the real widgets,
  where before this change it read `groupby='sample', method='sample',
  n_genes=1`.
- **Live in shape** is now a unit test
  (`test_toggling_a_real_widget_changes_the_previews_shape`), which drives the
  actual magicgui checkboxes — it needs a `QApplication`, not a dataset.
- **End-to-end.** Fork the `tail` block, re-run: 19 → 34 clusters, the recorded
  node's `code` contains the edit, `template_origin == "user+builtin"`, and the
  unedited run's node stays `builtin`. `customisation_banner` fires.
- **Marker plots honour overrides** (the fix in §1) — an override of the `head`
  block reaches the executed text with `origin="user+builtin"`, where the old
  `builtin_assemble` path demonstrably did not see it.
- **Failure paths on a real launch.** A syntax error and an undeclared `$token`,
  in two different templates: the viewer launched, warned **once per template**,
  badged both `✕ not used`, and both fell back to byte-identical builtin text.
  A third template was unaffected, and the rejected ones still preview from the
  builtin. `--no-user-templates` restores stock resolution.
- **`stock_templates`.** `verify_notebook.template_provenance` on the customised
  graph reports `stock_templates: false, n_customised: 1`, naming the node, its
  template id and hash.
- **Round-trip — the whole claim, on a customised session.** A real session was
  saved (Leiden with a forked `tail` block, persisted as `clustering_leiden_verify`,
  plus rank genes), then `scripts/verify_notebook.py <dataset> --out report.json`:

  ```
  Replayed in 26.2s
    ✓ leiden_verify: ARI = 1.000000 (28 clusters, 4104 cells)
    ✓ rank genes: top-10 identical in all 28 groups
    ⚠ 1 of 5 step(s) used customised templates
        clustering:leiden_verify  (user+builtin, clustering.leiden)
  ```

  Exit 0, `stock_templates: false`. The exported notebook carries
  `resolution=1.0 * 2` — the edit — so the customised template travelled into
  the notebook and reproduced the customised result exactly, which is the point:
  ARI 1.0 against a *stock* replay would have meant the override was ignored.

This run found one defect, now fixed: `verify_notebook.build_notebook` took its
cells from `prov_graph.graph_to_cells` (the one carrying `node_id`, which is what
turns nbclient's "cell 4 failed" into a named step) and so never got the
customisation banner that `notebook_export.graph_to_cells` adds. The `--work-dir`
notebook — the artifact most likely to be kept or forwarded — was the only output
of a customised session that did not say so. Fixed by prepending it explicitly,
with a matching `None` in `node_ids`: `node_of()` indexes that list by absolute
cell index, so inserting a cell without shifting it would have misattributed
every timing in the report.

## 3. Phase 4b — optional, not started

Only relevant if the per-user-only scope decision changes:

- **Per-dataset override scope**, so a customised analysis travels with the data.
  Site it at `<data_path>/palms_templates/` — *outside* `viewer_cache/`,
  whose documented invariant (`store_inventory.py:672-674`) is that everything
  in it is deletable viewer output. Otherwise `_BLOCKED_SIDECAR_PREFIX`
  (`store_inventory.py:158-166`) has to be generalised to `(prefix, reason)`
  pairs. Note this is where the upgrade problem multiplies: an override can
  arrive on a machine running a different app version, where nothing knows it is
  stale.
- **Template export / share** for a lab.
- **Dry-run on a subsample** — execute a step against `adata[::100]` in a scratch
  namespace via the Notebook tab's kernel, report exceptions, discard results.

## 4. Loose ends

- **No PR** for `feat/template-registry`.
- **Galaxy's tool-XML and tool-versioning docs were never read.** The design came
  from the dpkg conffile model instead, which covered the upgrade case well. It
  becomes relevant again only for Phase 4b's sharing work — Galaxy has long
  experience with tool definitions travelling between installations.
- **Testing blind spot to keep in mind.** `tests/conftest.py` sets
  `PALMS_TEMPLATE_PATH` for the whole suite, so the branch taken when it
  is *unset* — every real user — runs nowhere by default. That gap shipped a
  viewer that could not start. `tests/test_template_overrides.py::no_env` now
  covers it; anything new that is reachable at launch belongs there too.
