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

## 2. Verification that needs a real dataset

Plan steps 1–3 (unit, migration equivalence, notebook replay) are done and green.
Steps 4–6 have **not** been run, because they need a dataset. Given that a
launch-blocking recursion bug shipped past 807 passing tests, these matter:

- **The preview is live in shape, not only in value.** Open Tools → Templates →
  `clustering.leiden`, then toggle "Use HVGs" and "Scale" in the Clustering tab
  and confirm the *blocks* appear and vanish in the pane. Unit tests assert this
  through a stub provider; nothing yet observes it through real Qt widgets in a
  running viewer. Spot-check `roi.deg`, `genes.marker_plot` and
  `genes.cnv_infercnv` for live values.
- **End-to-end.** Launch `xenium-viewer <dataset>`, run Leiden, confirm the
  Templates preview matches the Notebook tab's `clustering:leiden_*` cell. Then
  edit the `hvg` block, Validate, Save, re-run Leiden, and confirm **(a)** the
  GUI result changes, **(b)** the Notebook cell shows the edited source, **(c)**
  the node is stamped `origin="user"`.
- **Marker plots honour overrides** (the fix in §1). Edit `genes.marker_plot`'s
  `head` block, Save, generate a dotplot, and confirm the Notebook cell shows
  the edit and the node carries `template_id="genes.marker_plot"`,
  `origin="user"`.
- **Round-trip.** `scripts/verify_notebook.py <dataset> --out report.json` on
  that customised session: it must replay against the raw Xenium output, and the
  report must carry `"stock_templates": false` naming the customised steps.
- **Failure paths on a real launch.** Hand-write a broken override (syntax
  error; undeclared `$X`; missing required param) into
  `~/.config/xenium-viewer/templates/`, then confirm the viewer **launches**,
  warns once, uses the builtin, and badges the template `✕ not used`. Then
  confirm `--no-user-templates` restores stock behaviour.

## 3. Phase 4b — optional, not started

Only relevant if the per-user-only scope decision changes:

- **Per-dataset override scope**, so a customised analysis travels with the data.
  Site it at `<data_path>/xenium_viewer_templates/` — *outside* `viewer_cache/`,
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
  `XENIUM_VIEWER_TEMPLATE_PATH` for the whole suite, so the branch taken when it
  is *unset* — every real user — runs nowhere by default. That gap shipped a
  viewer that could not start. `tests/test_template_overrides.py::no_env` now
  covers it; anything new that is reachable at launch belongs there too.
