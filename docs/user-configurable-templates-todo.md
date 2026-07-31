# TODO: user-configurable analysis templates — remaining work

Status: **Phases 0–3 implemented and committed** on `feat/template-registry`
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
| 4 | Per-dataset scope, export/share, dry-run | **not started (optional)** |

813 passed, 24 skipped, ruff clean. Branch has 5 commits and **no PR yet**.

---

## 1. Preview providers — only 1 of 14 tabs has one

**The main functional gap, and the most visible one.**

Tools → Templates shows a live preview of the exact string a step would `exec`.
It is real only for Leiden, because `tab_clustering.py` is the only tab that
registers a provider:

```python
ctx.state.setdefault("template_preview_params", {})[TEMPLATE_ID] = _leiden_params
```

Every other template falls back to synthesised literals of the right type, so
the pane reads:

```
# preview — sample values
# Rank genes: groupby='sample', method='sample', n_genes=1
```

That is labelled and working as designed, but it makes the tab's headline
question — *what will this button actually run?* — honest for one template and
merely illustrative for thirteen.

**The fix, per tab** (~20 lines each, `tab_clustering.py:156-181` is the worked
example): extract the params dict the callback already builds into a named
closure, register it, and have the callback call it. The point is that the run
and the preview then read **one** expression of "the current settings" — a
second copy is exactly the drift `Step` exists to rule out.

Tabs needing it: `tab_gene_analysis`, `tab_nhood`, `tab_co_occurrence`,
`tab_ligrec`, `tab_marker_genes`, `tab_gene_correlation`, `tab_roi` (×3),
`tab_cnv`, plus `normalize` / `spatial_neighbors` in `_helpers.py` (these two
have no widgets, so `sample_params` in the `.tmpl` header may be the better
answer for them than a provider).

`tests/test_tab_templates.py::test_the_clustering_tab_registers_a_preview_provider`
is the pattern for guarding each one.

## 2. Verification that needs a real dataset

Plan steps 1–3 (unit, migration equivalence, notebook replay) are done and green.
Steps 4–6 have **not** been run, because they need a dataset. Given that a
launch-blocking recursion bug shipped past 807 passing tests, these matter:

- **End-to-end.** Launch `xenium-viewer <dataset>`, run Leiden, confirm the
  Templates preview matches the Notebook tab's `clustering:leiden_*` cell. Then
  edit the `hvg` block, Validate, Save, re-run Leiden, and confirm **(a)** the
  GUI result changes, **(b)** the Notebook cell shows the edited source, **(c)**
  the node is stamped `origin="user"`.
- **Round-trip.** `scripts/verify_notebook.py <dataset> --out report.json` on
  that customised session: it must replay against the raw Xenium output, and the
  report must carry `"stock_templates": false` naming the customised step.
- **Failure paths on a real launch.** Hand-write a broken override (syntax
  error; undeclared `$X`; missing required param) into
  `~/.config/xenium-viewer/templates/`, then confirm the viewer **launches**,
  warns once, uses the builtin, and badges the template `✕ not used`. Then
  confirm `--no-user-templates` restores stock behaviour.

## 3. Phase 4 — optional, not started

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
  becomes relevant again only for Phase 4's sharing work — Galaxy has long
  experience with tool definitions travelling between installations.
- **Testing blind spot to keep in mind.** `tests/conftest.py` sets
  `XENIUM_VIEWER_TEMPLATE_PATH` for the whole suite, so the branch taken when it
  is *unset* — every real user — runs nowhere by default. That gap shipped a
  viewer that could not start. `tests/test_template_overrides.py::no_env` now
  covers it; anything new that is reachable at launch belongs there too.
