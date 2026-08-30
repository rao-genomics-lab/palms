# Design note: user-authored analyses, and an LLM that helps write them

Status: **design only — nothing here is implemented.** Written 2026-08-26 against
`main` (`6a0493b`). Continues `user-configurable-templates-todo.md`, whose Phase 4b
this partly subsumes.

The question this answers: can a user add a *new* analysis and a GUI element for
it, with the app linting the analysis and generating tests for it — and can a small
model fine-tuned on spatial-omics APIs guide them through doing so.

**Yes, and most of the machinery already exists.** But not the part it looks like.
The template registry gives us a declarative contract, a static validator, a
registry-wide test gate, an editor, provenance stamping and notebook replay. What
it does not give us is a way for a template id to **exist** that this package did
not ship. The `~/.config/palms/templates/` directory is a *shadowing*
mechanism, not a plugin mechanism — a subtlety that is easy to miss, because
everything about it reads like a plugin directory.

---

## 1. What already exists

| The ask | What covers it today |
|---|---|
| A declarative definition of an analysis | The `.tmpl` header — `id`, `params`, `requires`, `outputs`, `assemblies`, `frozen-blocks`, `sample-params` — parsed by `loader.parse_template` into a `TemplateSpec` (`utils/step_templates/spec.py`) |
| Lint | `utils/step_templates/validate.py::validate` — **fourteen distinct findings, thirteen of them errors that block activation**, all static: render with `synth_params()`, `ast.parse`, walk. Nothing is executed, which is why the Validate button answers in milliseconds instead of after a ten-minute run |
| Tests | `tests/test_template_registry.py` — parametrized over `builtin_ids() × assemblies`, **seven properties per rendering**. A template that reaches the registry inherits all of them for free |
| Safe execution | `utils/steps.py::StepExecutor.run` — literal-only params (`ast.literal_eval(repr(v)) == v`), statement-at-a-time `exec` with line numbers preserved, records **only on success**, raises if a declared output was never bound |
| Reproducibility | `utils/notebook_export.py` + `scripts/verify_notebook.py`. `template_origin` / `template_hash` already mark a non-stock run; `customisation_banner` puts it at the top of the exported notebook |
| An editor GUI | `tabs/tab_templates.py` — Default \| Yours panes, Validate / Save / Revert, a live preview that is `Step.render()` rather than a reconstruction, `overrides.json` fork tracking, stale-block diff and "Take new default" |
| LLM plumbing | `utils/gene_analysis.py:304-369` — `build_llm_annotation_prompt` / `parse_llm_annotation_response` / `run_llm_annotation`: subprocess to an installed `claude` / `gemini` / `codex`, tolerant JSON extraction, **no API keys anywhere**, driven from a `thread_worker` with `errored.connect` |

Two facts here are counter-intuitive enough to state plainly, because both change
what the remaining work is.

**The builtin registry is already a directory scan.** `loader._builtin_registry()`
(`loader.py:326`) iterates `resources.files(_BUILTIN_PACKAGE)` for `*.tmpl` and
raises if two files claim one id. Adding a *shipped* template is a file, not a code
change. So the barrier to a new analysis is not "the registry is a hardcoded list"
— it isn't one.

**The seven registry gates are properties, not per-analysis tests.** They assert
that every rendering parses, is self-contained against
`EXECUTOR_BASE_NAMES | requires`, binds its declared outputs, never imports
`palms`, uses every declared param, mentions every *required* param in
every assembly, and contains no comment-only block. That is most of what "generate
tests for the new analysis" means, and it is already written. What a new analysis
needs *generated* is the one thing those gates deliberately do not do: **run it**.

---

## 2. The four things that are missing

### (a) A user template id cannot exist

Three places, all in `utils/step_templates/loader.py`:

1. `resolve()` (line 492) opens with `base = builtin_spec(template_id)`.
   `builtin_spec` (line 346) raises `TemplateError("unknown template …")` for an id
   with no shipped file. So `resolve("custom.my_analysis")` raises rather than
   loading a user file.
2. `_override_files()` (line 434) builds **one exact path** per directory,
   `directory / f"{template_id}.tmpl"`. There is no `glob`/`iterdir` over
   `search_path()` anywhere in the package — a user file is only ever *looked for
   by a known id*, never *discovered*.
3. `tab_templates._populate()` (line 275) iterates `builtin_ids()`, so a user-only
   `.tmpl` would not appear in the tree even if the loader could read it.

**Proposed:** `user_ids()` (scan `search_path()` for `*.tmpl`, read only the
`# id:` header), `registry_ids()` (builtin ∪ user, builtin winning on collision),
and `resolve()` tolerating a user-only base — validate with `builtin=None`, which
`validate()` already supports and documents as "omit it to validate a builtin
against itself".

**The one place the safety model genuinely changes.** Today `_merge` (line 451)
takes only *block text* from an override and keeps the builtin's `params`,
`requires`, `outputs` and `assemblies`: "it is what the call site and the viewer
agree on, not something a template file gets to redefine." A bad header in an
override is therefore harmless — it cannot even be expressed. For a user-authored
template **the header is the contract**, with nothing shipped standing behind it.
Two consequences, and they are the design's load-bearing pair:

- A user-only template must pass the **full** validator, not the override subset,
  and is **inert until it does** — listed in the tree, dimmed, with its problems
  shown. Never silently falling back, because there is nothing to fall back to.
- `_check_structure`'s three checks (schema version, unknown blocks, frozen blocks)
  are meaningless without a builtin and are correctly skipped. That leaves the
  eleven that matter, including the one `validate.py`'s docstring singles out: a
  required param the template never mentions, which would run, succeed, and ignore
  the user's setting.

### (b) No GUI exists for a template with no call site

Every run site hardcodes its id as a module constant — twelve in `tabs/`, plus
`normalize` and `spatial_neighbors` bound at import in `tabs/_helpers.py:43,51`. A
template with no call site has no button. And `app.py` builds the panel from a
hardcoded list repeated **four times**: imports (360-385), phased construction
(388-434), `addTab` (482-526), and `all_exports` (529-536). There is no tab
registry, no entry point, no discovery of any kind — `pyproject.toml` declares only
console scripts.

**Proposed: one** new tab, `tabs/tab_analyses.py`, generic rather than generated.
It does not write Qt code; it *renders* a form from the spec:

- a `ParamSpec.type` → magicgui widget table, mirroring `spec._SYNTH`'s type list
  (`str int float bool list dict tuple`);
- a combo over `spec.assemblies` for block selection;
- a Run button that builds
  `Step(id=…, **step_template(tid, blocks), params=form_values)` and runs it in a
  `thread_worker` — with a failure path, which most but not all tabs have. Eleven
  connect `worker.errored`; `tab_nhood.py:114` instead catches `StepError` inside
  the worker and returns a payload carrying the message. `tab_clustering.py` does
  neither, so a failing step there leaves the Run button disabled with no
  explanation. A user-authored analysis fails far more often than a shipped one,
  so this tab must handle it explicitly;
- and a `Preview(blocks, params)` registered into `ctx.state["template_preview"]`.

That last line is the whole trick. `ctx.state["template_preview"]` is the one
genuinely dynamic registry in the codebase — a dict tabs populate with `setdefault`
and the Templates tab consumes. Registering into it means **the Templates tab's
preview pane works for user analyses with no new preview machinery at all**, and
`tab_templates._preview` needs no change.

**Block selection is the real design cost, and the note should not soften it.** At
every existing call site the branch structure *is* what the widgets mean
(`if use_hvg or do_scale` → the `pca` block), and keeping that in Python is
deliberate: "the branch structure *is* what the widgets mean, while the statements
are what someone actually wants to change" (`spec.py`). A generic form cannot infer
it. So a user template gets a **chosen** assembly rather than a derived one, and a
template that wants real branching declares one assembly per branch. This is the
price of the declarative-only scope, and it is the right price: the alternative is
importing user-written Qt into the process, which puts arbitrary code in the widget
tree with no static gate on it whatsoever.

**Dependencies need a home.** `spec.requires` names which *values* must be bound,
not which *step* binds them — and the prov-graph needs node ids, because
`ProvGraph.upsert` errors on a missing dependency at record time. Proposed: an
additive `# deps:` header field, with a default table for the names the app
already produces (`adata_norm` → `normalize`, and call `ctx.ensure_normalized()`
first; `roi_polygons` → `rois`).

One hazard to record: `parse_template` reads only the header fields it knows and
**silently ignores the rest**. A `.tmpl` carrying `# deps:` opened by an older
viewer loses its dependencies without a word — the analysis would still run and
still record, but as a root node, and the notebook's topological sort would be
wrong. That is an argument for bumping `schema-version` when the field lands:
`_check_structure`'s "written for a newer release; upgrade rather than guessing"
is exactly the message this case wants.

### (c) Test generation — two tiers, mirroring what the repo already does

- **In-app "Check".** Run the seven registry properties in-process against the
  user template and show them in the existing problems list. Free: same functions,
  no new assertions.
- **A generated `test_<id>.py`** written beside the template, using the
  session-scoped `replay_adata` fixture from `tests/conftest.py` (200 cells ×
  60 genes, two seeded populations, two spatial blobs). **Templated, not
  LLM-authored** — the properties are known, so rendering them from a jig is both
  cheaper and more trustworthy than asking a model to invent them. The user can
  run it with plain `pytest`.
- **Dry-run execution** against `adata[::N]` in a scratch namespace. This is the
  only check that catches a template that lints clean and raises at run time, and
  it is the one thing the seven gates cannot do. It is already scoped as Phase 4b
  in `user-configurable-templates-todo.md` — **for user-authored code it stops
  being optional and becomes the core of the feature.**
- **Notebook replay**, opt-in. `notebook_export.execute_notebook` in a clean
  kernel is the real proof, but it is ~35 s and needs `nbclient`/`ipykernel`
  (the `test` extra). It belongs behind a button, not in the save path.

### (d) `EXECUTOR_BASE_NAMES` needs no change

Worth a paragraph because it is what everyone assumes is the blocker. The set is
exactly-checked in both directions — `namespace.check_base_namespace` is called on
the dict `_helpers._get_executor` actually builds, so "the set validation checks
against and the set execution provides cannot drift." Widening it for user analyses
would be the wrong move twice over: it would weaken that check, and the notebook
preamble binds the same names, so an extra one replays as a `NameError` in a clean
kernel.

A template that needs `scipy` simply imports it — `genes.correlation`'s `head`
block already does exactly this, and `validate.module_level_bindings` counts
imports as bindings. Leave it alone.

---

## 3. Constraints the design must not break

Each is enforced by an existing test, so each is a hard rule rather than a
preference.

- **`builtin_*` must keep never seeing an override path.** Six template-pinning
  test modules depend on it, and it is what makes them immune to a developer's own
  config. It is also what silently disabled customisation for `genes.marker_plot`
  for a whole phase — the one run site that used `builtin_assemble`. Run sites use
  `step_template`; `tests/test_tab_templates.py::test_every_step_resolves_user_overrides`
  parses every `Step(...)` in every tab.
- **`tests/conftest.py` blanks `PALMS_TEMPLATE_PATH` suite-wide**, so the
  user-scope path — which is every real user — runs nowhere by default. That gap
  once shipped a viewer that could not start (`search_path()` and
  `user_template_dir()` recursed into each other). Everything in §2(a) is reachable
  at launch, so it belongs in `tests/test_template_overrides.py::no_env`, not in a
  module that inherits the blanked variable.
- **A new GUI tab needs its documentation in the same commit.** `TAB_PAGES` in
  `tests/test_docs_links.py` reconciles the abbreviated tab label against a
  `Tab-*.md` page, and `test_every_published_page_is_in_the_mkdocs_nav` /
  `test_every_published_page_is_linked_from_the_sidebar` require a nav entry and a
  sidebar link. A tab with no page turns CI red — deliberately.
- **`build_tab` must return a `make_tab()` / `scrollable()`-wrapped widget**, and
  a button row outside the scroll area must use `toolbar_row()`. A stacked widget's
  minimum size is the maximum over all its pages, hidden ones included, so one
  unwrapped page becomes the floor for the whole dock. `tests/test_control_panel_size.py`
  measures it.
- **The new tab must be added to `all_exports`** in `app.py`, or its
  `restore_session` silently never runs. That list is the fourth hardcoded copy of
  the tab set and the easiest of the four to forget.
- **`validate.py` is not a security boundary, and this feature must not be sold as
  though it were.** Its own docstring says so, and the reasoning is right: a user
  can already execute arbitrary Python in the Notebook tab's free-form cells, and
  the notebook this system exports is a file of code they are meant to run. A user
  analysis is `exec`'d at the same trust level. Saying "the app lints it" invites
  the wrong kind of review and gives false assurance. The lint is a **correctness**
  gate aimed at code that runs and quietly produces the wrong number.

---

## 4. Staging

| Stage | Content | Touches |
|---|---|---|
| A | User-scope ids: `user_ids()`, `registry_ids()`, `resolve()` with a user-only base, full-validator gate, "New template…" in the Templates tab | `loader.py`, `overrides.py`, `tab_templates.py` |
| B | `tabs/tab_analyses.py` — generated form, assembly combo, Run, `Preview` registration, `# deps:` field | new tab, `app.py` ×4, `Tab-Analyses.md` + nav + sidebar + `TAB_PAGES` |
| C | Result kinds and persistence — **split into C1–C4 in §7**, which is where most of the remaining work turned out to be | `adata_persistence.py`, `loader._USER_OBS_PREFIXES`, `tab_cell_coloring.py`, new `step_templates/scaffold.py` |
| D | LLM authoring assistant (§5) | new `utils/analysis_assistant.py` |

All of it behind one experimental flag, `--experimental-analyses`, off by default —
following the `--no-user-templates` precedent that a process-wide switch "must not
require editing files."

Stages A–C are useful on their own: they are the manual authoring path, and the
LLM in stage D only writes files that path already accepts.

---

## 5. The LLM half

### The validator is what makes any model usable here

Fourteen findings, all static, all offline, all sub-millisecond. That turns
generation into **rejection sampling**: generate → validate → feed the `Problem`
messages back → retry → dry-run → accept. And `Problem.__str__` is already written
as remediation prose aimed at a human —

> `assembly head+tail reads ['adata_nrom'], which nothing provides. Either a
> dependency step binds them — add them to '# requires:' — or this is a typo that
> would replay as a NameError.`

— so the error channel a repair loop needs exists and is well phrased for it. **A
weak model inside this loop beats a strong one outside it**, which is the fact that
makes a small local model a reasonable target at all.

### Why fine-tuning is the last stage, not the first

The corpus is the whole problem. What the repo has: **fourteen** templates. What a
7B-class code model needs to internalise a DSL and its API conventions: on the
order of **1–5k** instruction→template pairs. Nothing closes that gap by hand.

What this codebase does have, and what is genuinely unusual: a **free, exact,
automated reward signal**. A generation that passes `validate` *and* executes as a
dry-run *and* replays through the notebook is verified — no human labelling, no
judge model, no rubric. That makes the bootstrap below produce training data as a
by-product of the feature working.

**Stage 0 — prompt + retrieval.** Reuse the `run_llm_annotation` shell-out
pattern: subprocess to an installed CLI, no API key, no new dependency, degrade to
a status-bar message if the binary is absent. Context to supply: the header
grammar, `EXECUTOR_BASE_NAMES`, rule (e) from `CLAUDE.md` (prefer a library API
over a hand-rolled numpy/pandas equivalent — `sc.get.obs_df` over indexing `.X`,
`sd.polygon_query` over a point-in-polygon loop, `shapely.make_valid` over
`buffer(0)`), and the two or three shipped templates nearest the request.

One gap to close first: `run_llm_annotation` has **no availability probe** — a
missing CLI surfaces only as a `FileNotFoundError` at click time. An assistant that
is the primary path for a feature should `shutil.which` up front and say so.

**Stage 1 — harvest.** Log every (instruction, generation, verdict, repair round)
tuple. Accepted pairs are the supervised corpus. Rejected ones *with their
`Problem` lists and the accepted repair* are the more valuable half: repairing a
template against a named validator error is a better-specified task than writing
one from scratch, and it is exactly what the interactive loop does.

**Stage 2 — LoRA fine-tune** a small open code model on the harvested set,
evaluated against a **held-out benchmark whose metric is validator-pass +
dry-run-executes + notebook-replays**. That metric is objective, already
implemented, and needs no reference answer — rare for a code-generation eval, and
the reason this is worth doing rather than just prompting forever.

**Stage 3 — serve locally** (ollama / llama.cpp) so the feature works with no
account and no network, which is the actual argument for a small model over a
frontier one. Not speed, and not quality: *availability offline*, in a viewer whose
whole point is being the open alternative that runs where the vendor tool does not.

### Hardware reality, as measured on this box

- No NVIDIA GPU (`nvidia-smi` absent), no `ollama`.
- `claude`, `gemini` and `codex` are all on `PATH`.

So: **stage 0 works today with no new dependency**; the stage-2 fine-tune has to
happen on rented GPU; and stage-3 CPU inference of a quantised 7B is workable for a
~300-token template but will not feel interactive. Worth recording rather than
assuming a GPU that isn't there.

### The failure the model cannot be trusted to avoid

An analysis that is syntactically perfect, passes all fourteen checks, executes
cleanly — and is scientifically wrong. Wrong test for the design, unnormalised
input, a neighbourhood radius that means nothing at this resolution. `validate.py`
names this as the failure that matters most, and no static check reaches it.

So the position is: **the assistant produces a draft for a human to read**, and the
app's job is to make sure nobody downstream mistakes it for a shipped method. That
part is already built, which is the quiet strength of the existing design — a
generated template has `template_origin != "builtin"`, so it already:

- prepends `notebook_export.customisation_banner` to the exported notebook,
- reports `stock_templates: false` in `scripts/verify_notebook.py`'s report,
  naming the node, its template id and its hash,
- and badges in the Templates tab.

An LLM-authored analysis therefore arrives in a reader's hands already labelled as
not-ours, by machinery that exists for a different reason and happens to be exactly
right for this one. Any implementation of stage D must keep that true — a generated
template that stamps as `builtin` would be the single worst outcome of this whole
feature.

---

## 6. Worked example: spatially variable genes (`sq.gr.spatial_autocorr`)

**This transcript is measured, not imagined.** Every validator message, every
rendering and the executed result below were produced by running the real
`parse_template` / `validate` / `StepExecutor` / `notebook_export` against a
synthetic AnnData shaped like the `replay_adata` fixture (200 cells × 60 genes,
two seeded populations, two spatial blobs) in the `palms` env on
2026-08-26. Only the GUI mock-ups are drawn by hand.

### The gap

`sq.gr.spatial_autocorr` — Moran's I / Geary's C — ships in the pinned squidpy
(1.8.2) and appears **nowhere** in `src/`, `tests/` or `docs/`. It answers a
question the viewer currently cannot: *which genes are spatially structured at
all*, independent of any clustering. Every spatial tab we have starts from a
cluster label; this one does not, which is exactly why someone would want it and
why it never got built.

It is also an unusually good stress test of the design, because it needs all four
of the pieces §2 says are missing:

- it depends on a step the app already provides (`spatial_neighbors`, built on
  `adata_norm`) — so it exercises `requires` **and** `deps`;
- it takes real parameters (`mode`, `n_perms`, a result cut-off);
- it returns a **table**, not an obs column — so nothing in the existing
  clustering plumbing applies;
- and "significant genes only" is a genuine optional step — so it needs two
  assemblies.

### Beat 1 — a new template

Tools → Templates → **New analysis…**, id `spatial.autocorr`. The editor opens on
a skeleton of the header fields with an empty `main` block. Nothing is active yet:
a user-only template is inert until it validates (§2a).

### Beat 2 — the assistant drafts it

Retrieval picks the two nearest shipped templates by id prefix and shared
`requires` — `spatial_neighbors.tmpl` and `spatial.nhood.tmpl` — and the prompt
carries the header grammar, `EXECUTOR_BASE_NAMES`, and rule (e) from `CLAUDE.md`.
The user's instruction is one sentence: *"rank genes by Moran's I on the spatial
graph, with permutation p-values, and let me keep only the significant ones."*

### Beat 3 — the first draft fails, and the message is the fix

The draft declares `n_perms` but hardcodes `100` in the call — the single most
consequential mistake this validator exists to catch, because the analysis would
run, succeed, and ignore the user's setting:

```
warning: declares params ['n_perms'] that no block uses
error: assembly main never uses required param(s) ['n_perms'], so those settings
       would be silently ignored. Mark them optional with '?' in the header if
       that is genuinely intended.
```

Note what is *not* reported. `_check_assembly` returns on the first cause, so the
draft's two other defects stay hidden this round — by design, the same reasoning
as `_check_structure`'s early return: report the cause, not the cascade.

### Beat 4 — round two reveals the rest

With `n_perms=$n_perms` substituted, the next two surface together:

```
error: assembly main reads ['adata_norm'], which nothing provides. Either a
       dependency step binds them — add them to '# requires:' — or this is a
       typo that would replay as a NameError.
error: assembly main declares output(s) ['autocorr_df'] that it never binds;
       the viewer reads results by those names.
```

The first is the important one and it is worth dwelling on: `adata` *is* in the
base namespace, so a template can read it with no complaint — but the k-NN graph
is built on `adata_norm` (`spatial_neighbors.tmpl`), and `adata` has no `.obsp`.
A draft written against `adata` lints clean and fails at run time. This one failed
at lint only because it declared the output it meant to produce and reached for
`adata_norm` by name. **The lint catches the honest mistake; the dry-run in Beat 8
is what catches the plausible one.**

### Beat 5 — clean

```python
# palms template
# id: spatial.autocorr
# schema-version: 1
# assemblies: main+tail | main+significant+tail
# doc: Spatially variable genes by Moran's I / Geary's C on the k-NN graph.
# requires: sq, adata_norm
# outputs: autocorr_df
# params: mode:str, n_perms:int, n_genes:int, alpha:float?

#--- block main
# Spatially variable genes ($mode, n_perms=$n_perms)
autocorr_df = sq.gr.spatial_autocorr(
    adata_norm, mode=$mode, n_perms=$n_perms, corr_method='fdr_bh',
    seed=0, copy=True,
)

#--- block significant
autocorr_df = autocorr_df[autocorr_df['pval_sim_fdr_bh'] < $alpha]

#--- block tail
autocorr_df = autocorr_df.head($n_genes)
```

→ `No problems — safe to save and activate.`

`alpha:float?` is optional for the same reason `n_top_genes:int?` is in
`clustering.leiden`: only one of the two assemblies mentions it, and check #8
would otherwise reject the template. `copy=True` keeps the step from writing
`adata_norm.uns['moranI']` — a template must not mutate state it did not declare.

### Beat 6 — the GUI, generated from the header

`tabs/tab_analyses.py` reads `spec.params` and `spec.assemblies` and renders:

```
┌─ Spatial ▸ Analyses ─────────────────────────────┐
│ Analysis  [ spatial.autocorr            ▾]  ● user│
│ Spatially variable genes by Moran's I / Geary's C │
│                                                   │
│ Variant   [ main+significant+tail       ▾]        │
│ mode      [ moran                        ]        │
│ n_perms   [ 100        ⇅]                         │
│ n_genes   [  10        ⇅]                         │
│ alpha     [ 0.050      ⇅]                         │
│                                                   │
│ [ Check ]  [ Dry run ]            [ Run analysis ]│
└───────────────────────────────────────────────────┘
```

Two things to see here. `Variant` is the assembly combo — block selection is a
*choice* for a user template, not a derived thing (§2b). And `mode` is a bare text
field, because `ParamSpec` has a type but no **choices**: `spec._SYNTH` knows
`str`, not `str in {moran, geary}`. Typing `moran ` with a trailing space passes
every static check and fails inside squidpy. That is a real gap this example
surfaced, and it is now Open Question 5.

### Beat 7 — Run

The tab calls `ctx.ensure_spatial_neighbors(k)` first — idempotent per
`(adata_norm, n_neighs)`, and it calls `ensure_normalized()` itself — then builds:

```python
Step(id="autocorr:moran", **step_template("spatial.autocorr", blocks),
     params={"mode": "moran", "n_perms": 100, "n_genes": 10, "alpha": 0.05},
     deps=["spatial_neighbors"], outputs=["autocorr_df"])
```

Executed and recorded — the same string, by construction:

```python
# Spatially variable genes ('moran', n_perms=100)
autocorr_df = sq.gr.spatial_autocorr(
    adata_norm, mode='moran', n_perms=100, corr_method='fdr_bh',
    seed=0, copy=True,
)
autocorr_df = autocorr_df[autocorr_df['pval_sim_fdr_bh'] < 0.05]
autocorr_df = autocorr_df.head(10)
```

The declared output comes back through `StepExecutor.run`'s contract and is
rendered in a table:

```
              I   pval_sim_fdr_bh
gene22   0.856079        0.013815
gene23   0.850322        0.013815
gene7    0.846926        0.013815
gene19   0.842480        0.013815
```

### Beat 8 — the checks

The seven registry properties, run against the **user** template rather than a
shipped one, over both assemblies:

```
assembly main+tail                parses=OK  self-contained=True  outputs bound=True
assembly main+significant+tail    parses=OK  self-contained=True  outputs bound=True
```

Then the dry-run on `adata_norm[::10]`, which is the only check that would have
caught a draft written against `adata`: it raises where the graph is missing, in
seconds, before the user spends ten minutes on 500k cells with `n_perms=1000`.
The generated `test_spatial_autocorr.py` pins the same two properties plus the
executed result against the `replay_adata` fixture.

### Beat 9 — the notebook

The graph topologically sorts to

```
['preamble', 'normalize', 'spatial_neighbors', 'autocorr:moran']
```

— the new node slotting in behind a dependency it never had to think about
ordering for. And because `template_origin == "user"`, the exported notebook opens
with the banner that already exists, unchanged, for a template id that did not
exist when it was written:

> ## ⚠ This analysis did not use the shipped templates
>
> | step | template | origin | template hash |
> |---|---|---|---|
> | `autocorr:moran` | `spatial.autocorr` | user | `dadf9f98d3a1` |

`scripts/verify_notebook.py` reports `stock_templates: false` for the same reason.
**Nothing had to be added for this to work** — which is the strongest evidence the
design is sitting on the right foundation.

### What this example does *not* solve

Three gaps the exercise exposed, each real:

1. **No param choices.** Beat 6. `mode` should be a combo over `{moran, geary}`,
   and cannot be.
2. **No persistence contract.** `autocorr_df` lives in the executor namespace and
   dies with the session. Every shipped analysis writes its result somewhere —
   `save_clustering_to_adata`, a sidecar via `adata_persistence.sidecar_write_path()`,
   a session zarr array. A user analysis has none of that, so it reruns on every
   launch and its `restore_session` shows nothing. The `outputs` declaration is the
   natural hook (a `pandas` output → a parquet sidecar named for the node id), but
   it does not exist yet.
3. **No terminals.** A user would immediately want `sc.pl.spatial` coloured by the
   top gene. A `plot:*` terminal needs a save path that comes from a dialog, which
   is what `Preview.note` exists for on the preview side — but the generic tab has
   no convention for "this analysis also produces a figure".

None of the three blocks Stage B. §7 answers (2) and (3) properly, because they
turned out to be one question rather than two.

## 7. Result kinds and persistence

The worked example produced a DataFrame, which is the *least* interesting thing a
user analysis can produce. The three that matter in a viewer are a **file**, a
**figure**, and an **obs column that colours the cells** — and the third is the
whole reason to build an analysis inside the viewer instead of in a notebook.

**The plan as written did not account for these.** §2's design carried results back
through `outputs` and stopped there. This section is the correction.

### The key idea: `outputs` already names them; declare their *kind*

`outputs:` says *what comes back*. One more header field says *what it is*, and a
single dispatcher in `tab_analyses.py` does the right thing per kind:

```
# outputs: autocorr_df, svg_score
# produces: autocorr_df = table, svg_score = obs.continuous
```

That is the same shape as `requires` / `outputs` / `assemblies` — a declaration the
call site reads, not logic in a file. Five kinds cover everything the shipped tabs
do:

| kind | example | where the bytes go | comes back via | exists today? |
|---|---|---|---|---|
| `table` | `autocorr_df` | parquet sidecar, `sidecar_write_path(ctx, f"{node_id}.parquet")` | `find_sidecar` at startup | mechanism yes, wiring no |
| `obs.categorical` | a domain / cluster assignment | `adata.obs['clustering_<key>']` via `save_clustering_to_adata` → `_persist_table` → `safe_write_element` | `load_custom_clusterings_from_adata` | **entirely, already** |
| `obs.continuous` | a per-cell score | `adata.obs['analysis_<key>']` + `_persist_table` | a new loader read | storage yes, **colouring no** |
| `figure` | a matplotlib `Figure` | `<data_path>/plots/<stem>.<fmt>` via `ctx.auto_save_plot` | not restored, by design | yes |
| `file` | a CSV export | a path from a save dialog | not restored, by design | yes |

### Persistence is three separate problems

Conflating them is why it looks harder than it is. Each already has a working,
crash-safe implementation; none of them needs inventing.

**1. Where the bytes go.** Anything living on `adata` goes through
`_persist_table` (`utils/adata_persistence.py:162`) → `safe_write_element`, which
is the staged-write-then-two-renames path from `zarr_safe.py`. Anything not
AnnData-shaped goes to `viewer_cache/` via `sidecar_write_path` — *never* the store
root, where files make zarr's hierarchy walk warn on every consolidation and a
rebuild would delete them. Figures go through `ctx.auto_save_plot`. Nothing new.

**2. How it comes back.** A `restore_session` in the tab's exports dict — with the
reminder from §3 that `app.py`'s `all_exports` list is a fourth hardcoded copy of
the tab set, and a tab missing from it restores nothing, silently. The generic tab
needs **one** restore that loops the recorded analyses and rehydrates by kind.

**3. How the app knows it is deletable.** This is the one that will actually bite,
and it is worth being precise, because the obvious guess is wrong in both
directions:

- A cache **rebuild does not lose a user column.** `loader.py:509-521` is a
  *deny-list* — "Anything a freshly-built table does not already have is, by
  definition, something the user's session added." It was an allow-list once and
  silently dropped CNV scores; that is fixed. So an `adata.obs['my_score']`
  survives a rebuild today, with no change.
- But it is **permanently undeletable.** `store_inventory.py:482` reads
  `deletable = column.startswith(loader._USER_OBS_PREFIXES)`, and that tuple is
  `("clustering_", "cluster_labels_", "cnv_score", "copykat_leiden_res")`. The
  documented policy is that **unrecognised defaults to not deletable** — right for
  raw Xenium columns, wrong for a column the user's own analysis just wrote. Tools →
  Dataset would list it, with a size, dimmed and un-tickable, forever.

The fix is one reserved prefix: write user results as `analysis_<slug>` and add
`"analysis_"` to `_USER_OBS_PREFIXES` (and the `uns` equivalent to
`_USER_UNS_KEYS`). Note what does *not* have to change: `assert_deletable` stays
the single choke point, and `tests/test_store_inventory.py` asserts the containment
**property over every node the inventory produces** rather than a remembered list —
so it validates the extension automatically instead of needing to be taught about
it.

### The obs column → colouring answer, in two halves

This splits sharply, and the split is the most useful thing in this section.

**Categorical is nearly free — it already works end to end.** A user analysis that
produces a cluster-like `pd.Series` calls `save_clustering_to_adata(ctx, key,
series)`, sets `ctx.clusterings[key]` and `state["custom_clusterings"][key]`, and
calls `ctx.refresh_clustering_choices()`. From that one sequence it inherits: the
Coloring tab's clustering combo, glasbey categorical colours via
`get_cluster_colors` and `DirectLabelColormap`, the per-cluster filter checkboxes,
the label editor, selection as a `groupby` in **every** analysis tab (nine combos
are refreshed), persistence to the zarr, and rehydration on next launch by
`load_custom_clusterings_from_adata`. None of that needs writing.

Two rules come attached, both enforced:
- The producer **must record a `clustering:<key>` node**, because analysis tabs
  declare `deps=["clustering:<key>"]`. `ctx.record_clustering()` is only a backstop,
  and `tests/test_clustering_recording.py` fails a producer that persists a column
  without recording one. So this kind does **not** get a `analysis:*` node id — the
  id convention is forced by the graph.
- It must go through `save_clustering_to_adata`, which writes the `clustering_<key>`
  prefixed column that `store_inventory._clustering_twin_of` pairs with the bare
  one, so "delete this clustering" removes both instead of leaving an identical
  copy behind.

**Continuous is the real gap, and it is small but not free.**
`CellColorManager.get_continuous_colors` already exists, is already cached, and its
docstring already says it builds colours "for an arbitrary continuous per-cell
score" — but it has **exactly one caller**, `tabs/tab_cnv.py:1179`, which wires its
own UI to reach it. The Coloring tab has exactly two modes, `"Gene Expression"` and
`"Cluster"` (`tab_cell_coloring.py:110-119, 197`). A numeric obs column has no route
to the canvas.

So the work is a **third mode** — "Obs column" — with a combo over numeric
`analysis_*` / `cnv_score*` columns, dispatching to the function that is already
there and already proven by the CNV tab. That is the highest-value item in this
whole design: it is what turns "I wrote an analysis" into "I can see it on the
tissue", and the CNV tab is the existence proof that the hard part is done.

### Node ids follow from the kind

One run may record several nodes — already the shipped pattern, where ROI
expression and its CSV export are two templates and two nodes:

| kind | node id | `kind=` |
|---|---|---|
| `obs.categorical` | `clustering:<key>` (forced, above) | `ARTIFACT` |
| `table`, `obs.continuous` | `analysis:<template_id>:<key>` | `ARTIFACT` |
| `figure` | `plot:<template_id>:<key>` | `TERMINAL` |
| `file` | `export:<template_id>:<key>` | `TERMINAL` |

The terminal kinds matter for a reason beyond tidiness: `notebook_export` can drop
terminals (`graph_to_cells(include_terminals=False)`), and a terminal's cell must
still be **real code** — `sc.pl.*`, `to_csv` — not prose.
`tests/test_recorded_code_is_code.py` parses every recorder call site and fails on
a comment-only node that is not declared `NOTE`. A generated terminal that wrote
`# figure saved to plots/foo.svg` would be caught by an existing test, which is the
right outcome.

### Revised staging

Stage C grows and splits, because persistence is what makes the feature usable
rather than a demo:

| | |
|---|---|
| **C1** | `# produces:` field; the `table` and `obs.categorical` kinds; `analysis_` added to `_USER_OBS_PREFIXES`; one `restore_session` |
| **C2** | Third Coloring mode for `obs.continuous`, over the existing `get_continuous_colors` |
| **C3** | `figure` and `file` terminals via `auto_save_plot` and a save dialog, with `Preview.note` for the dialog path |
| **C4** | Check button, generated `test_<id>.py`, subsample dry-run (the original Stage C) |

## 8. Open questions

- **Where does a user template live?** §2(a) assumes the existing user scope
  (`platformdirs.user_config_dir("palms") / "templates"`). Phase 4b's
  per-dataset scope (`<data_path>/palms_templates/`, *outside*
  `viewer_cache/` whose documented invariant is that everything in it is deletable)
  is more attractive for a new analysis than for an override, because a novel
  analysis is more likely to be dataset-specific and more likely to be worth
  sharing. It also multiplies the upgrade problem: an analysis can arrive on a
  machine running a different app version, where nothing knows it is stale.
- **Does a user analysis get an id namespace?** `custom.*` would make
  `tab_templates._group_of` sort them into their own group instead of "Setup", and
  would make a builtin/user id collision impossible rather than merely resolved.
- **Should the generated `test_<id>.py` be discoverable by the repo's own
  `pytest`?** It must not be — CI asserts what the package *ships*, and
  `conftest.py` blanks the template path precisely so a developer's own files
  cannot change what the suite says.
- **`ParamSpec` needs choices.** Surfaced by §6 beat 6: `mode:str` renders as a
  free text field, and `'moran '` passes every static check before failing inside
  squidpy. A `mode:str{moran,geary}` annotation would give the generated form a
  combo *and* give `validate` a new check worth having. It also improves the
  shipped templates, several of which take a param that is really an enum
  (`method` in `genes.rank_genes`, `flavor` in `clustering.leiden`).
- **What is the unit of "an analysis result" in the session?** §7 answers where
  each kind is *stored*, but not what the session attrs record so a restore knows
  which analyses to rehydrate. `_build_session_attrs` / `_json_safe` already carry a
  JSON-safe dict, and the prov graph already lists every node with its
  `template_id` — so the list may need no new storage at all, only a read. Worth
  settling before C1, because it decides whether `restore_session` walks the graph
  or a separate index.
- **Should `obs.continuous` colouring reuse the cluster-filter machinery?** The
  filter checkboxes are inherently categorical. A continuous score wants a range
  slider, which does not exist. Simplest honest answer for C2: no filter in that
  mode, and say so in the tab.
- **Galaxy's tool XML and tool-versioning docs are still unread.** Flagged in
  `user-configurable-templates-todo.md` as relevant once templates start
  travelling between installations. A user-authored analysis is that case, arriving
  earlier than expected.
