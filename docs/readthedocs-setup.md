# Publishing the docs on Read the Docs

Repo-only note (lower-case filename ⇒ not a wiki page, and excluded from the built
site by `exclude_docs` in `mkdocs.yml`).

## What is already done

`mkdocs.yml` had shipped a complete 44-page nav for months, and `site_url` already
pointed at `https://xenium-viewer.readthedocs.io` — but **the site had never been
built**, and would not have worked if it had. `docs/` is authored as GitHub Wiki
source, where links are extensionless (`[Clustering](Tab-Clustering)`); mkdocs reads
those as paths to files that do not exist. Measured: **158 dead cross-links**.

In place now:

- `mkdocs_hooks.py` — appends `.md` to bare wiki link targets at build time, so the
  source keeps one link convention instead of two.
- `mkdocs.yml` — `hooks:`, `exclude_docs:` (internal notes and `_Sidebar.md` stay out
  of the published site), and `validation:` turning dead links into warnings so
  `--strict` fails on them.
- `.readthedocs.yaml` + `requirements-docs.txt` — Ubuntu 24.04, Python 3.12, mkdocs,
  `fail_on_warning: true`.
- CI builds the site with `--strict` on every push, in the **lint** job.

The build needs no conda env, no scanpy and no Qt: `Analysis-Templates.md` and
`API-Reference.md` are generated and checked in (`scripts/generate_docs.py`, guarded
by `tests/test_generated_docs.py`). Keep it that way — a docs build that needs the
application is a docs build that breaks whenever the application does.

Build it locally with:

```bash
pip install -r requirements-docs.txt
mkdocs build --strict     # or: mkdocs serve
```

## What is left, and what blocks it

**Connecting the Read the Docs project requires the repository to be public.**
RTD Community only builds public repos; `sraorao/xenium_viewer` is still private.
This is the same step the preprint plan's Phase A needs, so the sequence is: make the
repo public, then connect RTD — not pay for RTD Business to work around a step that
is happening anyway.

Once public:

1. Sign in at <https://readthedocs.org> with GitHub and grant access to the repo.
2. **Import a Project** → `sraorao/xenium_viewer`. RTD detects `.readthedocs.yaml`;
   no settings need changing.
3. Confirm the first build is green, then check a cross-link on the rendered site —
   that is the thing that was broken and is the only part CI cannot prove.
4. Add the RTD webhook (RTD offers it during import) so pushes rebuild.
5. Leave `site_url` as it is. If the project ends up under a different slug, update
   `site_url` to match rather than leaving it pointing at a site that does not exist.

The GitHub Wiki stays as-is; `scripts/push_to_wiki.sh` is unaffected by any of this,
and `docs/` remains the single source for both.
