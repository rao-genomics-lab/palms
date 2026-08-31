# Publishing the docs on Read the Docs

Repo-only note (lower-case filename ⇒ not a wiki page, and excluded from the built
site by `exclude_docs` in `mkdocs.yml`).

## What is already done

`mkdocs.yml` had shipped a complete 44-page nav for months, and `site_url` already
pointed at `https://palms.readthedocs.io` — but **the site had never been
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

## The site root, and why it is a hook

A GitHub wiki's landing page must be `Home.md`; mkdocs's must be `index.md`. `docs/` is
wiki source, so it has the first and not the second — and **mkdocs builds a site with no
page at its root without complaint, `--strict` included**. RTD does not serve one: the
first build failed with *"Index file is not present in HTML output directory"* after this
repo's CI had been building the site green on every push.

Renaming is not available in either direction. `Home.md` is what `push_to_wiki.sh`
publishes as the wiki home, and a lower-case `docs/index.md` would be excluded from the
site by `exclude_docs` (`[a-z]*.md`) — silently, because that pattern is what keeps
internal notes like this file off the published site. So `mkdocs_hooks.on_files` retargets
`Home.md` to the root at build time, and CI checks `index.html` exists afterwards, since
`--strict` cannot.

## Connecting the project

The repository is public (2026-08-31), which was the one blocker — RTD Community does not
build private repos, and that is why this was never worth paying RTD Business to work
around.

1. Sign in at <https://readthedocs.org> with GitHub and grant access to the repo. A repo
   that has been *renamed* usually needs **Settings → Connected Services → GitHub →
   Resync** before it appears in the import list.
2. **Import a Project** → `rao-genomics-lab/palms`. RTD detects `.readthedocs.yaml`;
   no settings need changing. The slug must come out as **`palms`**, since `site_url`
   already says so — if an older `xenium-viewer` project is still in the dashboard, delete
   it rather than renaming: on RTD Community the slug is fixed at creation, and that
   project has never built.
3. Confirm the first build is green, then check a cross-link on the rendered site —
   that is the thing that was broken and is the only part CI cannot prove.
4. Add the RTD webhook (RTD offers it during import) so pushes rebuild.
5. Leave `site_url` as it is. If the project ends up under a different slug, update
   `site_url` to match rather than leaving it pointing at a site that does not exist.

The GitHub Wiki stays as-is; `scripts/push_to_wiki.sh` is unaffected by any of this,
and `docs/` remains the single source for both.
