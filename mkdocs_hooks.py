"""Make the GitHub-Wiki link style resolve under mkdocs.

``docs/`` is authored as wiki source — that is what ``scripts/push_to_wiki.sh``
publishes and what ``tests/test_docs_links.py`` treats as the single source of
truth for what a page is. Wiki links are extensionless: ``[Clustering](Tab-Clustering)``.

GitHub's wiki resolves those. mkdocs does not — it reads them as relative paths to
files that do not exist, so **every cross-link on the built site 404s**, which is
why the site was never buildable despite ``mkdocs.yml`` having shipped a complete
nav for months.

Rewriting the ~40 links in ``docs/`` would fix mkdocs and introduce a second link
convention for authors to remember, one the wiki publisher does not check. So the
rewrite happens here, at build time, and the source keeps one convention.

This file lives at the repo root rather than in ``docs/``: the publish script
selects wiki pages by filename convention, and a stray module inside ``docs/``
is exactly the kind of thing that has leaked to the public wiki before.
"""

from __future__ import annotations

import os
import posixpath
import re

#: ``[label](target)`` but not ``![alt](target)`` — an image target already
#: carries an extension, so it fails the test below anyway; the lookbehind just
#: keeps the intent legible.
_LINK = re.compile(r"(?<!!)(\[[^\]]*\])\(([^)\s]+)\)")

#: Anything with a scheme, an anchor-only target, or an explicit path is left
#: exactly as written — those already mean what they say to mkdocs.
_LEAVE_ALONE = ("http://", "https://", "mailto:", "ftp://", "//", "#", "/", ".")


def _rewrite(target: str) -> str:
    if target.startswith(_LEAVE_ALONE):
        return target

    path, sep, anchor = target.partition("#")
    if not path:                       # pure "#anchor"
        return target
    if posixpath.splitext(path)[1]:    # already has an extension (.md, .png, …)
        return target

    return f"{path}.md{sep}{anchor}"


def on_page_markdown(markdown: str, **_kwargs) -> str:
    """Append ``.md`` to bare wiki-style link targets."""
    return _LINK.sub(lambda m: f"{m.group(1)}({_rewrite(m.group(2))})", markdown)


def on_files(files, config, **_kwargs):
    """Serve the wiki's ``Home.md`` at the site root, as ``index.html``.

    A GitHub wiki's landing page must be called ``Home.md``; mkdocs's is
    ``index.md``. ``docs/`` is wiki source, so it has the first and not the
    second — and mkdocs is happy to build a site with no page at its root. Read
    the Docs is not: the first build after the repo went public failed with
    "Index file is not present in HTML output directory", after ``mkdocs build
    --strict`` had passed locally and in CI. Nothing in mkdocs itself checks it.

    Renaming the file is not an option in either direction: ``Home.md`` is what
    ``scripts/push_to_wiki.sh`` publishes as the wiki home, and a lower-case
    ``index.md`` would be dropped by ``exclude_docs`` (``[a-z]*.md``) — silently,
    since that pattern is how internal notes are kept out of the site.

    So the mapping happens here, with the same reasoning as the link rewrite
    above: the source keeps the wiki's one convention. Only ``File.name`` is
    set; ``dest_uri`` and ``url`` are lazy and mkdocs derives both from it —
    the same route it already takes for ``README.md``. Setting them by hand
    would hard-code the two places mkdocs decides directory-URL layout.
    """
    for file in files:
        if file.src_uri == "Home.md":
            file.name = "index"
            file.dest_uri = "index.html"
            file.url = "./"
            file.abs_dest_path = os.path.join(file.dest_dir, "index.html")
    return files
