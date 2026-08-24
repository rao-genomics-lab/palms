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
