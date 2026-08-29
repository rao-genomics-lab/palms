#!/bin/bash
# Push docs/ content to the GitHub Wiki.
#
# Prerequisites:
#   1. Initialize the wiki by visiting https://github.com/sraorao/palms/wiki
#      and creating any first page (e.g. titled "Home", any content).
#   2. Run this script from the repo root.
#
# Usage:
#   bash scripts/push_to_wiki.sh

set -e

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WIKI_DIR="/tmp/palms.wiki"

echo "Cloning wiki repo..."
rm -rf "$WIKI_DIR"
git clone https://github.com/sraorao/palms.wiki.git "$WIKI_DIR"

echo "Copying docs to wiki..."
# docs/ doubles as the wiki source, so every .md added there would be published
# by default. Internal planning/design notes live there too but are not wiki
# pages. Selecting by *convention* rather than by an exclusion list, because the
# list only ever grew after a note had already leaked: it was extended for
# pyqt6-migration.md and reproducible_notebook_plan.md, and still missed
# user-configurable-templates-todo.md, which reached the public wiki as an
# orphan page.
#
# The convention: a wiki page is Title-Case-With-Hyphens.md (Home, Installation,
# Tab-CNV, Tutorial-Getting-Started) or the _Sidebar.md nav file. Internal notes
# are lower-case, so a new one is repo-only without anyone remembering to say so.
# tests/test_docs_links.py parses this pattern, so it is the single source of
# truth for what is a wiki page.
WIKI_PAGE_RE='^(_Sidebar|[A-Z][A-Za-z0-9-]*)\.md$'

for src in "$REPO_ROOT"/docs/*.md; do
    name="$(basename "$src")"
    if ! [[ "$name" =~ $WIKI_PAGE_RE ]]; then
        echo "  skipping $name (repo-only)"
        continue
    fi
    cp "$src" "$WIKI_DIR"/
done

# Copy screenshots if they exist
if ls "$REPO_ROOT"/docs/screenshots/*.png 2>/dev/null; then
    mkdir -p "$WIKI_DIR"/screenshots
    cp "$REPO_ROOT"/docs/screenshots/*.png "$WIKI_DIR"/screenshots/
fi

echo "Committing and pushing..."
cd "$WIKI_DIR"
git add -A
git commit -m "Update documentation from docs/ in main repo" || echo "Nothing to commit."
git push

echo "Done. Wiki updated: https://github.com/sraorao/palms/wiki"
