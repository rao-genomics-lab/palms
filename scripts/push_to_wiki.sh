#!/bin/bash
# Push docs/ content to the GitHub Wiki.
#
# Prerequisites:
#   1. Initialize the wiki by visiting https://github.com/sraorao/xenium_viewer/wiki
#      and creating any first page (e.g. titled "Home", any content).
#   2. Run this script from the repo root.
#
# Usage:
#   bash scripts/push_to_wiki.sh

set -e

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WIKI_DIR="/tmp/xenium_viewer.wiki"

echo "Cloning wiki repo..."
rm -rf "$WIKI_DIR"
git clone https://github.com/sraorao/xenium_viewer.wiki.git "$WIKI_DIR"

echo "Copying docs to wiki..."
# docs/ doubles as the wiki source, so every .md added there is published by
# default. Internal planning/design notes live there too but are not wiki pages;
# list them here to keep them repo-only.
WIKI_EXCLUDE=(
    pyqt6-migration.md
    reproducible_notebook_plan.md
)

for src in "$REPO_ROOT"/docs/*.md; do
    name="$(basename "$src")"
    skip=""
    for excluded in "${WIKI_EXCLUDE[@]}"; do
        [ "$name" = "$excluded" ] && skip=1 && break
    done
    if [ -n "$skip" ]; then
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

echo "Done. Wiki updated: https://github.com/sraorao/xenium_viewer/wiki"
