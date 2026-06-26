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
# Copy all markdown files (but not mkdocs.yml)
cp "$REPO_ROOT"/docs/*.md "$WIKI_DIR"/

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
