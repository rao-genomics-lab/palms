"""Keep ``docs/`` consistent with itself, with mkdocs, and with the GUI.

``docs/`` is three things at once: the GitHub Wiki source (published by
``scripts/push_to_wiki.sh``), the mkdocs site source (``mkdocs.yml``), and the
place a reader is sent from the app. Nothing kept those three in step, and they
drifted: three shipped tabs had no page and appeared in no navigation,
``mkdocs.yml`` omitted two pages that had existed for weeks, an internal
planning note was published to the public wiki as an orphan, and three links
pointed at a page that has never existed.

Every check here is derived from a source of truth rather than a remembered
list — the wiki-page pattern is parsed out of the publish script, the tab labels
out of ``app.py`` — so the checks cannot themselves go stale in the way the docs
did. Pure stdlib: no Qt, no dataset, no network.
"""
from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DOCS = REPO / "docs"
SCREENSHOTS = DOCS / "screenshots"
PUBLISH_SCRIPT = REPO / "scripts" / "push_to_wiki.sh"
MKDOCS = REPO / "mkdocs.yml"
SIDEBAR = DOCS / "_Sidebar.md"
APP = REPO / "src" / "xenium_viewer" / "app.py"

# ``![alt](target)`` and ``[label](target)``. The negative lookbehind keeps the
# link pattern from also matching the image one.
IMAGE_RE = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")
LINK_RE = re.compile(r"(?<!!)\[[^\]]*\]\(([^)]+)\)")

# Pages that are wiki nav constructs rather than content, so they are not
# expected to appear in the mkdocs nav or to be linked from the sidebar.
NAV_ONLY = {"_Sidebar.md"}


def _md_files() -> list[Path]:
    return sorted(DOCS.glob("*.md"))


def _wiki_page_pattern() -> re.Pattern[str]:
    """The publish script's own definition of what counts as a wiki page.

    Parsed rather than duplicated: if the convention changes there, these tests
    follow it instead of contradicting it.
    """
    text = PUBLISH_SCRIPT.read_text()
    match = re.search(r"^WIKI_PAGE_RE='([^']+)'", text, re.MULTILINE)
    assert match, "push_to_wiki.sh no longer defines WIKI_PAGE_RE — update this test"
    return re.compile(match.group(1))


def _published_pages() -> set[str]:
    """The set of docs/*.md that push_to_wiki.sh would copy to the wiki."""
    pattern = _wiki_page_pattern()
    return {p.name for p in _md_files() if pattern.match(p.name)}


def _mkdocs_nav_pages() -> set[str]:
    """Every ``*.md`` referenced from mkdocs.yml's nav.

    Read as text on purpose — pyyaml is not a test dependency, and the nav is a
    flat set of ``Title: Page.md`` leaves, so a regex is honest here.
    """
    text = MKDOCS.read_text()
    nav = text.split("nav:", 1)[1] if "nav:" in text else ""
    return set(re.findall(r"([A-Za-z0-9_-]+\.md)\s*$", nav, re.MULTILINE))


def _sidebar_targets() -> set[str]:
    return {f"{t}.md" for t in _link_targets(SIDEBAR)}


def _link_targets(path: Path) -> list[str]:
    """Internal wiki-style link targets in one file.

    Skips external URLs, anchors, and explicit ``.md`` links; what remains is
    the extensionless form the GitHub Wiki resolves.
    """
    targets = []
    for raw in LINK_RE.findall(path.read_text()):
        target = raw.split(" ", 1)[0].strip()
        if not target or "://" in target or target.startswith(("#", "mailto:")):
            continue
        if target.endswith(".md") or "/" in target:
            continue
        targets.append(target.split("#", 1)[0])
    return targets


def test_internal_links_resolve():
    """Every wiki-style link points at a page that exists.

    The case that motivated this: ``[Gene Analysis](Tab-Gene-Analysis)`` in two
    files, for a page that has never existed under that name.
    """
    broken = []
    for md in _md_files():
        for target in _link_targets(md):
            if not (DOCS / f"{target}.md").exists():
                broken.append(f"{md.name} -> {target}")
    assert not broken, "links to non-existent pages:\n  " + "\n  ".join(broken)


def test_referenced_screenshots_exist():
    missing = []
    for md in _md_files():
        for target in IMAGE_RE.findall(md.read_text()):
            target = target.split(" ", 1)[0].strip()
            if "://" in target:
                continue
            if not (DOCS / target).exists():
                missing.append(f"{md.name} -> {target}")
    assert not missing, "referenced images that do not exist:\n  " + "\n  ".join(missing)


def test_no_internal_notes_are_published():
    """Planning notes stay in the repo.

    ``user-configurable-templates-todo.md`` reached the public wiki because the
    old exclusion list had to be extended by hand for each new note. Assert the
    convention holds: a lower-case filename is repo-only.
    """
    published = _published_pages()
    leaked = sorted(n for n in published if n[0].islower())
    assert not leaked, f"internal notes would be published to the wiki: {leaked}"


def test_every_published_page_is_in_the_mkdocs_nav():
    """The check that would have caught every navigation gap in one go."""
    expected = _published_pages() - NAV_ONLY
    missing = sorted(expected - _mkdocs_nav_pages())
    assert not missing, f"pages absent from the mkdocs.yml nav: {missing}"


def test_every_published_page_is_linked_from_the_sidebar():
    expected = _published_pages() - NAV_ONLY
    missing = sorted(expected - _sidebar_targets())
    assert not missing, f"pages not linked from docs/_Sidebar.md: {missing}"


def test_mkdocs_nav_has_no_dangling_entries():
    dangling = sorted(p for p in _mkdocs_nav_pages() if not (DOCS / p).exists())
    assert not dangling, f"mkdocs.yml nav points at missing files: {dangling}"


# ── GUI tabs ↔ reference pages ───────────────────────────────────────────────

# The GUI label is abbreviated to fit the tab bar; the page keeps the full name.
# This mapping is the reconciliation between the two, and the reason a new tab
# with no documentation turns CI red. Keep it in the app's own tab order.
TAB_PAGES = {
    "Clustering": "Tab-Clustering.md",
    "Coloring": "Tab-Cell-Coloring.md",
    "Transcripts": "Tab-Transcripts.md",
    "UMAP": "Tab-UMAP.md",
    "Rank Genes": "Tab-Rank-Genes.md",
    "Markers": "Tab-Markers.md",
    "Correlation": "Tab-Gene-Correlation.md",
    "CNV": "Tab-CNV.md",
    "ROI DEG": "Tab-ROI-Analysis.md",
    "Lig-Rec": "Tab-Ligand-Receptor.md",
    "Nhood Enrich": "Tab-Neighborhood-Enrichment.md",
    "Co-occur": "Tab-Co-occurrence.md",
    "Domains": "Tab-Domains.md",
    "Annot Nhood": "Tab-Annot-Nhood.md",
    "Annot Dist": "Tab-Annot-Distance.md",
    "H&E": "Tab-HE-Registration.md",
    "ARMS": "Tab-ARMS-Overlay.md",
    "Ext Images": "Tab-External-Images.md",
    "Patches": "Tab-Patches.md",
    "Annotations": "Tab-Annotations.md",
    "Segmentation": "Tab-Segmentation.md",
    "Crop Dataset": "Tab-Crop-Dataset.md",
    "Notebook": "Tab-Notebook.md",
    "Dataset": "Tab-Dataset.md",
    "Cache": "Tab-Cache.md",
    "Templates": "Tab-Templates.md",
}

# Outer group labels, which are containers rather than documented tabs.
GROUP_LABELS = {"Cells", "Genes", "Spatial", "Images", "Tools"}


def _gui_tab_labels() -> set[str]:
    """``addTab(widget, "Label")`` literals from app.py, read with ``ast``.

    Parsed rather than imported: importing app.py would pull in Qt, napari and
    spatialdata for what is a question about a string literal.
    """
    tree = ast.parse(APP.read_text())
    labels = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == "addTab"):
            continue
        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                labels.add(arg.value)
    return labels


def test_every_gui_tab_has_a_reference_page():
    labels = _gui_tab_labels() - GROUP_LABELS
    assert labels, "found no addTab() labels in app.py — has the parsing broken?"

    undocumented = sorted(lbl for lbl in labels if lbl not in TAB_PAGES)
    assert not undocumented, (
        f"GUI tabs with no entry in TAB_PAGES: {undocumented}. "
        "Write the page, then add it here."
    )

    missing = sorted(
        f"{lbl} -> {page}" for lbl, page in TAB_PAGES.items()
        if lbl in labels and not (DOCS / page).exists()
    )
    assert not missing, f"reference pages that do not exist: {missing}"


def test_tab_pages_mapping_has_no_stale_entries():
    """A removed tab should not leave a mapping behind claiming it exists."""
    labels = _gui_tab_labels() - GROUP_LABELS
    stale = sorted(lbl for lbl in TAB_PAGES if lbl not in labels)
    assert not stale, f"TAB_PAGES names tabs that are no longer in app.py: {stale}"


# ── House style ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("page", sorted(DOCS.glob("Tab-*.md")), ids=lambda p: p.name)
def test_tab_pages_follow_the_house_style(page: Path):
    """Every reference page is ``# Title`` then Controls / Workflow / Notes.

    All 23 pages that predate this test already comply; pinning it keeps the
    next one consistent without depending on a reviewer noticing.
    """
    lines = page.read_text().splitlines()
    assert lines and lines[0].startswith("# "), (
        f"{page.name} must open with a level-1 heading"
    )
    headings = [ln for ln in lines if ln.startswith("## ")]
    assert headings == ["## Controls", "## Workflow", "## Notes"], (
        f"{page.name} has sections {headings}, expected "
        "['## Controls', '## Workflow', '## Notes']"
    )
