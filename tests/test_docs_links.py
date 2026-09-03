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
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DOCS = REPO / "docs"

# Imported for its anchor rule only, so the check here and the anchors the
# generator emits cannot disagree. Stdlib at import time; the package is only
# imported inside its render functions, so this stays dependency-free.
sys.path.insert(0, str(REPO / "scripts"))
import generate_docs  # noqa: E402
SCREENSHOTS = DOCS / "screenshots"
PUBLISH_SCRIPT = REPO / "scripts" / "push_to_wiki.sh"
MKDOCS = REPO / "mkdocs.yml"
SIDEBAR = DOCS / "_Sidebar.md"
APP = REPO / "src" / "palms" / "app.py"

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


def _anchors_in(page: Path) -> set[str]:
    """Every same-page link target a reader could land on, by GitHub's rules."""
    anchors = set()
    fenced = False
    for line in page.read_text().splitlines():
        if line.lstrip().startswith("```"):
            fenced = not fenced
            continue
        if not fenced and line.startswith("#"):
            anchors.add(generate_docs._anchor(line.lstrip("#").strip()))
    return anchors


def test_link_anchors_resolve():
    """A ``Page#anchor`` link must land on a heading that exists.

    ``test_internal_links_resolve`` splits the anchor off and checks only the
    page, so a link to a heading that was renamed — or never existed — passed
    while silently dropping the reader at the top of the page. Cross-references
    into the generated pages are written by hand and use anchors GitHub derives
    by a rule that is easy to get wrong (underscores are kept, dots are not), so
    they need checking.
    """
    broken = []
    for md in _md_files():
        for raw in LINK_RE.findall(md.read_text()):
            target = raw.split(" ", 1)[0].strip()
            if "://" in target or "#" not in target or target.startswith("mailto:"):
                continue
            page_part, _, anchor = target.partition("#")
            if not anchor:
                continue
            page = md if not page_part else DOCS / f"{page_part}.md"
            if not page.exists():
                continue                      # test_internal_links_resolve owns this
            if f"#{anchor}" not in _anchors_in(page):
                broken.append(f"{md.name} -> {target}")
    assert not broken, "links to non-existent anchors:\n  " + "\n  ".join(broken)


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
    "Publish": "Tab-Publish.md",
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


def _addtab_order() -> dict[int, list[str]]:
    """Tab labels per group, in ``addTab()`` order, keyed by group index.

    The inner tabs are added to ``<group>_tabs`` before the groups are added to
    the outer ``tab_widget``, so the receiver name identifies the group and the
    outer order comes from the ``tab_widget.addTab`` calls.
    """
    tree = ast.parse(APP.read_text())
    per_receiver: dict[str, list[str]] = {}
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "addTab"
                and isinstance(node.func.value, ast.Name)):
            continue
        label = next((a.value for a in node.args
                      if isinstance(a, ast.Constant) and isinstance(a.value, str)), None)
        widget = next((a.id for a in node.args if isinstance(a, ast.Name)), None)
        if label is None:
            continue
        per_receiver.setdefault(node.func.value.id, []).append((label, widget))

    outer = per_receiver.get("tab_widget", [])
    order = {}
    for index, (_group_label, group_widget) in enumerate(outer):
        order[index] = [lbl for lbl, _ in per_receiver.get(group_widget, [])]
    return order


def test_screenshot_indices_point_at_the_tab_they_are_named_for():
    """``TAB_SHOTS`` positions must match ``app.py``'s ``addTab`` order.

    The indices in ``scripts/capture_screenshots.py`` are positional, so
    inserting a tab shifts every entry below it — and the failure is silent: the
    run succeeds and writes a picture of the wrong tab under the right name.
    That happened once already (``tab-notebook.png`` was a picture of Crop
    Dataset for as long as the list was one short), and the script carries a
    comment about it. This turns that comment into a check.
    """
    capture = (REPO / "scripts" / "capture_screenshots.py").read_text()
    tree = ast.parse(capture)
    shots = None
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign)
                and any(getattr(t, "id", None) == "TAB_SHOTS" for t in node.targets)):
            shots = ast.literal_eval(node.value)
    assert shots, "TAB_SHOTS not found in capture_screenshots.py"

    order = _addtab_order()
    # A screenshot file name is the lower-cased page name: Tab-Cell-Coloring.md
    # is illustrated by tab-cell-coloring.png.
    label_for_stem = {page[: -len(".md")].lower(): label
                      for label, page in TAB_PAGES.items()}

    wrong, checked = [], 0
    for outer, inner, fname in shots:
        stem = "tab-" + fname[len("tab-"):-len(".png")]
        expected_label = label_for_stem.get(stem)
        if expected_label is None:
            wrong.append(f"{fname}: no page in TAB_PAGES is named for it")
            continue
        labels = order.get(outer, [])
        actual = labels[inner] if inner < len(labels) else "<out of range>"
        checked += 1
        if actual != expected_label:
            wrong.append(f"{fname}: ({outer},{inner}) is {actual!r}, "
                         f"expected {expected_label!r}")
    assert not wrong, "capture_screenshots.py indices are stale:\n  " + "\n  ".join(wrong)
    # Without this the whole test passes by matching nothing, which is how the
    # first version of it "passed" against 26 entries it never looked at.
    assert checked == len(shots), (
        f"only {checked} of {len(shots)} TAB_SHOTS entries were checked"
    )


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


# ── mkdocs: the wiki link style, rendered ────────────────────────────────────
#
# docs/ is authored for the GitHub Wiki, where links are extensionless. mkdocs
# reads those as paths to files that do not exist, so the built site had 158 dead
# cross-links — the reason `mkdocs.yml` could ship a complete nav for months while
# nothing on the rendered site actually worked. `mkdocs_hooks.py` rewrites them at
# build time so the source keeps one convention.
#
# These tests exercise the hook directly rather than running mkdocs: the rewrite
# rule is what can be wrong, and CI builds the site for real in the lint job.

sys.path.insert(0, str(REPO))
import mkdocs_hooks  # noqa: E402


def test_bare_wiki_links_gain_an_extension():
    out = mkdocs_hooks.on_page_markdown("see [Clustering](Tab-Clustering) for more")
    assert "(Tab-Clustering.md)" in out


def test_an_anchor_survives_the_rewrite():
    """`Analysis-Templates#normalize` is a real link in API-Reference.md."""
    out = mkdocs_hooks.on_page_markdown("[normalize](Analysis-Templates#normalize)")
    assert "(Analysis-Templates.md#normalize)" in out


@pytest.mark.parametrize("target", [
    "https://example.com/x",
    "http://example.com/x",
    "mailto:someone@example.com",
    "#a-heading-on-this-page",
    "screenshots/tab-cnv.png",
    "Tab-CNV.md",
    "./relative.md",
    "/absolute",
])
def test_links_that_already_mean_what_they_say_are_untouched(target):
    src = f"[x]({target})"
    assert mkdocs_hooks.on_page_markdown(src) == src


def test_images_are_left_alone():
    """An image target carries an extension, and must never gain `.md`."""
    src = "![Crop Dataset](screenshots/tab-crop-dataset.png)"
    assert mkdocs_hooks.on_page_markdown(src) == src


def test_every_bare_link_in_docs_rewrites_onto_a_real_page():
    """The rewrite has to land on a file, or the site 404s as before.

    Derived from the pages actually present, not from a remembered list — the
    same rule the rest of this module follows.
    """
    pages = {p.name for p in DOCS.glob("*.md")}
    missing = set()
    for md in sorted(DOCS.glob("*.md")):
        for target in LINK_RE.findall(md.read_text()):
            rewritten = mkdocs_hooks._rewrite(target)
            if rewritten == target:
                continue                       # untouched: external, anchor, has an extension
            page = rewritten.split("#", 1)[0]
            if page not in pages:
                missing.add(f"{md.name} -> {target}")
    assert not missing, f"links that rewrite onto a page that does not exist: {sorted(missing)}"


# ── The site needs a page at its root ────────────────────────────────────────
# A GitHub wiki's landing page is `Home.md`; mkdocs's is `index.md`. `docs/` is
# wiki source, so it has the first and not the second, and mkdocs builds a site
# with nothing at its root without complaint — `--strict` passed here and in CI
# while Read the Docs refused the result ("Index file is not present in HTML
# output directory"). `mkdocs_hooks.on_files` maps one to the other.
#
# Exercised through a stand-in for mkdocs's `File`, since mkdocs is not in the
# test env (it is installed in the lint job, which builds the site for real and
# then checks that `index.html` exists).

class _FakeFile:
    def __init__(self, src_uri):
        self.src_uri = src_uri
        self.name = src_uri.rsplit(".", 1)[0]
        self.dest_uri = f"{self.name}/index.html"
        self.url = f"{self.name}/"
        self.dest_dir = "/site"
        self.abs_dest_path = f"/site/{self.dest_uri}"


def test_the_wiki_home_page_becomes_the_site_root():
    home = _FakeFile("Home.md")
    mkdocs_hooks.on_files([home], config=None)

    assert home.dest_uri == "index.html"
    assert home.abs_dest_path == "/site/index.html"
    assert home.url == "./"


def test_no_other_page_is_moved():
    page = _FakeFile("Tab-Clustering.md")
    mkdocs_hooks.on_files([page], config=None)

    assert page.dest_uri == "Tab-Clustering/index.html"


def test_the_home_page_the_hook_looks_for_exists():
    """The hook matches by name, so a rename in `docs/` must not pass silently:
    it would leave the site with no root page and fail only at Read the Docs."""
    assert (DOCS / "Home.md").exists()


def test_an_index_md_would_be_excluded_from_the_site():
    """Why the mapping is a hook and not a file. `exclude_docs` treats lower-case
    as "internal note", which is how planning notes are kept off the site — so
    adding `docs/index.md` would drop it from the build without saying so."""
    assert not (DOCS / "index.md").exists()
    assert "[a-z]*.md" in MKDOCS.read_text()
