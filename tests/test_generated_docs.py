"""The two generated reference pages must describe the code that exists.

``docs/Analysis-Templates.md`` and ``docs/API-Reference.md`` carry *code* — every
template's shipped body, and the signature of every function a notebook is
invited to call. Hand-maintained, both go stale the moment a template is edited
or an argument renamed, and a documented call that no longer exists is worse than
no documentation at all: the reader trusts it.

So they are generated from the registry and from the live objects, and these
tests are what make that claim true rather than aspirational. Each one is aimed
at a specific way this rots — a page not regenerated after a template change, a
new template shipped without prose, a function renamed out from under an entry.

Pure stdlib plus the package itself: no Qt, no dataset, no network.
"""
from __future__ import annotations

import importlib
import inspect
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DOCS = REPO / "docs"
SCRIPTS = REPO / "scripts"

if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

generate_docs = pytest.importorskip("generate_docs")


# ── Freshness ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("page", sorted(generate_docs.PAGES), ids=lambda p: p.name)
def test_page_is_not_stale(page: Path):
    """The checked-in page equals what the generator produces right now."""
    assert page.exists(), f"{page.name} missing — run python scripts/generate_docs.py"
    produced = generate_docs.PAGES[page]()
    assert page.read_text() == produced, (
        f"{page.name} is out of date. Regenerate with:\n"
        f"    python scripts/generate_docs.py"
    )


# ── Templates ────────────────────────────────────────────────────────────────

def _builtin_ids() -> list[str]:
    from palms.utils.step_templates import loader as registry
    return registry.builtin_ids()


def test_every_template_has_prose():
    """A new template cannot ship undocumented.

    The generator raises rather than emitting a section with nothing in it, so
    this asserts the same property at the point where it is cheap to read.
    """
    missing = sorted(set(_builtin_ids()) - set(generate_docs.TEMPLATE_NOTES))
    assert not missing, (
        f"template(s) {missing} have no TemplateNote in scripts/generate_docs.py"
    )


def test_no_prose_for_templates_that_no_longer_exist():
    """The other direction: a deleted template must not keep its section."""
    extra = sorted(set(generate_docs.TEMPLATE_NOTES) - set(_builtin_ids()))
    assert not extra, (
        f"TEMPLATE_NOTES documents {extra}, which the registry does not ship"
    )


def test_every_template_group_is_known():
    """A typo'd group would silently drop a template out of the page."""
    unknown = {tid: note.group
               for tid, note in generate_docs.TEMPLATE_NOTES.items()
               if note.group not in generate_docs.GROUP_ORDER}
    assert not unknown, f"unknown group(s): {unknown}"


def test_every_template_body_appears_verbatim():
    """The published code is the shipped code, block for block.

    Not a re-run of the generator — this reads the page as a reader would and
    looks for the actual template text in it, so a formatting change that
    mangles the source (a stray indent, a swallowed line) fails here even though
    the freshness test would still pass.
    """
    from palms.utils.step_templates import loader as registry

    page = (DOCS / "Analysis-Templates.md").read_text()
    for tid in _builtin_ids():
        spec = registry.builtin_spec(tid)
        for name, block in spec.blocks.items():
            body = block.text.strip("\n")
            assert body in page, (
                f"block {name!r} of template {tid!r} does not appear verbatim in "
                f"Analysis-Templates.md"
            )


# ── API allow-list ───────────────────────────────────────────────────────────

def _api_entries():
    for section in generate_docs.API_SECTIONS:
        for entry in section.entries:
            yield section.title, entry


@pytest.mark.parametrize(
    "section_title,entry",
    list(_api_entries()),
    ids=lambda v: v.name if hasattr(v, "name") else str(v),
)
def test_documented_object_exists(section_title, entry):
    """Every documented callable imports and is callable.

    This is the check that earns the page its keep: a rename or a removal turns
    CI red instead of leaving behind a plausible-looking call that raises
    ``AttributeError`` in a reader's notebook.
    """
    module = importlib.import_module(entry.module)
    assert hasattr(module, entry.name), (
        f"{entry.module}.{entry.name} is documented under {section_title!r} but "
        f"does not exist"
    )
    obj = getattr(module, entry.name)
    assert callable(obj), f"{entry.module}.{entry.name} is not callable"


def test_documented_signatures_match_the_page():
    """The signature printed on the page is the live one."""
    page = (DOCS / "API-Reference.md").read_text()
    for _, entry in _api_entries():
        sig, _ = generate_docs._signature(entry.module, entry.name)
        assert f"{entry.name}{sig}" in page, (
            f"{entry.name}'s signature on the page is not the live one"
        )


def test_api_entries_are_not_duplicated():
    seen = [(e.module, e.name) for _, e in _api_entries()]
    dupes = {x for x in seen if seen.count(x) > 1}
    assert not dupes, f"documented twice: {sorted(dupes)}"


# ── Nothing machine-specific ─────────────────────────────────────────────────

# Absolute paths reach a generated page through default arguments — several
# resolve against the working directory at import time, so a raw repr bakes in
# whatever machine ran the generator. That is wrong for every reader, and it
# already happened once (``TranscriptLoader``'s cache_dir/parquet_path).
_MACHINE_PATH = re.compile(r"(/home/[a-z]|/Users/|PosixPath\('/)")


@pytest.mark.parametrize("page", sorted(generate_docs.PAGES), ids=lambda p: p.name)
def test_no_machine_specific_paths(page: Path):
    hits = [ln for ln in page.read_text().splitlines() if _MACHINE_PATH.search(ln)]
    assert not hits, f"{page.name} contains machine-specific path(s): {hits[:3]}"


def test_generated_pages_declare_that_they_are_generated():
    """A reader who edits one by hand should be told not to, at the top."""
    for page in generate_docs.PAGES:
        head = page.read_text()[:400]
        assert "generate_docs.py" in head, (
            f"{page.name} does not say it is generated"
        )


# ── The one place the override rule inverts ──────────────────────────────────

def test_the_catalogue_reads_builtin_text_not_resolved_text():
    """The catalogue documents the *shipped* defaults, deliberately.

    Everywhere else, reading a template with ``builtin_*`` instead of
    ``step_template`` is the bug that silently disabled customisation for
    ``genes.marker_plot``. Here it is correct — a published page must not vary
    with the config of whoever ran the generator — so this pins the intent
    rather than leaving it to be "fixed" later.
    """
    source = (SCRIPTS / "generate_docs.py").read_text()
    render = source.split("def render_templates_page")[1].split("\ndef ")[0]
    assert "builtin_spec" in render and "resolve(" not in render, (
        "render_templates_page must read the shipped templates, not the "
        "user-resolved ones"
    )


def test_generator_is_importable_without_side_effects():
    """Importing it must not write anything — the test suite imports it."""
    source = inspect.getsource(generate_docs)
    body = source.split('if __name__ == "__main__":')[0]
    assert "write_text(" in body, "sanity: the generator does write, inside main()"
    for line in body.splitlines():
        stripped = line.strip()
        if stripped.startswith("#") or not stripped:
            continue
        if "write_text(" in stripped:
            assert line.startswith(" "), (
                f"module-level write in generate_docs.py: {stripped}"
            )
