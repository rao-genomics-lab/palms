"""What happens to a customised template when the shipped one moves on.

This is the ``dpkg`` conffile problem, and most of it is already solved by the
shape of the storage rather than by logic: an override records **only** the
blocks the user changed, so every other block resolves against whatever the
current release ships. The "unmodified file, replace silently" case that
``dpkg`` has to detect simply cannot arise here.

What is left is the genuinely hard case — a block the user *did* change, whose
shipped version has since changed too. Their edit still applies (silently
dropping it would be far worse), but they are told, because the upstream change
they never saw might be a correctness fix their edit is now shadowing. Telling
that apart from an ordinary customisation is the whole job of
``overrides.json``: it records the hash of the **shipped** text each block was
forked from, so a later release is compared against what the user diverged
*from* rather than against their own edit.

These tests simulate a release by patching the builtin registry, which is the
only way to exercise this without waiting for an actual version bump.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from xenium_viewer.utils.step_templates import (
    TEMPLATE_PATH_ENV,
    builtin_spec,
    clear_cache,
    resolve,
)
from xenium_viewer.utils.step_templates import loader
from xenium_viewer.utils.step_templates.overrides import (
    read_manifest,
    remove_override,
    save_override,
)
from xenium_viewer.utils.step_templates.spec import BlockSpec

LEIDEN = "clustering.leiden"
CUSTOM_SCALE = "\nsc.pp.scale(adata_leiden, max_value=99)"


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv(TEMPLATE_PATH_ENV, str(tmp_path))
    clear_cache()
    yield tmp_path
    clear_cache()


def ship_a_new_version(monkeypatch, template_id: str, *,
                       blocks: dict = None, schema_version: int = None):
    """Patch the builtin registry as if a release changed the template."""
    original = dict(loader._builtin_registry())
    spec = original[template_id]
    new_blocks = dict(spec.blocks)
    for name, text in (blocks or {}).items():
        new_blocks[name] = BlockSpec(name=name, text=text,
                                     editable=spec.blocks[name].editable)
    updated = replace(spec, blocks=new_blocks,
                      schema_version=schema_version or spec.schema_version)
    monkeypatch.setattr(loader, "_builtin_registry",
                        lambda: {**original, template_id: updated})
    clear_cache()


# ── the manifest ─────────────────────────────────────────────────────────────

def test_saving_records_what_was_forked_from(home):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})

    record = read_manifest(home)[LEIDEN]
    assert record["based_on_schema"] == builtin_spec(LEIDEN).schema_version
    assert record["forked_from"] == {
        "scale": builtin_spec(LEIDEN).blocks["scale"].hash()
    }, "the hash must be of the SHIPPED text, not the user's edit"


def test_reverting_forgets_the_fork_record(home):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    remove_override(LEIDEN)
    assert LEIDEN not in read_manifest(home)


def test_a_corrupt_manifest_costs_the_warning_not_the_override(home):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    (home / "overrides.json").write_text("{ not json")

    resolved = resolve(LEIDEN)
    assert resolved.is_customised, "losing bookkeeping must not lose the edit"
    assert not resolved.needs_review


# ── the untouched case, which needs no logic ─────────────────────────────────

def test_an_untouched_block_silently_adopts_the_new_version(home, monkeypatch):
    """The dpkg 'unmodified conffile' case — structural, not detected."""
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    ship_a_new_version(monkeypatch, LEIDEN, blocks={
        "hvg": "\nsc.pp.highly_variable_genes(adata_leiden, n_top_genes=$n_top_genes, flavor='cell_ranger')",
    })

    resolved = resolve(LEIDEN)
    assert "cell_ranger" in resolved.spec.blocks["hvg"].text, (
        "a block the user never customised must track the new release"
    )
    assert resolved.spec.blocks["scale"].text == CUSTOM_SCALE
    assert not resolved.needs_review, "nothing the user changed was affected"


# ── the conflict case ────────────────────────────────────────────────────────

def test_a_customised_block_that_upstream_changed_needs_review(home, monkeypatch):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    ship_a_new_version(monkeypatch, LEIDEN, blocks={
        "scale": "\nsc.pp.scale(adata_leiden, max_value=10, zero_center=False)",
    })

    resolved = resolve(LEIDEN)
    assert resolved.stale_blocks == ("scale",)
    assert resolved.needs_review


def test_a_conflicted_override_still_applies(home, monkeypatch):
    """Flagged, not dropped. Silently reverting a user's method is worse."""
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    ship_a_new_version(monkeypatch, LEIDEN, blocks={
        "scale": "\nsc.pp.scale(adata_leiden, max_value=10, zero_center=False)",
    })

    resolved = resolve(LEIDEN)
    assert resolved.is_customised
    assert not resolved.rejected
    assert "max_value=99" in resolved.spec.blocks["scale"].text


def test_a_conflict_that_no_longer_validates_is_deactivated(home, monkeypatch):
    """Review is for edits that still work; a broken one is refused outright.

    Simulates the upgrade that most deserves refusing: a release drops a
    parameter, and the user's forked block still references it. Rendering would
    fail outright, so there is nothing to review — the shipped template runs.
    """
    save_override(LEIDEN, {
        "scale": "\nsc.pp.scale(adata_leiden, max_value=$n_pcs)",
    })

    original = dict(loader._builtin_registry())
    spec = original[LEIDEN]
    without_n_pcs = replace(
        spec, params=tuple(p for p in spec.params if p.name != "n_pcs"))
    monkeypatch.setattr(loader, "_builtin_registry",
                        lambda: {**original, LEIDEN: without_n_pcs})
    clear_cache()

    resolved = resolve(LEIDEN)
    assert resolved.rejected
    assert not resolved.is_customised
    assert any("does not declare" in p.message for p in resolved.problems)


def test_a_schema_bump_is_flagged(home, monkeypatch):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    ship_a_new_version(monkeypatch, LEIDEN, schema_version=2)

    resolved = resolve(LEIDEN)
    assert resolved.schema_moved
    assert resolved.needs_review


def test_re_saving_clears_the_review_flag(home, monkeypatch):
    """Saving is the act of reviewing: it re-records what was forked from."""
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    ship_a_new_version(monkeypatch, LEIDEN, blocks={
        "scale": "\nsc.pp.scale(adata_leiden, max_value=10, zero_center=False)",
    })
    assert resolve(LEIDEN).needs_review

    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    assert not resolve(LEIDEN).needs_review


def test_taking_the_new_default_ends_the_customisation(home, monkeypatch):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    new_text = "\nsc.pp.scale(adata_leiden, max_value=10, zero_center=False)"
    ship_a_new_version(monkeypatch, LEIDEN, blocks={"scale": new_text})

    # "Take new default" == save the shipped text, which resolves to no override.
    save_override(LEIDEN, {"scale": new_text})

    resolved = resolve(LEIDEN)
    assert not resolved.is_customised
    assert not resolved.needs_review
    assert LEIDEN not in read_manifest(home)


def test_an_override_with_no_manifest_entry_is_not_flagged(home):
    """A hand-placed .tmpl has no fork record; prompting on it teaches nothing."""
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    (home / "overrides.json").unlink()
    clear_cache()

    resolved = resolve(LEIDEN)
    assert resolved.is_customised
    assert not resolved.needs_review


# ── the review UI ────────────────────────────────────────────────────────────

def _tab(home):
    from types import SimpleNamespace
    from xenium_viewer.tabs.tab_templates import build_tab
    widget, _ = build_tab(SimpleNamespace(state={}))
    return widget


def _row(widget, template_id):
    from qtpy.QtCore import Qt
    from qtpy.QtWidgets import QTreeWidget
    tree = widget.findChild(QTreeWidget)
    for g in range(tree.topLevelItemCount()):
        group = tree.topLevelItem(g)
        for c in range(group.childCount()):
            if group.child(c).data(0, Qt.UserRole) == template_id:
                return tree, group.child(c)
    raise AssertionError(f"{template_id} not listed")


def _button(widget, text):
    from qtpy.QtWidgets import QPushButton
    return next(b for b in widget.findChildren(QPushButton)
                if b.text().replace("&", "") == text)


def test_a_stale_override_is_badged_for_review(qapp, home, monkeypatch):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    ship_a_new_version(monkeypatch, LEIDEN, blocks={
        "scale": "\nsc.pp.scale(adata_leiden, max_value=10, zero_center=False)",
    })

    widget = _tab(home)
    _, item = _row(widget, LEIDEN)
    assert item.text(1) == "⚠ review"


def test_the_diff_shows_both_versions(qapp, home, monkeypatch):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    ship_a_new_version(monkeypatch, LEIDEN, blocks={
        "scale": "\nsc.pp.scale(adata_leiden, max_value=10, zero_center=False)",
    })

    from xenium_viewer.tabs.tab_templates import _diff_text
    diff = _diff_text(resolve(LEIDEN))
    # Each on its own line: block text has no trailing newline, and keeping line
    # endings ran the '-' and '+' together in the one pane whose job is showing
    # the user exactly which line changed.
    lines = diff.splitlines()
    assert "-sc.pp.scale(adata_leiden, max_value=99)" in lines
    assert "+sc.pp.scale(adata_leiden, max_value=10, zero_center=False)" in lines


def test_take_new_default_keeps_other_customisations(qapp, home, monkeypatch):
    """Accepting one block's update must not silently discard the others."""
    other = "\nsc.pp.pca(adata_leiden, n_comps=7)"
    save_override(LEIDEN, {"scale": CUSTOM_SCALE, "pca": other})
    new_scale = "\nsc.pp.scale(adata_leiden, max_value=10, zero_center=False)"
    ship_a_new_version(monkeypatch, LEIDEN, blocks={"scale": new_scale})

    widget = _tab(home)
    tree, item = _row(widget, LEIDEN)
    tree.setCurrentItem(item)
    _button(widget, "Take new default for changed blocks").click()

    resolved = resolve(LEIDEN)
    assert not resolved.needs_review
    assert resolved.spec.blocks["scale"].text == new_scale, "took the update"
    assert resolved.spec.blocks["pca"].text == other, "kept the other edit"
    assert resolved.changed_blocks() == ["pca"]


def test_the_take_button_is_hidden_when_there_is_nothing_to_review(qapp, home):
    save_override(LEIDEN, {"scale": CUSTOM_SCALE})
    widget = _tab(home)
    tree, item = _row(widget, LEIDEN)
    tree.setCurrentItem(item)
    assert not _button(widget, "Take new default for changed blocks").isVisible()


# ── making it visible downstream ─────────────────────────────────────────────

def _graph_with(origin: str, **kw):
    from xenium_viewer.utils.prov_graph import ProvGraph
    graph = ProvGraph()
    graph.upsert("preamble", "import scanpy as sc", kind="setup")
    graph.upsert("clustering:k", "x = 1", deps=["preamble"],
                 label="Clustering: k", template_origin=origin, **kw)
    return graph


def test_a_stock_notebook_gets_no_banner():
    """Don't clutter the ordinary case with a note saying nothing happened."""
    from xenium_viewer.utils import notebook_export

    graph = _graph_with("builtin")
    assert notebook_export.customisation_banner(graph) is None
    assert all(t != "markdown" or "customised" not in s.lower()
               for t, s in notebook_export.graph_to_cells(graph))


def test_a_customised_notebook_says_so_at_the_top():
    from xenium_viewer.utils import notebook_export

    graph = _graph_with("user", template_id="clustering.leiden",
                        template_hash="abc123def456789")
    cells = notebook_export.graph_to_cells(graph)

    kind, source = cells[0]
    assert kind == "markdown"
    assert "did not use the shipped templates" in source
    assert "clustering:k" in source and "clustering.leiden" in source
    assert "abc123def456" in source


def test_the_banner_does_not_disturb_the_code_cells():
    """The replay test asserts code cells are node.code verbatim."""
    from xenium_viewer.utils import notebook_export

    graph = _graph_with("user", template_id="clustering.leiden")
    code = [s for t, s in notebook_export.graph_to_cells(graph) if t == "code"]
    assert code == [graph.get(n).code for n in graph.topo_sort()]


def test_a_hand_edited_cell_is_called_out_specifically():
    """It is the one origin whose code may not describe what produced the result."""
    from xenium_viewer.utils import notebook_export

    banner = notebook_export.customisation_banner(_graph_with("hand-edited"))
    assert "not re-executed" in banner


def test_the_verification_report_states_whether_templates_were_stock():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "verify_notebook",
        Path(__file__).resolve().parent.parent / "scripts" / "verify_notebook.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    stock = module.template_provenance(_graph_with("builtin"))
    assert stock["stock_templates"] is True
    assert stock["n_customised"] == 0

    custom = module.template_provenance(
        _graph_with("user", template_id="clustering.leiden"))
    assert custom["stock_templates"] is False
    assert [s["node"] for s in custom["customised"]] == ["clustering:k"]
    assert len(custom["steps"]) == 2
