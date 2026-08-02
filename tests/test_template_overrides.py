"""User overrides: what is applied, what is refused, and what is said about it.

Three properties carry the safety of this feature, and each has been a real bug
in some other tool:

* **A broken override must never make the viewer unlaunchable.** The user would
  have no way in to fix the file that is breaking it.
* **A refused override must be loud.** Silently running the shipped template
  while the user believes their edit is in effect is worse than either extreme.
* **Resolution is per block.** A user who customises one block keeps receiving
  upstream fixes to the others. Whole-file override freezes the whole template
  at the version it was forked from, which is how someone quietly misses a
  correctness fix.

``conftest.py`` empties ``XENIUM_VIEWER_TEMPLATE_PATH`` for the suite, so every
test here sets it explicitly to a tmp_path.
"""

from __future__ import annotations

import pytest

from xenium_viewer.utils.step_templates import (
    ERROR,
    TEMPLATE_PATH_ENV,
    builtin_spec,
    clear_cache,
    resolve,
    set_overrides_enabled,
    validate,
)
from xenium_viewer.utils.step_templates.namespace import EXECUTOR_BASE_NAMES

LEIDEN = "clustering.leiden"


@pytest.fixture
def override_dir(tmp_path, monkeypatch):
    """A search path containing only *tmp_path*, with caches cleared around it."""
    monkeypatch.setenv(TEMPLATE_PATH_ENV, str(tmp_path))
    clear_cache()
    yield tmp_path
    clear_cache()


def write_override(directory, template_id: str, blocks: dict, *,
                   schema_version: int = 1, extra_header: str = "") -> None:
    """Write a partial override supplying only *blocks*."""
    spec = builtin_spec(template_id)
    lines = [
        "# xenium-viewer template",
        f"# id: {template_id}",
        f"# schema-version: {schema_version}",
    ]
    if extra_header:
        lines.append(extra_header)
    lines.append("")
    for name, text in blocks.items():
        lines.append(f"#--- block {name}")
        lines.append(text.lstrip("\n"))
        lines.append("")
    (directory / f"{template_id}.tmpl").write_text("\n".join(lines))
    return spec


# ── the happy path ───────────────────────────────────────────────────────────

def test_no_override_resolves_to_the_shipped_template(override_dir):
    resolved = resolve(LEIDEN)
    assert not resolved.is_customised
    assert resolved.origin == "builtin"
    assert resolved.problems == ()
    assert resolved.spec.blocks == builtin_spec(LEIDEN).blocks


def test_an_override_replaces_only_the_block_it_supplies(override_dir):
    """The property that makes upgrades survivable."""
    builtin = builtin_spec(LEIDEN)
    write_override(override_dir, LEIDEN, {
        "scale": "\nsc.pp.scale(adata_leiden, max_value=99)",
    })

    resolved = resolve(LEIDEN)
    assert resolved.is_customised
    assert resolved.changed_blocks() == ["scale"]
    assert "max_value=99" in resolved.spec.blocks["scale"].text
    for name in ("head", "hvg", "pca", "tail"):
        assert resolved.spec.blocks[name].text == builtin.blocks[name].text, (
            f"{name} was not customised and must keep tracking the shipped text"
        )


def test_a_partial_override_is_reported_as_blended(override_dir):
    write_override(override_dir, LEIDEN, {
        "scale": "\nsc.pp.scale(adata_leiden, max_value=99)",
    })
    assert resolve(LEIDEN).origin == "user+builtin"


def test_the_resolved_text_is_what_would_run(override_dir):
    from xenium_viewer.utils.step_templates import resolved_text

    write_override(override_dir, LEIDEN, {
        "scale": "\nsc.pp.scale(adata_leiden, max_value=99)",
    })
    text = resolved_text(LEIDEN, ["head", "scale", "pca", "tail"])
    assert "max_value=99" in text
    assert "adata_leiden = adata_norm.copy()" in text     # untouched head


# ── refusal ──────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("bad_block,reason", [
    ("\nsc.pp.scale(adata_leiden,", "not valid Python"),
    ("\nsc.pp.scale(adata_leiden, max_value=$nonsense)", "does not declare"),
    ("\nsc.pp.scale(mystery_object)", "reads"),
    ("\nfrom xenium_viewer.utils import gene_analysis", "xenium_viewer"),
])
def test_an_invalid_override_is_refused_and_the_builtin_is_used(
        override_dir, bad_block, reason):
    builtin = builtin_spec(LEIDEN)
    write_override(override_dir, LEIDEN, {"scale": bad_block})

    resolved = resolve(LEIDEN)
    assert resolved.rejected
    assert not resolved.is_customised
    assert resolved.spec.blocks == builtin.blocks, "must fall back to shipped text"
    assert any(reason in p.message for p in resolved.problems), (
        f"expected a problem mentioning {reason!r}, got "
        f"{[p.message for p in resolved.problems]}"
    )


def test_an_unknown_block_is_refused(override_dir):
    """Blocks are chosen by name, so an unknown one would never render."""
    write_override(override_dir, LEIDEN, {"not_a_block": "\nx = 1"})
    resolved = resolve(LEIDEN)
    assert resolved.rejected
    assert any("does not have" in p.message for p in resolved.problems)


def test_a_frozen_block_cannot_be_customised(override_dir):
    """The Arrow shim: no gate could warn when editing it breaks silently."""
    tid = "genes.cnv_infercnv"
    write_override(override_dir, tid, {"arrow_shim": "\n_old_infer = None"})
    resolved = resolve(tid)
    assert resolved.rejected
    assert any("workaround" in p.message for p in resolved.problems)


def test_an_override_from_a_newer_release_is_a_hard_stop(override_dir):
    """Guessing at a contract this build does not know is worse than refusing."""
    write_override(override_dir, LEIDEN,
                   {"scale": "\nsc.pp.scale(adata_leiden)"}, schema_version=99)
    resolved = resolve(LEIDEN)
    assert resolved.rejected
    assert any("newer release" in p.message for p in resolved.problems)


def test_dropping_a_required_param_is_refused(override_dir):
    """The silently-wrong-answer case, and the reason this gate exists.

    A template that stops mentioning ``$n_neighbors`` still runs, still reports
    success, and quietly ignores the setting the user chose.
    """
    builtin = builtin_spec(LEIDEN)
    tail = builtin.blocks["tail"].text.replace(
        "n_neighbors=$n_neighbors, ", "")
    write_override(override_dir, LEIDEN, {"tail": tail})

    resolved = resolve(LEIDEN)
    assert resolved.rejected
    assert any("silently ignored" in p.message for p in resolved.problems)


def test_an_unreadable_override_does_not_raise(override_dir):
    """A directory where a file should be — the viewer must still start."""
    (override_dir / f"{LEIDEN}.tmpl").mkdir()
    resolved = resolve(LEIDEN)
    assert resolved.spec.blocks == builtin_spec(LEIDEN).blocks
    assert resolved.rejected


def test_garbage_in_the_config_dir_does_not_raise(override_dir):
    (override_dir / f"{LEIDEN}.tmpl").write_text("this is not a template at all")
    resolved = resolve(LEIDEN)
    assert resolved.spec.blocks == builtin_spec(LEIDEN).blocks


# ── the off switches ─────────────────────────────────────────────────────────

def test_no_user_templates_ignores_a_valid_override(override_dir):
    """The first thing to try when a result is in doubt."""
    write_override(override_dir, LEIDEN, {
        "scale": "\nsc.pp.scale(adata_leiden, max_value=99)",
    })
    assert resolve(LEIDEN).is_customised

    set_overrides_enabled(False)
    try:
        assert not resolve(LEIDEN).is_customised
        assert resolve(LEIDEN).origin == "builtin"
    finally:
        set_overrides_enabled(True)


def test_the_suite_itself_runs_with_overrides_disabled():
    """conftest empties the search path; a dev's own config must not leak in."""
    import os

    from xenium_viewer.utils.step_templates import search_path
    assert os.environ.get(TEMPLATE_PATH_ENV) == ""
    assert search_path() == []


def test_builtin_accessors_never_see_an_override(override_dir):
    """Why the template-pinning tests are immune by construction, not by care."""
    from xenium_viewer.utils.step_templates import builtin_assemble

    write_override(override_dir, LEIDEN, {
        "scale": "\nsc.pp.scale(adata_leiden, max_value=99)",
    })
    assert resolve(LEIDEN).is_customised          # the override *is* active
    assert "max_value=99" not in builtin_assemble(LEIDEN, ["head", "scale"])
    assert "max_value=10" in builtin_assemble(LEIDEN, ["head", "scale"])


# ── validate() on its own ────────────────────────────────────────────────────

def test_a_shipped_template_validates_clean():
    """Guard the guard: a gate that rejected everything would pass the tests above."""
    for spec in (builtin_spec(t) for t in ("clustering.leiden", "roi.deg")):
        problems = validate(spec, builtin=spec,
                            available=EXECUTOR_BASE_NAMES | spec.requires)
        assert [p for p in problems if p.severity == ERROR] == []


def test_a_bare_dollar_is_reported_clearly():
    from xenium_viewer.utils.step_templates import placeholders
    from xenium_viewer.utils.steps import StepError

    with pytest.raises(StepError, match="literal dollar"):
        placeholders("cost = 5 $ each")


# ── saying so out loud ───────────────────────────────────────────────────────

def test_a_rejected_override_is_reported_not_swallowed(override_dir):
    """The failure this feature most needs to avoid causing.

    A user who edited a template and is silently getting the shipped one will
    attribute the resulting numbers to their own method.
    """
    from xenium_viewer.utils import reporting

    reporting.clear_template_rejections()
    write_override(override_dir, LEIDEN, {"scale": "\nsc.pp.scale(adata_leiden,"})

    assert resolve(LEIDEN).rejected
    rejections = reporting.template_rejections()
    assert LEIDEN in rejections
    assert any("not valid Python" in line for line in rejections[LEIDEN])
    assert "customised template" in reporting.failure_summary()
    reporting.clear_template_rejections()


def test_a_valid_override_reports_nothing(override_dir):
    from xenium_viewer.utils import reporting

    reporting.clear_template_rejections()
    write_override(override_dir, LEIDEN, {
        "scale": "\nsc.pp.scale(adata_leiden, max_value=99)",
    })
    assert resolve(LEIDEN).is_customised
    assert reporting.template_rejections() == {}


def test_reporting_failure_never_breaks_resolution(override_dir, monkeypatch):
    """Resolution runs at analysis time; a reporting bug must not take it down."""
    from xenium_viewer.utils import reporting

    monkeypatch.setattr(reporting, "report_template_rejected",
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    write_override(override_dir, LEIDEN, {"scale": "\nsc.pp.scale(adata_leiden,"})
    assert resolve(LEIDEN).spec.blocks == builtin_spec(LEIDEN).blocks


# ── the provenance stamp ─────────────────────────────────────────────────────

def test_a_stock_run_is_stamped_builtin(override_dir):
    from xenium_viewer.utils.step_templates import step_template

    stamp = step_template(LEIDEN, ["head", "tail"])
    assert stamp["template_id"] == LEIDEN
    assert stamp["template_origin"] == "builtin"
    assert stamp["template_hash"] == builtin_spec(LEIDEN).hash_of(["head", "tail"])
    assert "sc.tl.leiden(" in stamp["template"]


def test_a_customised_run_is_stamped_and_hashes_differently(override_dir):
    from xenium_viewer.utils.step_templates import step_template

    stock = step_template(LEIDEN, ["head", "scale", "pca", "tail"])
    write_override(override_dir, LEIDEN, {
        "scale": "\nsc.pp.scale(adata_leiden, max_value=99)",
    })
    clear_cache()
    custom = step_template(LEIDEN, ["head", "scale", "pca", "tail"])

    assert custom["template_origin"] == "user+builtin"
    assert custom["template_hash"] != stock["template_hash"], (
        "the hash is what lets a reader tell a stock run from a customised one"
    )
    assert "max_value=99" in custom["template"]


def test_the_stamp_reaches_the_provenance_node(override_dir):
    """End to end: a customised template is visible in the recorded graph."""
    from xenium_viewer.utils.step_templates import step_template
    from xenium_viewer.utils.steps import Step, StepExecutor

    write_override(override_dir, "roi.polygons", {
        "main": "\nroi_polygons = [np.array(_p) for _p in $polygons]  # customised",
    })
    clear_cache()

    import numpy as np
    ex = StepExecutor(namespace={"np": np})
    ex.run(Step(id="rois", **step_template("roi.polygons", ["main"]),
                params={"polygons": [[[0, 0], [1, 1]]]}, kind="setup"))

    node = ex.graph.get("rois")
    assert node.template_origin == "user"
    assert node.template_id == "roi.polygons"
    assert node.template_hash
    assert "# customised" in node.code, "the recorded code is the customised code"


# ── the configuration a real user actually has ───────────────────────────────
#
# Everything above sets XENIUM_VIEWER_TEMPLATE_PATH, and conftest sets it for
# the whole suite. That is what keeps a developer's own customisations out of
# the tests — and it meant the *default* path, with the variable unset, was
# never executed by anything. It was infinitely recursive, and the viewer would
# not start. These tests run with the variable deleted.

@pytest.fixture
def no_env(monkeypatch, tmp_path):
    """As a real user has it: no XENIUM_VIEWER_TEMPLATE_PATH at all."""
    from xenium_viewer.utils.step_templates import loader

    monkeypatch.delenv(TEMPLATE_PATH_ENV, raising=False)
    # Redirect the platform config dir so nothing touches the real one.
    monkeypatch.setattr(loader, "_default_user_dir", lambda: tmp_path)
    clear_cache()
    yield tmp_path
    clear_cache()


def test_the_search_path_resolves_with_no_env_var(no_env):
    from xenium_viewer.utils.step_templates import search_path, user_template_dir

    assert search_path() == [no_env]
    assert user_template_dir() == no_env


def test_resolving_works_with_no_env_var(no_env):
    """The launch path. This raised RecursionError and the viewer never opened."""
    resolved = resolve(LEIDEN)
    assert not resolved.is_customised
    assert resolved.spec.blocks == builtin_spec(LEIDEN).blocks


def test_every_template_resolves_with_no_env_var(no_env):
    """The tab populates by resolving all of them, so all of them must work."""
    from xenium_viewer.utils.step_templates import builtin_ids

    for template_id in builtin_ids():
        assert resolve(template_id).spec.blocks


def test_an_override_still_applies_with_no_env_var(no_env):
    """The default location is a real location, not just a non-crashing one."""
    from xenium_viewer.utils.step_templates.overrides import save_override

    save_override(LEIDEN, {"scale": "\nsc.pp.scale(adata_leiden, max_value=99)"})
    assert (no_env / f"{LEIDEN}.tmpl").is_file()
    assert resolve(LEIDEN).is_customised


def test_disabling_overrides_works_with_no_env_var(no_env):
    resolve(LEIDEN)
    set_overrides_enabled(False)
    try:
        from xenium_viewer.utils.step_templates import search_path
        assert search_path() == []
        assert not resolve(LEIDEN).is_customised
    finally:
        set_overrides_enabled(True)


def test_the_two_entry_points_do_not_delegate_to_each_other():
    """Source guard: the shape that made this recursive.

    Both must read the configuration themselves. One calling the other reads
    naturally, passes every test that sets the variable, and hangs the viewer
    for everyone who does not.
    """
    import ast
    import inspect

    from xenium_viewer.utils.step_templates import loader

    for name, forbidden in (("search_path", "user_template_dir"),
                            ("user_template_dir", "search_path")):
        tree = ast.parse(inspect.getsource(getattr(loader, name)))
        called = {n.func.id for n in ast.walk(tree)
                  if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
        assert forbidden not in called, (
            f"{name}() calls {forbidden}(); with the env var unset these two "
            f"recurse forever"
        )
