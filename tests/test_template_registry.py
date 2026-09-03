"""Every shipped template, in every legal assembly, is well-formed code.

``check_step`` and ``free_names`` have existed in ``utils/steps.py`` since the
Step system landed, described there as "the CI template lint" — but they were
only ever called from five hand-written tests covering five templates. This
module is that lint, applied to the registry: 14 templates and every assembly
each of them declares, which is 40 concrete renderings rather than 5.

What it establishes, per rendering:

* it parses as Python;
* every name it *reads* is one the executor namespace guarantees or its own
  declared ``requires`` — so no template can quietly depend on a leftover from a
  previous step in a long session, which would replay as a ``NameError``;
* every name it declares as an output is actually bound at module level;
* it never reaches back into ``palms``, because the notebook has to run
  without this package installed.

Pure Python — no Qt, no napari. The block *selection* logic lives in the tab
modules, but the legal sequences are declared in the ``.tmpl`` headers, which is
what lets this run without importing a tab.
"""

from __future__ import annotations

import ast

import pytest

from palms.utils.step_templates import (
    EXECUTOR_BASE_NAMES,
    builtin_ids,
    builtin_spec,
)
from palms.utils.steps import Step, check_step

#: (template id, assembly) for every legal rendering the registry declares.
ALL_ASSEMBLIES = [
    (tid, assembly)
    for tid in builtin_ids()
    for assembly in builtin_spec(tid).assemblies
]


def _ids(case):
    tid, assembly = case
    return f"{tid}[{'+'.join(assembly)}]"


def _step(tid: str, assembly) -> Step:
    spec = builtin_spec(tid)
    return Step(id=f"test:{tid}", template=spec.assemble(assembly),
                params=spec.synth_params(), template_id=tid)


def _module_level_bindings(code: str) -> set[str]:
    """Names bound by a top-level statement of *code*.

    Conservative in the safe direction: a name bound only inside an ``if`` or a
    loop body counts as bound, so this can pass something that fails at run time
    — where ``StepExecutor.run`` catches it — but it never fails something that
    would have worked.
    """
    bound: set[str] = set()
    for node in ast.parse(code).body:
        for sub in ast.walk(node):
            if isinstance(sub, ast.Name) and isinstance(sub.ctx, ast.Store):
                bound.add(sub.id)
            elif isinstance(sub, (ast.Import, ast.ImportFrom)):
                for alias in sub.names:
                    bound.add((alias.asname or alias.name).split(".")[0])
            elif isinstance(sub, (ast.FunctionDef, ast.ClassDef)):
                bound.add(sub.name)
    return bound


def test_the_registry_is_not_empty():
    """Guard the guard: a loader returning nothing would pass every test below."""
    assert len(builtin_ids()) >= 14
    assert len(ALL_ASSEMBLIES) >= 40


@pytest.mark.parametrize("case", ALL_ASSEMBLIES, ids=_ids)
def test_every_assembly_renders_and_parses(case):
    tid, assembly = case
    code = _step(tid, assembly).render()
    tree = ast.parse(code)
    assert tree.body, f"{tid} rendered to nothing executable"


@pytest.mark.parametrize("case", ALL_ASSEMBLIES, ids=_ids)
def test_every_assembly_is_self_contained(case):
    """The promoted ``check_step`` lint, over the whole registry."""
    tid, assembly = case
    spec = builtin_spec(tid)
    available = EXECUTOR_BASE_NAMES | spec.requires
    missing = check_step(_step(tid, assembly), available)
    assert missing == set(), (
        f"{tid}[{'+'.join(assembly)}] reads {sorted(missing)}, which neither the "
        f"executor namespace nor its declared 'requires' provides. Either add "
        f"them to the .tmpl's '# requires:' line (if a dependency step binds "
        f"them) or stop reading them."
    )


@pytest.mark.parametrize("case", ALL_ASSEMBLIES, ids=_ids)
def test_declared_outputs_are_bound(case):
    """A declared output the template never binds is a StepError at run time.

    Catching it here means a broken template fails in CI rather than after a
    ten-minute analysis has already run.
    """
    tid, assembly = case
    spec = builtin_spec(tid)
    if not spec.outputs:
        pytest.skip(f"{tid} declares no outputs")
    bound = _module_level_bindings(_step(tid, assembly).render())
    assert set(spec.outputs) <= bound, (
        f"{tid}[{'+'.join(assembly)}] declares outputs "
        f"{sorted(set(spec.outputs) - bound)} that it does not bind"
    )


#: Templates allowed to import ``palms``, and why. Kept here as well as in each
#: ``.tmpl`` header so that granting the allowance is a change to a shared list a
#: reviewer sees, not a line added to one file. Every one of them is an H&E
#: registration step: those algorithms are this package's own and no scanpy,
#: squidpy or spatialdata call computes them, so the choice is between a cell
#: that imports palms and a cell that states a matrix as a literal — and the
#: literal reproduces nothing, which is what these templates exist to prevent.
_MAY_IMPORT_PALMS = {
    "he.load", "he.coarse_align", "he.nuclei_align", "he.landmark_align",
}


@pytest.mark.parametrize("tid", builtin_ids())
def test_no_template_imports_the_viewer(tid):
    """The notebook replays from raw Xenium output, without this package.

    Declaration rather than prohibition since 2026-09-03: a template may import
    ``palms`` when its header says why, and every other template is guarded
    exactly as before. The two halves are asserted in both directions, because
    an exemption nobody uses and an import nobody declared are each a way for
    this list to stop describing the tree.
    """
    spec = builtin_spec(tid)
    imports = [b.name for b in spec.blocks.values() if "palms" in b.text]
    if tid not in _MAY_IMPORT_PALMS:
        assert not imports, (
            f"{tid}/{imports} reaches back into palms; the exported notebook "
            f"must run without it installed. If the computation genuinely has "
            f"no scverse equivalent, declare it with a '# palms:' header and "
            f"add the id to _MAY_IMPORT_PALMS here."
        )
        assert not spec.palms_reason, (
            f"{tid} declares a '# palms:' reason but is not in _MAY_IMPORT_PALMS"
        )
        return
    assert imports, f"{tid} is exempted but no block imports palms"
    assert spec.palms_reason, (
        f"{tid} imports palms but its header declares no '# palms:' reason; "
        f"the validator is what turns that into a hard stop"
    )


def test_the_palms_allowance_cannot_be_granted_by_a_user_override():
    """A user's own ``.tmpl`` must not be able to authorise a palms import.

    The contract comes from the builtin on merge, so ``palms_reason`` is the
    shipped template's answer to a question a user's text does not get to
    re-answer. Asserted here rather than assumed, because the whole guard would
    be decorative if an override could switch it off.
    """
    from palms.utils.step_templates.loader import _merge, parse_template
    from palms.utils.step_templates.validate import ERROR, validate

    base = builtin_spec("clustering.leiden")
    forged = parse_template(
        "# palms template\n"
        "# id: clustering.leiden\n"
        "# palms: I would like to import whatever I please\n"
        "\n#--- block head\n"
        "from palms.utils.registration import compute_coarse_affine\n"
    )
    merged = _merge(base, forged, source="<forged>")

    assert merged.palms_reason == base.palms_reason == ""
    problems = validate(merged, builtin=base)
    assert any(p.severity == ERROR and "palms" in p.message for p in problems)


@pytest.mark.parametrize("tid", builtin_ids())
def test_every_declared_param_is_used_somewhere(tid):
    """A param no assembly references is dead weight in the contract."""
    spec = builtin_spec(tid)
    all_text = "".join(b.text for b in spec.blocks.values())
    unused = {p.name for p in spec.params if f"${p.name}" not in all_text}
    assert not unused, f"{tid} declares unused params: {sorted(unused)}"


@pytest.mark.parametrize("tid", builtin_ids())
def test_required_params_appear_in_every_assembly(tid):
    """"Required" has to mean something: optional params are marked with '?'.

    This is the check that most directly protects against a future param being
    silently ignored — an analysis that runs with a setting having no effect is
    worse than one that fails.
    """
    spec = builtin_spec(tid)
    for assembly in spec.assemblies:
        text = spec.assemble(assembly)
        missing = {n for n in spec.required_params if f"${n}" not in text}
        assert not missing, (
            f"{tid}[{'+'.join(assembly)}] omits required params {sorted(missing)}; "
            f"mark them optional with '?' in the .tmpl header if that is intended"
        )


@pytest.mark.parametrize("tid", builtin_ids())
def test_no_block_is_comment_only(tid):
    """A block that renders to comments is a step that silently does nothing.

    ``allow_errors=False`` cannot catch it — the cell runs fine. Only NOTE nodes
    are allowed to be prose, and none of them come from a template.
    """
    spec = builtin_spec(tid)
    for block in spec.blocks.values():
        stripped = [ln for ln in block.text.splitlines()
                    if ln.strip() and not ln.strip().startswith("#")]
        assert stripped, f"{tid}/{block.name} contains no executable line"
