"""Whether a template is safe to run — the gate a user override must pass.

This is deliberately **not a security boundary**, and saying so plainly matters
more than it might seem: a user can already execute arbitrary Python in the
Notebook tab's free-form cells, and the notebook this system exports is a file
of code they are meant to run. Treating template validation as a sandbox would
invite the wrong kind of review and give false assurance.

What it *is* is a correctness gate, aimed at one failure mode above all others.
A template that crashes is fine — the user sees a ``StepError`` and fixes it.
A template that **runs and quietly produces the wrong number** is the one that
ends up in a paper. So the checks are weighted toward that: a required parameter
the template no longer mentions is a hard stop, because the analysis would run
with that setting having no effect and look entirely successful.

Everything here is static — render with sentinel values, parse, walk the AST.
Nothing is executed. That is what lets the Templates tab offer a "Validate"
button that answers in milliseconds instead of after a ten-minute run.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from string import Template
from typing import Iterable, Optional

from palms.utils.step_templates.spec import TemplateSpec
from palms.utils.steps import StepError, free_names, render

ERROR = "error"
WARNING = "warning"


@dataclass(frozen=True)
class Problem:
    """One validation finding. ``ERROR`` blocks activation; ``WARNING`` does not."""

    message: str
    severity: str = ERROR
    block: Optional[str] = None
    line: Optional[int] = None

    def __str__(self) -> str:
        where = f" [{self.block}]" if self.block else ""
        at = f" (line {self.line})" if self.line else ""
        return f"{self.severity}{where}{at}: {self.message}"


def placeholders(text: str) -> set[str]:
    """Names ``Template.substitute`` would demand from *text*.

    Hand-rolled off ``Template.pattern`` rather than using
    ``Template.get_identifiers()``, which is Python 3.11+ while the package
    declares ``requires-python = ">=3.10"``. The stdlib method would work in the
    dev env and fail for anyone on the declared floor — the worst kind of
    version bug, because local testing never sees it.

    A malformed placeholder (a bare ``$``) raises here rather than surfacing
    later as ``Template.substitute``'s much vaguer ValueError.
    """
    names: set[str] = set()
    for match in Template.pattern.finditer(text):
        if match.group("invalid") is not None:
            raise StepError(
                f"malformed placeholder {match.group(0)!r}; write '$$' for a "
                f"literal dollar sign"
            )
        name = match.group("named") or match.group("braced")
        if name is not None:
            names.add(name)
    return names


def module_level_bindings(code: str) -> set[str]:
    """Names bound by a top-level statement of *code*.

    Conservative in the safe direction: a name bound only inside an ``if`` body
    counts as bound, so this can pass something that fails at run time — where
    ``StepExecutor.run`` catches it — but never fails something that works. A
    false accusation would be worse: it would block a template that is fine.
    """
    bound: set[str] = set()
    for node in ast.parse(code).body:
        for sub in ast.walk(node):
            if isinstance(sub, ast.Name) and isinstance(sub.ctx, ast.Store):
                bound.add(sub.id)
            elif isinstance(sub, (ast.Import, ast.ImportFrom)):
                for alias in sub.names:
                    bound.add((alias.asname or alias.name).split(".")[0])
            elif isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                bound.add(sub.name)
            elif isinstance(sub, ast.withitem) and isinstance(sub.optional_vars, ast.Name):
                bound.add(sub.optional_vars.id)
    return bound


def validate(
    spec: TemplateSpec,
    *,
    builtin: Optional[TemplateSpec] = None,
    available: Optional[Iterable[str]] = None,
) -> list[Problem]:
    """Every problem with *spec*. Empty means it is safe to activate.

    *builtin* is the shipped spec to compare against (structure, frozen blocks,
    schema version); omit it to validate a builtin against itself. *available*
    is the set of names the executor namespace guarantees.
    """
    problems: list[Problem] = []
    available = set(available or ())

    if builtin is not None:
        problems += _check_structure(spec, builtin)
        # A structural mismatch makes per-assembly checks meaningless — they
        # would report a cascade of consequences instead of the one cause.
        if any(p.severity == ERROR for p in problems):
            return problems

    problems += _check_params(spec)
    for assembly in spec.assemblies:
        problems += _check_assembly(spec, assembly, available)

    for block in spec.blocks.values():
        if "palms" in block.text:
            problems.append(Problem(
                "reaches back into palms; the exported notebook has to "
                "run without this package installed",
                block=block.name,
            ))
    return problems


# ── individual checks ────────────────────────────────────────────────────────

def _check_structure(spec: TemplateSpec, builtin: TemplateSpec) -> list[Problem]:
    problems = []
    if spec.schema_version > builtin.schema_version:
        problems.append(Problem(
            f"declares schema-version {spec.schema_version}, but this version "
            f"of the viewer only understands {builtin.schema_version}. It was "
            f"written for a newer release; upgrade rather than guessing."
        ))
    unknown = sorted(set(spec.blocks) - set(builtin.blocks))
    if unknown:
        problems.append(Problem(
            f"defines block(s) {unknown} that this template does not have. "
            f"Blocks are selected by name, so these would never render. "
            f"Known blocks: {sorted(builtin.blocks)}"
        ))
    for name, block in builtin.blocks.items():
        if block.editable or name not in spec.blocks:
            continue
        if spec.blocks[name].text != block.text:
            problems.append(Problem(
                "this block is a version-compatibility workaround and cannot be "
                "customised — no check could warn you when it silently stops "
                "working against a future dependency release",
                block=name,
            ))
    return problems


def _check_params(spec: TemplateSpec) -> list[Problem]:
    problems = []
    all_text = "".join(b.text for b in spec.blocks.values())
    try:
        used = placeholders(all_text)
    except StepError as exc:
        return [Problem(str(exc))]

    undeclared = sorted(used - spec.param_names)
    if undeclared:
        problems.append(Problem(
            f"uses ${{{'}, ${'.join(undeclared)}}}, which the header does not "
            f"declare. Add them to '# params:' or remove the references — "
            f"rendering fails outright on an undeclared placeholder."
        ))
    unused = sorted(spec.param_names - used)
    if unused:
        problems.append(Problem(
            f"declares params {unused} that no block uses", severity=WARNING,
        ))
    return problems


def _check_assembly(spec: TemplateSpec, assembly, available: set) -> list[Problem]:
    label = "+".join(assembly)
    try:
        text = spec.assemble(assembly)
    except KeyError as exc:
        return [Problem(f"assembly {label}: {exc}")]

    # The check that most directly prevents a silently-wrong result: a required
    # param the template never mentions means the analysis runs with that
    # setting having no effect, and reports success.
    missing = sorted(n for n in spec.required_params if f"${n}" not in text)
    if missing:
        return [Problem(
            f"assembly {label} never uses required param(s) {missing}, so those "
            f"settings would be silently ignored. Mark them optional with '?' "
            f"in the header if that is genuinely intended."
        )]

    try:
        code = render(text, spec.synth_params())
    except StepError as exc:
        return [Problem(f"assembly {label}: {exc}")]

    try:
        tree = ast.parse(code)
    except SyntaxError as exc:
        return [Problem(f"assembly {label} is not valid Python: {exc.msg}",
                        line=exc.lineno)]
    if not tree.body:
        return [Problem(f"assembly {label} renders to nothing executable")]

    problems = []
    unresolved = sorted(free_names(code) - available - set(spec.requires))
    if unresolved:
        problems.append(Problem(
            f"assembly {label} reads {unresolved}, which nothing provides. "
            f"Either a dependency step binds them — add them to '# requires:' — "
            f"or this is a typo that would replay as a NameError."
        ))
    unbound = sorted(set(spec.outputs) - module_level_bindings(code))
    if unbound:
        problems.append(Problem(
            f"assembly {label} declares output(s) {unbound} that it never "
            f"binds; the viewer reads results by those names."
        ))
    return problems
