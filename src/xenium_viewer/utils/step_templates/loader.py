"""Reading ``.tmpl`` files into :class:`TemplateSpec`.

File format — a header of comment lines, then one or more blocks::

    # xenium-viewer template
    # id: clustering.leiden
    # schema-version: 1
    # requires: sc, adata_norm
    # outputs: leiden_labels
    # params: key:str, resolution:float, n_top_genes:int?
    # doc: Leiden community detection on the normalised copy.

    #--- block head
    adata_leiden = adata_norm.copy()

    #--- block hvg
    sc.pp.highly_variable_genes(adata_leiden, n_top_genes=$n_top_genes)

A ``?`` after a param's type marks it optional (not every assembly uses it).
Everything structural is a comment, so a ``.tmpl`` is valid Python: the syntax
highlighter the Templates tab reuses needs no special case, and a stray edit
shows up as a comment rather than as a parse error somewhere else.

``builtin_*`` reads only files shipped inside this package, via
``importlib.resources``, and never consults any override path. That is load
bearing rather than incidental: the tab modules bind their private template
constants through it, so the six test modules that pin template text stay
immune to whatever a developer has in their own config. A test asserting
override-immunity would be checking a promise; routing through a function that
*cannot* see overrides makes it structural.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field, replace
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Iterable, Optional

import platformdirs

from xenium_viewer.utils.prov_graph import TEMPLATE_BLENDED, TEMPLATE_USER
from xenium_viewer.utils.step_templates.namespace import EXECUTOR_BASE_NAMES
from xenium_viewer.utils.step_templates.spec import (
    BLOCK_MARKER,
    BlockSpec,
    ParamSpec,
    TemplateSpec,
)
from xenium_viewer.utils.step_templates.validate import (
    ERROR,
    WARNING,
    Problem,
    validate,
)

_BUILTIN_PACKAGE = "xenium_viewer.utils.step_templates.builtin"

#: ``name:type`` or ``name:type?`` — the trailing ``?`` means "optional".
_PARAM_RE = re.compile(r"^(?P<name>\w+)\s*:\s*(?P<type>\w+)(?P<optional>\?)?$")


class TemplateError(ValueError):
    """A ``.tmpl`` file is malformed or references something unknown."""


@dataclass(frozen=True)
class ResolvedTemplate:
    """What a template id resolves to, and the honest account of why.

    ``spec`` is always usable — it falls back to ``builtin`` when an override
    cannot be trusted, so callers never have to handle "no template". The
    problems and ``rejected`` flag are how the GUI reports that fallback instead
    of silently pretending the user's edit is in effect.
    """

    spec: TemplateSpec
    builtin: TemplateSpec
    problems: tuple = ()
    path: Optional[Path] = None
    rejected: bool = False

    @property
    def is_customised(self) -> bool:
        return self.path is not None

    @property
    def origin(self) -> str:
        return self.spec.origin

    def hash_of(self, block_names) -> str:
        return self.spec.hash_of(block_names)

    def changed_blocks(self) -> list[str]:
        """Blocks whose text differs from the shipped version."""
        return sorted(
            name for name, block in self.spec.blocks.items()
            if name not in self.builtin.blocks
            or block.text != self.builtin.blocks[name].text
        )


# ── parsing ──────────────────────────────────────────────────────────────────

def _split_header_and_blocks(text: str) -> tuple[list[str], list[tuple[str, str]]]:
    """Return (header comment lines, [(block name, block text), ...]).

    Block text keeps its leading newline, so ``"".join`` of consecutive blocks
    reproduces the original template constants exactly.
    """
    header: list[str] = []
    blocks: list[tuple[str, list[str]]] = []

    for line in text.splitlines():
        if line.startswith(BLOCK_MARKER):
            name = line[len(BLOCK_MARKER):].strip()
            if not name:
                raise TemplateError(f"block marker with no name: {line!r}")
            blocks.append((name, []))
        elif blocks:
            blocks[-1][1].append(line)
        elif line.startswith("#"):
            # Strip the '#' and at most one following space, keeping any further
            # indentation: that indentation is what marks a continuation line,
            # and long 'params:'/'assemblies:' lists rely on wrapping. Stripping
            # it made 'flavor:str, ...' on a continuation look like a new field
            # called 'flavor', silently truncating the param list.
            body = line[1:]
            header.append(body[1:] if body.startswith(" ") else body)
        elif line.strip():
            raise TemplateError(
                f"content before the first {BLOCK_MARKER!r} marker: {line!r}"
            )

    # Each block is stored as "\n" + its lines, and trailing blank lines (the
    # separation between blocks in the file) are dropped — they are formatting,
    # not part of the recorded source.
    out: list[tuple[str, str]] = []
    for name, lines in blocks:
        while lines and not lines[-1].strip():
            lines.pop()
        out.append((name, "\n" + "\n".join(lines)))
    return header, out


def _parse_header(lines: Iterable[str]) -> dict[str, str]:
    """``field: value`` lines; an indented line continues the one above it."""
    fields: dict[str, str] = {}
    key: Optional[str] = None
    for line in lines:
        indented = line[:1].isspace()
        if not indented and ":" in line:
            candidate, _, rest = line.partition(":")
            candidate = candidate.strip().lower().replace("-", "_")
            if candidate.isidentifier():
                key = candidate
                fields[key] = rest.strip()
                continue
        if indented and key is not None and line.strip():
            fields[key] = f"{fields[key]} {line.strip()}".strip()
        else:
            # An unindented line that is not 'key: value' is prose. It must not
            # extend the field above it — a explanatory comment in a saved
            # override was being appended to 'schema-version', which then failed
            # to parse as an int and rejected the user's own file.
            key = None
    return fields


def _csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def _parse_params(value: str) -> tuple[ParamSpec, ...]:
    specs = []
    for entry in _csv(value):
        match = _PARAM_RE.match(entry)
        if match is None:
            raise TemplateError(
                f"malformed param {entry!r}; expected 'name:type' or 'name:type?'"
            )
        specs.append(ParamSpec(
            name=match["name"],
            type=match["type"],
            required=match["optional"] is None,
        ))
    return tuple(specs)


def _parse_assemblies(value: str, blocks: dict, template_id: str) -> tuple:
    """``head+hvg+tail | head+tail`` -> the legal block sequences.

    Declared as data even though the *choosing* stays in Python, because the
    validator and the Templates tab both need to enumerate what can be produced
    without importing a tab module (and therefore Qt). The duplication is
    checkable: a test asserts the Python selectors only ever emit a sequence
    listed here, so the two cannot drift silently.

    A template with a single block sequence need not declare it — that case is
    unambiguous, and spelling it out in every flat file would be noise.
    """
    if not value.strip():
        return (tuple(blocks),)
    out = []
    for alternative in value.split("|"):
        names = [n.strip() for n in alternative.split("+") if n.strip()]
        unknown = [n for n in names if n not in blocks]
        if unknown:
            raise TemplateError(
                f"template {template_id!r} declares an assembly using unknown "
                f"block(s) {unknown} (known: {sorted(blocks)})"
            )
        out.append(tuple(names))
    return tuple(out)


def parse_template(text: str, *, source: Optional[str] = None) -> TemplateSpec:
    """Build a :class:`TemplateSpec` from the contents of a ``.tmpl`` file."""
    header_lines, blocks = _split_header_and_blocks(text)
    fields = _parse_header(header_lines)

    template_id = fields.get("id")
    if not template_id:
        raise TemplateError(f"{source or 'template'} has no '# id:' header")
    if not blocks:
        raise TemplateError(f"template {template_id!r} declares no blocks")

    seen: set[str] = set()
    block_specs: dict[str, BlockSpec] = {}
    for name, block_text in blocks:
        if name in seen:
            raise TemplateError(f"template {template_id!r} repeats block {name!r}")
        seen.add(name)
        block_specs[name] = BlockSpec(name=name, text=block_text)

    for name in _csv(fields.get("frozen_blocks", "")):
        if name not in block_specs:
            raise TemplateError(
                f"template {template_id!r} freezes unknown block {name!r}"
            )
        block_specs[name] = BlockSpec(
            name=name, text=block_specs[name].text, editable=False,
        )

    assemblies = _parse_assemblies(fields.get("assemblies", ""), block_specs,
                                   template_id)

    return TemplateSpec(
        id=template_id,
        blocks=block_specs,
        assemblies=assemblies,
        params=_parse_params(fields.get("params", "")),
        requires=frozenset(_csv(fields.get("requires", ""))),
        outputs=tuple(_csv(fields.get("outputs", ""))),
        doc=fields.get("doc", ""),
        schema_version=int(fields.get("schema_version", "1")),
        source=source,
    )


# ── the builtin registry ─────────────────────────────────────────────────────

@lru_cache(maxsize=1)
def _builtin_registry() -> dict[str, TemplateSpec]:
    registry: dict[str, TemplateSpec] = {}
    for entry in resources.files(_BUILTIN_PACKAGE).iterdir():
        if not entry.name.endswith(".tmpl"):
            continue
        spec = parse_template(entry.read_text(encoding="utf-8"), source=entry.name)
        if spec.id in registry:
            raise TemplateError(
                f"two builtin templates claim id {spec.id!r}: "
                f"{registry[spec.id].source} and {entry.name}"
            )
        registry[spec.id] = spec
    return registry


def builtin_ids() -> list[str]:
    """Every registered builtin template id, sorted."""
    return sorted(_builtin_registry())


def builtin_spec(template_id: str) -> TemplateSpec:
    """The shipped spec for *template_id*. Never consults an override path."""
    try:
        return _builtin_registry()[template_id]
    except KeyError:
        raise TemplateError(
            f"unknown template {template_id!r} (known: {builtin_ids()})"
        ) from None


def builtin_assemble(template_id: str, block_names: Iterable[str]) -> str:
    """The shipped text for *template_id*, assembled from *block_names*."""
    return builtin_spec(template_id).assemble(block_names)


def builtin_text(template_id: str) -> str:
    """The shipped text of a single-block template, in file order."""
    spec = builtin_spec(template_id)
    return spec.assemble(spec.blocks)


# ── user overrides ───────────────────────────────────────────────────────────

#: Colon-separated search path replacing the user scope. Set to "" to disable
#: overrides entirely, which is what ``tests/conftest.py`` does: a developer's
#: own customisations must never change what CI-equivalent runs assert.
TEMPLATE_PATH_ENV = "XENIUM_VIEWER_TEMPLATE_PATH"

#: Set by ``--no-user-templates``. A process-wide off switch is the first thing
#: to reach for when a result is in doubt, so it must not require editing files.
_overrides_disabled = False


def set_overrides_enabled(enabled: bool) -> None:
    """Turn user overrides on or off for this process."""
    global _overrides_disabled
    _overrides_disabled = not enabled
    resolve.cache_clear()


def user_template_dir() -> Path:
    """Where a user's own templates are *written*. Not created until first save.

    Derived from the same search path used to read them, so a write can never
    land somewhere the reader does not look. Pointing
    ``XENIUM_VIEWER_TEMPLATE_PATH`` at a directory therefore redirects saving
    too — which is what makes the tests isolable with one environment variable
    instead of by patching, and what stopped them writing into the developer's
    real config directory.
    """
    path = search_path()
    if path:
        return path[0]
    return Path(platformdirs.user_config_dir("xenium-viewer")) / "templates"


def search_path() -> list[Path]:
    """Directories searched for overrides, lowest precedence first.

    Builtin is not on this list — it is not a directory on disk from the
    caller's point of view, and it can never be absent.
    """
    if _overrides_disabled:
        return []
    raw = os.environ.get(TEMPLATE_PATH_ENV)
    if raw is not None:
        return [Path(p).expanduser() for p in raw.split(os.pathsep) if p.strip()]
    return [user_template_dir()]


def _override_files(template_id: str) -> list[Path]:
    """Every override file for *template_id*, lowest precedence first."""
    found = []
    for directory in search_path():
        candidate = directory / f"{template_id}.tmpl"
        try:
            # `exists`, not `is_file`: a directory (or anything else) sitting
            # where a template should be is a mistake the user needs told about.
            # Skipping it silently would leave them believing an override is
            # active while the shipped template runs.
            if candidate.exists():
                found.append(candidate)
        except OSError:      # unreadable directory must not break startup
            continue
    return found


def _merge(base: TemplateSpec, override: TemplateSpec,
           source: str) -> TemplateSpec:
    """Base spec with *override*'s blocks substituted in, block by block.

    Per-block rather than whole-file, and that is the load-bearing choice. Most
    of a fork is text the user never touched; resolving per block means those
    parts keep tracking upstream automatically, so a release that fixes a
    template still reaches everyone who customised a *different* part of it.
    Whole-file override would freeze the entire template at the version it was
    forked from — which is how a user ends up quietly missing a correctness fix.

    Only the header fields that describe *text* are taken from the override.
    The contract (params, requires, outputs, assemblies) stays the builtin's:
    it is what the call site and the viewer agree on, not something a template
    file gets to redefine.
    """
    blocks = dict(base.blocks)
    changed = []
    for name, block in override.blocks.items():
        if name not in blocks:
            blocks[name] = block            # reported by validation, not here
            changed.append(name)
            continue
        if block.text != blocks[name].text:
            changed.append(name)
        blocks[name] = BlockSpec(
            name=name, text=block.text, editable=blocks[name].editable,
        )
    origin = TEMPLATE_USER if len(changed) == len(base.blocks) else TEMPLATE_BLENDED
    return replace(
        base, blocks=blocks, origin=origin if changed else base.origin,
        source=source if changed else base.source,
        # The override's schema-version, not the builtin's: it says which
        # contract the *user's* text was written against, which is the whole
        # question validation has to answer. Carrying the builtin's forward
        # would make a mismatch impossible to detect.
        schema_version=override.schema_version,
    )


@lru_cache(maxsize=None)
def resolve(template_id: str) -> ResolvedTemplate:
    """The spec actually used for *template_id*, plus why.

    Never raises on a bad override and never returns nothing: an unreadable,
    unparseable or invalid override degrades to the shipped template with the
    problems attached for the caller to surface. A broken file in a config
    directory must not make the viewer unlaunchable — the user would have no way
    in to fix it.
    """
    base = builtin_spec(template_id)
    files = _override_files(template_id)
    if not files:
        return ResolvedTemplate(spec=base, builtin=base)

    spec, problems, applied = base, [], None
    for path in files:
        try:
            candidate = parse_template(path.read_text(encoding="utf-8"),
                                       source=str(path))
        except (OSError, TemplateError, ValueError) as exc:
            problems.append(Problem(f"could not read {path}: {exc}"))
            continue
        merged = _merge(spec, candidate, source=str(path))
        found = validate(merged, builtin=base,
                         available=EXECUTOR_BASE_NAMES | base.requires)
        if any(p.severity == ERROR for p in found):
            problems += found
            continue
        problems += [p for p in found if p.severity == WARNING]
        spec, applied = merged, path

    if applied is None:
        # Say so. A user who edited a template and is silently getting the
        # shipped one will attribute the numbers to their own method.
        _report_rejection(template_id, problems)
        return ResolvedTemplate(spec=base, builtin=base, problems=tuple(problems),
                                rejected=True)
    return ResolvedTemplate(spec=spec, builtin=base, problems=tuple(problems),
                            path=applied)


def _report_rejection(template_id: str, problems) -> None:
    """Route to the reporter, but never let reporting break resolution."""
    try:
        from xenium_viewer.utils.reporting import report_template_rejected
        report_template_rejected(template_id, problems)
    except Exception:  # pragma: no cover - reporting is best-effort here
        pass


def resolved_text(template_id: str, block_names: Iterable[str]) -> str:
    """The text that will actually run for *template_id*, overrides included."""
    return resolve(template_id).spec.assemble(block_names)


def step_template(template_id: str, block_names: Iterable[str]) -> dict:
    """Step kwargs for *template_id*: the resolved text plus its provenance stamp.

    Returned as one dict, to be splatted straight into the constructor::

        Step(id=..., **step_template(TEMPLATE_ID, blocks), params=...)

    Text and stamp travel together on purpose. Fetching them separately would
    let a stamp describe a different resolution than the one that produced the
    string — the same class of drift, one level up, that ``Step`` exists to make
    impossible for code and parameters.
    """
    names = list(block_names)
    resolved = resolve(template_id)
    return {
        "template": resolved.spec.assemble(names),
        "template_id": template_id,
        "template_origin": resolved.origin,
        "template_hash": resolved.spec.hash_of(names),
    }


def clear_cache() -> None:
    """Forget resolved overrides — call after saving or reverting one."""
    resolve.cache_clear()
