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

import re
from functools import lru_cache
from importlib import resources
from typing import Iterable, Optional

from xenium_viewer.utils.step_templates.spec import (
    BLOCK_MARKER,
    BlockSpec,
    ParamSpec,
    TemplateSpec,
)

_BUILTIN_PACKAGE = "xenium_viewer.utils.step_templates.builtin"

#: ``name:type`` or ``name:type?`` — the trailing ``?`` means "optional".
_PARAM_RE = re.compile(r"^(?P<name>\w+)\s*:\s*(?P<type>\w+)(?P<optional>\?)?$")


class TemplateError(ValueError):
    """A ``.tmpl`` file is malformed or references something unknown."""


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
        if key is not None and line.strip():
            fields[key] = f"{fields[key]} {line.strip()}".strip()
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
