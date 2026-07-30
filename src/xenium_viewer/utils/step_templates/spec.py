"""What a template *is*, declaratively.

A template used to be a private module constant in a tab module — sometimes a
flat string, sometimes seven fragments concatenated by a private function keyed
on booleans. That works when the only reader is the call site next to it, and
stops working the moment anything else needs to know what a template takes, what
it binds, or which combinations of fragments are legal.

The model here says exactly that, and nothing more:

* A template is an **ordered dict of named blocks**.
* The **call site chooses which blocks**; the registry owns their *text*.
* An **assembly** is one legal block sequence. They are enumerated, not derived,
  because the choosing logic is real Python (``if use_hvg or do_scale``) and
  turning it into a mini-language would buy nothing.

That split matters beyond tidiness: the branch structure *is* what the widgets
mean, while the statements are what someone actually wants to change. Keeping
selection in Python and text in files lets the second be edited without letting
the first be broken.

Deliberately not modelled here: ``deps``. Dependencies are computed from runtime
state at the call site (which clustering is selected, whether a filter is on), so
they stay in Python. Making graph edges data would put the DAG's integrity —
which ``ProvGraph.upsert`` currently guarantees — into an editable file.

Pure Python: no Qt, no napari, no scanpy.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Optional

#: Marker introducing a block in a ``.tmpl`` file. A comment, so the whole file
#: stays valid Python and the editor's syntax highlighter needs no special case.
BLOCK_MARKER = "#--- block "


@dataclass(frozen=True)
class ParamSpec:
    """One ``$name`` a template expects, and what kind of literal it must be."""

    name: str
    type: str = "str"          # str | int | float | bool | list | dict | tuple
    required: bool = True
    doc: str = ""


@dataclass(frozen=True)
class BlockSpec:
    """One named, individually addressable piece of a template.

    ``editable=False`` marks text that must not be customised even once user
    overrides exist — a version workaround whose byte identity with a helper
    elsewhere is pinned by a test, where no validation gate could tell the user
    "changing this breaks on the next pandas".
    """

    name: str
    text: str
    editable: bool = True


@dataclass(frozen=True)
class TemplateSpec:
    """A registered template: its blocks, its contract, and where it came from."""

    id: str
    blocks: dict[str, BlockSpec]
    assemblies: tuple[tuple[str, ...], ...] = ()
    params: tuple[ParamSpec, ...] = ()
    requires: frozenset[str] = frozenset()
    outputs: tuple[str, ...] = ()
    doc: str = ""
    schema_version: int = 1
    sample_params: dict = field(default_factory=dict)
    origin: str = "builtin"
    source: Optional[str] = None

    # ── assembly ─────────────────────────────────────────────────────────────
    def assemble(self, block_names) -> str:
        """Concatenate the named blocks, in the order given.

        Every block's text begins with a newline, so plain concatenation
        reproduces the original constants byte for byte.
        """
        names = list(block_names)
        missing = [n for n in names if n not in self.blocks]
        if missing:
            raise KeyError(
                f"template {self.id!r} has no block(s) {missing} "
                f"(known: {sorted(self.blocks)})"
            )
        return "".join(self.blocks[n].text for n in names)

    def hash_of(self, block_names) -> str:
        """sha256 of the assembled *template* text, before substitution.

        Hashing the rendered code would add nothing — the provenance node
        already stores that verbatim. Hashing the template is what lets a reader
        say "this is the shipped ``clustering.leiden``".
        """
        return hashlib.sha256(self.assemble(block_names).encode()).hexdigest()

    # ── introspection ────────────────────────────────────────────────────────
    @property
    def param_names(self) -> set[str]:
        return {p.name for p in self.params}

    @property
    def required_params(self) -> set[str]:
        return {p.name for p in self.params if p.required}

    def param(self, name: str) -> Optional[ParamSpec]:
        for p in self.params:
            if p.name == name:
                return p
        return None

    def synth_params(self) -> dict:
        """A literal of the right type for each declared param.

        Enough to *render* a template and parse the result, which is what
        static validation needs — nothing here is meant to be executed, so the
        values only have to be well-typed, not meaningful. Overridden by
        ``sample_params`` where a realistic value reads better in a preview.
        """
        values = {name: _SYNTH.get(spec.type, "sample")
                  for name, spec in ((p.name, p) for p in self.params)}
        values.update(self.sample_params)
        return values


#: One literal per declared param type, for render-and-parse validation.
_SYNTH = {
    "str": "sample",
    "int": 1,
    "float": 1.0,
    "bool": True,
    "list": ["sample"],
    "tuple": ("sample",),
    "dict": {"group": ["sample"]},
}
