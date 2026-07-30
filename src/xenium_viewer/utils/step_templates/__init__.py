"""Declarative analysis templates: the source the GUI runs and the notebook records.

This package owns the *text* of every analysis step. The call sites in
``xenium_viewer.tabs`` own which blocks are selected and what the parameters
are; see :mod:`xenium_viewer.utils.steps` for how a template becomes a single
string handed to both ``exec`` and the provenance graph.

Kept free of Qt, napari and scanpy imports so that template validation can run
in isolation — the same discipline as ``utils.steps`` and ``utils.prov_graph``.
"""

from xenium_viewer.utils.step_templates.loader import (
    TemplateError,
    builtin_spec,
    builtin_assemble,
    builtin_ids,
    builtin_text,
    parse_template,
)
from xenium_viewer.utils.step_templates.namespace import (
    EXECUTOR_BASE_NAMES,
    NamespaceMismatch,
    check_base_namespace,
)
from xenium_viewer.utils.step_templates.spec import (
    BlockSpec,
    ParamSpec,
    TemplateSpec,
)

__all__ = [
    "EXECUTOR_BASE_NAMES",
    "NamespaceMismatch",
    "check_base_namespace",
    "BlockSpec",
    "ParamSpec",
    "TemplateSpec",
    "TemplateError",
    "builtin_spec",
    "builtin_assemble",
    "builtin_ids",
    "builtin_text",
    "parse_template",
]
