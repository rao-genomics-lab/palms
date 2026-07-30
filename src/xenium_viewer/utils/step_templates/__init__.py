"""Declarative analysis templates: the source the GUI runs and the notebook records.

This package owns the *text* of every analysis step. The call sites in
``xenium_viewer.tabs`` own which blocks are selected and what the parameters
are; see :mod:`xenium_viewer.utils.steps` for how a template becomes a single
string handed to both ``exec`` and the provenance graph.

Kept free of Qt, napari and scanpy imports so that template validation can run
in isolation — the same discipline as ``utils.steps`` and ``utils.prov_graph``.
"""

from xenium_viewer.utils.step_templates.loader import (
    TEMPLATE_PATH_ENV,
    ResolvedTemplate,
    TemplateError,
    builtin_assemble,
    builtin_ids,
    builtin_spec,
    builtin_text,
    clear_cache,
    parse_template,
    resolve,
    resolved_text,
    step_template,
    search_path,
    set_overrides_enabled,
    user_template_dir,
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
from xenium_viewer.utils.step_templates.validate import (
    ERROR,
    WARNING,
    Problem,
    placeholders,
    validate,
)

__all__ = [
    # namespace
    "EXECUTOR_BASE_NAMES",
    "NamespaceMismatch",
    "check_base_namespace",
    # spec
    "BlockSpec",
    "ParamSpec",
    "TemplateSpec",
    # loader — builtin scope only, never sees an override
    "TemplateError",
    "builtin_assemble",
    "builtin_ids",
    "builtin_spec",
    "builtin_text",
    "parse_template",
    # loader — full resolution, overrides included
    "TEMPLATE_PATH_ENV",
    "ResolvedTemplate",
    "clear_cache",
    "resolve",
    "resolved_text",
    "step_template",
    "search_path",
    "set_overrides_enabled",
    "user_template_dir",
    # validation
    "ERROR",
    "WARNING",
    "Problem",
    "placeholders",
    "validate",
]
