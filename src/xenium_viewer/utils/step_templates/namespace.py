"""The names a template may reach for without binding them itself.

``tabs._helpers._get_executor`` seeds the :class:`~xenium_viewer.utils.steps.StepExecutor`
namespace with these, and the exported notebook's preamble binds the same set —
that correspondence is what makes a recorded cell replayable in a clean kernel.

The list lives here, apart from the executor, for one reason: template
validation (``check_step``) has to know exactly which names are free-for-use,
and a validator working from its own copy of the list would eventually disagree
with the executor. So the executor builds its dict and then calls
:func:`check_base_namespace` on it — the two cannot drift without failing loudly
at startup.

This module imports nothing heavy on purpose. Naming the names must not cost a
scanpy import, or the validator could not run outside the GUI.
"""

from __future__ import annotations

from typing import Iterable

#: Names bound in a fresh executor namespace, before any step has run.
#:
#: Third-party modules (``sc``, ``sq``, ``sd``, ``pd``, ``np``, ``plt``) plus
#: ``Path`` are the imports the notebook preamble makes. ``sd`` is
#: ``spatialdata`` itself, which templates need for the operations that belong
#: to the container rather than to the table — ``sd.polygon_query`` and the
#: ``sd.transformations`` that let a template *declare* a coordinate frame
#: instead of applying it by hand. ``data_path``, ``sdata`` and ``adata`` are
#: the dataset the viewer already loaded.
#:
#: Names bound *by* a step (``adata_norm``, ``adata_leiden``, ``rank_df``,
#: ``roi_polygons``, …) are deliberately absent: a step that needs one must
#: declare the producing step in its ``deps``, and template validation checks it
#: against ``EXECUTOR_BASE_NAMES | spec.requires`` rather than against whatever
#: happens to be lying around in a long-lived session.
EXECUTOR_BASE_NAMES = frozenset({
    "sc", "sq", "sd", "pd", "np", "plt",
    "Path",
    "data_path", "sdata", "adata",
})


class NamespaceMismatch(RuntimeError):
    """The executor was seeded with a namespace that is not the declared one."""


def check_base_namespace(names: Iterable[str]) -> None:
    """Raise :class:`NamespaceMismatch` unless *names* is exactly the base set.

    Called on the freshly built executor namespace. Both directions matter: an
    extra name is one templates may silently come to rely on but the notebook
    never binds (it would replay as a ``NameError``), and a missing one is a
    name validation promises but execution does not provide.
    """
    actual = set(names)
    extra = actual - EXECUTOR_BASE_NAMES
    missing = EXECUTOR_BASE_NAMES - actual
    if extra or missing:
        raise NamespaceMismatch(
            "executor namespace does not match EXECUTOR_BASE_NAMES "
            f"(unexpected: {sorted(extra)}, missing: {sorted(missing)}). "
            "Update utils/step_templates/namespace.py and confirm the notebook "
            "preamble binds the same names."
        )
