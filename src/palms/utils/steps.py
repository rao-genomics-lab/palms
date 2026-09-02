"""
Analysis steps: the code the GUI *executes* is the code the notebook *records*.

Before this module, each tab callback ran an analysis and, separately, wrote a
string describing that analysis into the provenance graph. Two expressions of
the same computation, maintained by hand — which drifted (the GUI ran
``leidenalg.find_partition(..., seed=42)`` while the notebook recorded
``sc.tl.leiden(..., random_state=0)``; the GUI normalised with
``target_sum=1e4`` while the notebook recorded scanpy's median default).

Here there is only one expression. A :class:`Step` holds a *template* of plain
scverse source plus a dict of literal ``params``. Rendering produces a single
string, and that same string object is both handed to ``exec`` and stored on
the provenance node. Drift is not unlikely — it is impossible.

The invariant that makes this hold, and which review must enforce:

    A tab callback may never call an analysis function with a widget value.
    It may only build a ``params`` dict.

Templates use :class:`string.Template` (``$name``) rather than ``str.format``
so that ``{...}`` dict literals and f-strings inside templates are left alone.
Params are substituted via ``repr()``, so what lands in the source is a Python
literal that ``ast.literal_eval`` round-trips exactly.

This module is deliberately pure Python (no Qt/napari/scanpy imports) so it can
be unit-tested in isolation, like ``prov_graph``.

Usage:
    ex = StepExecutor(namespace={"adata": adata}, graph=graph)
    ex.run(Step(
        id="clustering:leiden_r1.0",
        template='sc.tl.leiden(adata, resolution=$resolution, key_added=$key)',
        params={"resolution": 1.0, "key": "leiden_r1.0"},
        deps=["normalize"],
        label="Leiden clustering",
    ))
"""

from __future__ import annotations

import ast
import builtins
import math
import threading
from dataclasses import dataclass, field
from string import Template
from typing import Any, Callable, Iterable, Optional

from palms.utils.prov_graph import ARTIFACT, TEMPLATE_BUILTIN, ProvGraph

# Params must render to a literal that ``ast.literal_eval`` accepts. Anything
# else (a numpy scalar, an AnnData, an ndarray) has to live in the namespace
# under a name that the template references instead.
_LITERAL_TYPES = (str, int, float, bool, type(None), list, tuple, dict, set)


class StepError(RuntimeError):
    """Raised when a step fails to validate, render, or execute."""


class ParamError(StepError):
    """Raised when a param cannot be rendered as a faithful Python literal."""


@dataclass
class Step:
    """One analysis action: the source that performs it and its parameters.

    Attributes
    ----------
    id
        Provenance node id — the artifact this step produces, e.g.
        ``"clustering:leiden_r1.0"``. Re-running with the same id revises the
        node in place (see :meth:`ProvGraph.upsert`).
    template
        Plain scverse source with ``$name`` placeholders. Must reference only
        (a) names declared in ``params`` and (b) names already bound in the
        executor's namespace.
    params
        Literal-rendering values. Validated by :func:`validate_params`.
    outputs
        Names to lift back out of the namespace after execution and hand to the
        caller. The GUI reads results from here rather than keeping its own copy.
    """

    id: str
    template: str
    params: dict = field(default_factory=dict)
    deps: list[str] = field(default_factory=list)
    kind: str = ARTIFACT
    label: Optional[str] = None
    outputs: list[str] = field(default_factory=list)
    # Which registered template ``template`` came from, and whether it is the
    # shipped text. Carried onto the provenance node so a reader can tell a
    # stock run from a customised one; see :mod:`palms.utils.prov_graph`.
    template_id: Optional[str] = None
    template_origin: str = TEMPLATE_BUILTIN
    template_hash: Optional[str] = None

    def render(self) -> str:
        """Return the single source string — executed *and* recorded."""
        return render(self.template, self.params)


# ── parameter handling ───────────────────────────────────────────────────────

def coerce(value: Any) -> Any:
    """Best-effort conversion of a widget/numpy value to a plain Python literal.

    Call this at the *widget boundary*, when building ``params`` — never at
    render time, which would reintroduce two different values.

    numpy scalars and 0-d arrays expose ``.item()``; numpy arrays and pandas
    Index/Series become lists. Containers are coerced recursively.
    """
    if isinstance(value, (str, bytes)) or value is None:
        return value
    if isinstance(value, bool):  # before int — bool is an int subclass
        return value
    if isinstance(value, dict):
        return {coerce(k): coerce(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return type(value)(coerce(v) for v in value)
    # numpy scalars / 0-d arrays, pandas scalars
    item = getattr(value, "item", None)
    if callable(item):
        try:
            if getattr(value, "ndim", 0) == 0:
                return coerce(item())
        except (ValueError, TypeError):
            pass
    # ndarray / Index / Series → list
    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        try:
            return coerce(tolist())
        except (ValueError, TypeError):
            pass
    return value


def validate_param(name: str, value: Any) -> None:
    """Raise :class:`ParamError` unless *value* renders as a faithful literal.

    The check is ``ast.literal_eval(repr(value)) == value``. That rejects numpy
    scalars (whose repr under NumPy 2 is ``np.float64(1.0)``, which needs numpy
    imported and is not stable across versions), objects with a default
    ``<... at 0x...>`` repr, and non-finite floats.
    """
    if isinstance(value, float) and not math.isfinite(value):
        raise ParamError(
            f"param {name!r} is {value!r}; non-finite floats cannot be a "
            f"reproducible parameter"
        )
    if not isinstance(value, _LITERAL_TYPES):
        raise ParamError(
            f"param {name!r} has type {type(value).__name__}, which has no "
            f"literal form. Coerce it at the widget boundary (utils.steps.coerce), "
            f"or bind it in the namespace and reference it by name in the template."
        )
    text = repr(value)
    try:
        parsed = ast.literal_eval(text)
    except (ValueError, SyntaxError) as exc:
        raise ParamError(
            f"param {name!r} does not round-trip through repr(): {text!r} ({exc})"
        ) from exc
    if parsed != value or type(parsed) is not type(value):
        raise ParamError(
            f"param {name!r} does not round-trip through repr(): "
            f"{value!r} -> {text!r} -> {parsed!r}"
        )


def validate_params(params: dict) -> None:
    """Validate every param. See :func:`validate_param`."""
    for name, value in params.items():
        validate_param(name, value)


def render(template: str, params: dict) -> str:
    """Substitute ``$name`` placeholders with ``repr()`` of each param.

    Raises :class:`StepError` if the template references an undeclared name —
    failing loudly at record time beats a ``NameError`` at replay.
    """
    validate_params(params)
    try:
        return Template(template).substitute(
            {name: repr(value) for name, value in params.items()}
        )
    except KeyError as exc:
        raise StepError(
            f"template references ${exc.args[0]}, which is not in params "
            f"(declared: {sorted(params)})"
        ) from exc
    except ValueError as exc:  # bad placeholder syntax, e.g. a bare '$'
        raise StepError(f"malformed template: {exc}") from exc


# ── template hygiene (the CI lint that makes the guarantee auditable) ────────

def free_names(code: str) -> set[str]:
    """Names *loaded* by ``code`` that it never binds and that aren't builtins.

    A conservative, module-level heuristic used to check that a rendered
    template only reaches for names the namespace guarantees. Any name bound
    anywhere in the module counts as bound, which avoids false positives at the
    cost of missing use-before-assignment — acceptable for a lint.
    """
    tree = ast.parse(code)
    bound: set[str] = set()
    loaded: list[str] = []

    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            if isinstance(node.ctx, ast.Load):
                loaded.append(node.id)
            else:  # Store | Del
                bound.add(node.id)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                bound.add((alias.asname or alias.name).split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                bound.add(alias.asname or alias.name)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bound.add(node.name)
        elif isinstance(node, ast.ExceptHandler) and node.name:
            bound.add(node.name)
        elif isinstance(node, ast.arg):
            bound.add(node.arg)
        elif isinstance(node, ast.Global) or isinstance(node, ast.Nonlocal):
            bound.update(node.names)

    builtin_names = set(dir(builtins))
    return {name for name in loaded if name not in bound and name not in builtin_names}


def check_step(step: Step, available: Iterable[str]) -> set[str]:
    """Return the names *step* needs that neither it nor the namespace provides.

    An empty set means the step is self-contained given *available*. Used by the
    CI template lint; an empty result for every registered step is what makes
    the exactness guarantee auditable rather than merely asserted.
    """
    return free_names(step.render()) - set(available)


# ── execution ────────────────────────────────────────────────────────────────

ProgressFn = Callable[[int, int, str], None]


class StepExecutor:
    """Runs :class:`Step` objects against one shared namespace and records them.

    The namespace mirrors the exported notebook's globals: ``data_path``,
    ``sdata``, ``adata`` and result objects all live here, so the GUI and the
    notebook operate on the same state rather than on parallel copies.

    Execution is serialised — steps mutate shared state, so two may never run
    concurrently even from different napari worker threads.
    """

    def __init__(self, namespace: Optional[dict] = None,
                 graph: Optional[ProvGraph] = None) -> None:
        self.ns: dict = namespace if namespace is not None else {}
        self.graph: ProvGraph = graph if graph is not None else ProvGraph()
        self._lock = threading.RLock()

    # -- introspection -------------------------------------------------------
    def names(self) -> set[str]:
        """Names currently bound in the namespace."""
        return set(self.ns)

    def get(self, name: str, default: Any = None) -> Any:
        return self.ns.get(name, default)

    # -- the one path from a Step to executed + recorded code ----------------
    def run(self, step: Step, progress: Optional[ProgressFn] = None) -> dict:
        """Render *step*, execute it, record it, and return its declared outputs.

        The rendered string is executed one top-level statement at a time so a
        long step can report progress, but the source compiled is byte-identical
        to the source recorded — only the granularity of the ``exec`` calls
        differs, and statement line numbers are preserved so tracebacks point
        into the recorded cell.

        Recording happens only on success: a step that raises leaves no node
        behind claiming to have produced an artifact that does not exist.
        """
        code = step.render()
        filename = f"<step:{step.id}>"

        with self._lock:
            _exec_statements(code, self.ns, filename, f"step {step.id!r}",
                             progress=progress)

            self.graph.upsert(
                step.id, code, deps=list(step.deps), kind=step.kind,
                label=step.label, params=dict(step.params),
                template_id=step.template_id,
                template_origin=step.template_origin,
                template_hash=step.template_hash,
            )

            missing = [name for name in step.outputs if name not in self.ns]
            if missing:
                raise StepError(
                    f"step {step.id!r} declared outputs {missing} that it did not bind"
                )
            return {name: self.ns[name] for name in step.outputs}

    # -- display only: executes a Step's source, records nothing ---------------
    def preview(self, step: Step, bindings: Optional[dict] = None) -> dict:
        """Render and execute *step* for display, and record nothing.

        This exists so display code can reuse a template's text instead of
        keeping a second implementation of the same computation — the drift
        risk that made CopyKAT's reconstruction cell a documented exception
        rather than a pattern. It is **not** a second way to run analysis, and
        recording is deliberately not a flag on :meth:`run`: recording is not a
        mode of ``run``, it is what ``run`` is.

        Two properties, both load-bearing:

        1. The provenance graph is never reached from here. There is no
           ``upsert`` in this method and nothing that turns one on.
        2. Execution goes into a **copy** of the namespace, not the namespace.
           That is a correctness requirement rather than tidiness: a preview
           that bound ``transcript_points`` in the shared namespace would let a
           later recorded step consume those values while recording source that
           says it read them from somewhere else — exactly the executed-versus-
           recorded drift this module exists to make impossible. The copy is
           shallow, so ``adata`` and ``sdata`` are the same objects and a
           template that *mutates* them would mutate them for real; that is why
           only read-only templates should be previewed.

        *bindings* supplies names the template needs that the shared namespace
        does not hold — the point of the whole exercise, since the preview's
        input is fetched by a different route than the recorded step's. It may
        not shadow a name already bound, or the preview would silently describe
        different data than the one it is standing in for.

        Note the absence of a ``progress`` parameter. If a preview is slow
        enough to need progress reporting it is not a preview, and leaving the
        parameter off keeps that pressure visible.

        Guarded by ``tests/test_preview_never_records.py``.
        """
        scratch = dict(self.ns)
        shadowed = sorted(set(bindings or {}) & set(self.ns))
        if shadowed:
            raise StepError(
                f"preview {step.id!r} would shadow namespace name(s) {shadowed}; "
                f"a preview may only bind names the shared namespace does not hold"
            )
        scratch.update(bindings or {})

        code = step.render()
        with self._lock:
            _exec_statements(code, scratch, f"<preview:{step.id}>",
                             f"preview {step.id!r}")

        missing = [name for name in step.outputs if name not in scratch]
        if missing:
            raise StepError(
                f"preview {step.id!r} declared outputs {missing} that it did not bind"
            )
        return {name: scratch[name] for name in step.outputs}


def _exec_statements(code: str, namespace: dict, filename: str, what: str,
                     progress: Optional[ProgressFn] = None) -> None:
    """Execute *code* one top-level statement at a time into *namespace*.

    The one exec loop, shared by :meth:`StepExecutor.run` and
    :meth:`StepExecutor.preview` so the two cannot drift. Statement-by-statement
    only so a long step can report progress: the source compiled is
    byte-identical to the source passed in, and line numbers are preserved so a
    traceback points into the cell.
    """
    try:
        tree = ast.parse(code, filename=filename)
    except SyntaxError as exc:
        raise StepError(f"{what} is not valid Python: {exc}") from exc

    statements = tree.body
    total = len(statements)
    for index, statement in enumerate(statements, start=1):
        if progress is not None:
            progress(index, total, _statement_label(code, statement))
        module = ast.Module(body=[statement], type_ignores=[])
        try:
            exec(compile(module, filename, "exec"), namespace)  # noqa: S102
        except Exception as exc:
            raise StepError(
                f"{what} failed at statement {index}/{total}: "
                f"{type(exc).__name__}: {exc}"
            ) from exc


def _statement_label(code: str, statement: ast.stmt) -> str:
    """A short human-readable label for a single statement, for the status bar."""
    segment = ast.get_source_segment(code, statement)
    if not segment:
        return type(statement).__name__
    first_line = segment.strip().splitlines()[0].strip()
    return first_line if len(first_line) <= 80 else first_line[:77] + "..."
