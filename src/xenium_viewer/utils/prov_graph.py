"""
Provenance graph for reproducible-code recording.

The viewer records each analysis step as a node in a directed acyclic graph
(DAG) of *artifacts* rather than as an append-only list of code strings. Each
node knows the code that regenerates its artifact and which other nodes it
depends on. The exported notebook / script is *derived* from the graph by a
topological sort, so:

- Re-running a step (same ``id``) revises that node in place and flags its
  descendants stale, instead of being silently dropped or appended out of order.
- A missing dependency is an error at record time, not a ``NameError`` at replay.
- The emitted cell order always respects dependencies, regardless of the
  wall-clock order actions were taken (even across sessions).

This module is deliberately pure Python (no Qt/napari/nbformat imports) so the
graph logic can be unit-tested in isolation. nbformat wrapping lives in
``notebook_export.py``; this module only produces an ordered list of ``Cell``.

Usage:
    g = ProvGraph()
    g.upsert("preamble", "import scanpy as sc\\n...", kind="setup")
    g.upsert("normalize", "sc.pp.normalize_total(adata)", deps=["preamble"],
             kind="setup")
    g.upsert("clustering:leiden_r1.0", "sc.tl.leiden(adata, ...)",
             deps=["normalize"], label="Leiden clustering")
    cells = graph_to_cells(g)          # topo-ordered Cell list
    script = graph_to_script(g)        # flat .py text
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

# ── Node kinds ───────────────────────────────────────────────────────────────
SETUP = "setup"        # preamble, normalize — always sorts first
ARTIFACT = "artifact"  # reusable state: clustering, neighbors, DEG, nhood, ...
TERMINAL = "terminal"  # side-effect only: plot/export/viewer-only; no dependents

# Sort priority among nodes that become ready simultaneously in the topo sort.
_KIND_ORDER = {SETUP: 0, ARTIFACT: 1, TERMINAL: 2}


@dataclass
class ProvNode:
    """A single recorded step: the code, its dependencies, and metadata."""
    id: str
    code: str
    deps: list[str] = field(default_factory=list)
    kind: str = ARTIFACT
    label: Optional[str] = None          # markdown header for the derived cell
    params: dict = field(default_factory=dict)
    stale: bool = False                  # an upstream input changed after this ran
    seq: int = 0                         # insertion order, for deterministic sort


@dataclass
class Cell:
    """A derived notebook cell (nbformat-agnostic)."""
    cell_type: str                       # "code" | "markdown"
    source: str
    node_id: Optional[str] = None
    stale: bool = False


class CycleError(ValueError):
    """Raised when an edge would introduce a cycle in the DAG."""


class ProvGraph:
    """A mutable DAG of :class:`ProvNode`, keyed by node id."""

    def __init__(self) -> None:
        self._nodes: dict[str, ProvNode] = {}
        self._counter = 0

    # ── introspection ────────────────────────────────────────────────────────
    def __contains__(self, node_id: str) -> bool:
        return node_id in self._nodes

    def __len__(self) -> int:
        return len(self._nodes)

    def get(self, node_id: str) -> Optional[ProvNode]:
        return self._nodes.get(node_id)

    def nodes(self) -> list[ProvNode]:
        return list(self._nodes.values())

    # ── mutation ──────────────────────────────────────────────────────────────
    def upsert(
        self,
        node_id: str,
        code: str,
        deps: Iterable[str] = (),
        kind: str = ARTIFACT,
        label: Optional[str] = None,
        params: Optional[dict] = None,
    ) -> ProvNode:
        """Insert a node, or revise it in place if ``node_id`` already exists.

        - New id → inserted.
        - Existing id with identical code/deps/kind → no-op (label/params may
          refresh); descendants are *not* marked stale.
        - Existing id with changed code/deps → revised, and every transitive
          descendant is flagged ``stale`` (its input changed).

        Raises ``KeyError`` if a dependency id is unknown, and ``CycleError`` if
        the edge set would create a cycle.
        """
        deps = list(deps)
        if node_id in deps:
            raise CycleError(f"node '{node_id}' cannot depend on itself")
        for d in deps:
            if d not in self._nodes:
                raise KeyError(
                    f"node '{node_id}' depends on unknown node '{d}' "
                    f"(dependencies must be recorded first)"
                )
        # A cycle forms if any dependency can already reach node_id via deps edges.
        for d in deps:
            if self._reaches(d, node_id):
                raise CycleError(
                    f"adding '{node_id}' with dependency '{d}' creates a cycle"
                )

        existing = self._nodes.get(node_id)
        if existing is not None:
            unchanged = existing.code == code and existing.deps == deps \
                and existing.kind == kind
            if label is not None:
                existing.label = label
            if params is not None:
                existing.params = dict(params)
            if unchanged:
                return existing
            existing.code = code
            existing.deps = deps
            existing.kind = kind
            existing.stale = False  # the node itself was just re-recorded → fresh
            self._mark_descendants_stale(node_id)
            return existing

        self._counter += 1
        node = ProvNode(
            id=node_id, code=code, deps=deps, kind=kind,
            label=label, params=dict(params or {}), seq=self._counter,
        )
        self._nodes[node_id] = node
        return node

    def remove(self, node_id: str) -> None:
        """Remove a node. Raises if any other node still depends on it."""
        dependents = [n.id for n in self._nodes.values() if node_id in n.deps]
        if dependents:
            raise ValueError(
                f"cannot remove '{node_id}'; still required by {dependents}"
            )
        self._nodes.pop(node_id, None)

    # ── graph queries ─────────────────────────────────────────────────────────
    def _reaches(self, start: str, target: str) -> bool:
        """True if ``target`` is reachable from ``start`` following deps edges."""
        seen: set[str] = set()
        stack = [start]
        while stack:
            n = stack.pop()
            if n == target:
                return True
            if n in seen:
                continue
            seen.add(n)
            node = self._nodes.get(n)
            if node:
                stack.extend(node.deps)
        return False

    def _children_map(self) -> dict[str, list[str]]:
        """Reverse edges: node id → list of node ids that depend on it."""
        kids: dict[str, list[str]] = {nid: [] for nid in self._nodes}
        for n in self._nodes.values():
            for d in n.deps:
                if d in kids:
                    kids[d].append(n.id)
        return kids

    def _mark_descendants_stale(self, node_id: str) -> None:
        kids = self._children_map()
        seen: set[str] = set()
        stack = list(kids.get(node_id, []))
        while stack:
            c = stack.pop()
            if c in seen:
                continue
            seen.add(c)
            self._nodes[c].stale = True
            stack.extend(kids.get(c, []))

    def ancestors_closure(self, keep_ids: Iterable[str]) -> set[str]:
        """All nodes in ``keep_ids`` plus every transitive dependency of them."""
        closure: set[str] = set()
        stack = list(keep_ids)
        while stack:
            n = stack.pop()
            if n in closure or n not in self._nodes:
                continue
            closure.add(n)
            stack.extend(self._nodes[n].deps)
        return closure

    def topo_sort(self) -> list[str]:
        """Return node ids in dependency order.

        Ties (nodes ready simultaneously) break by (kind, id) — NOT insertion
        order — so setup sorts before artifacts before terminals and the derived
        notebook is identical regardless of the order actions were recorded
        (ordering invariance). ``id`` is unique, so the sort is total.
        """
        indeg = {n.id: len(n.deps) for n in self._nodes.values()}
        kids = self._children_map()

        def sort_key(nid: str):
            n = self._nodes[nid]
            return (_KIND_ORDER.get(n.kind, 1), n.id)

        ready = sorted((nid for nid, d in indeg.items() if d == 0), key=sort_key)
        order: list[str] = []
        while ready:
            nid = ready.pop(0)
            order.append(nid)
            newly = []
            for c in kids.get(nid, []):
                indeg[c] -= 1
                if indeg[c] == 0:
                    newly.append(c)
            if newly:
                ready.extend(newly)
                ready.sort(key=sort_key)
        if len(order) != len(self._nodes):
            raise CycleError("cycle detected during topological sort")
        return order

    # ── serialization ─────────────────────────────────────────────────────────
    def to_list(self) -> list[dict]:
        """Serialize to a JSON-friendly list (order-preserving by seq)."""
        return [
            {
                "id": n.id, "code": n.code, "deps": list(n.deps),
                "kind": n.kind, "label": n.label, "params": n.params,
                "stale": n.stale, "seq": n.seq,
            }
            for n in sorted(self._nodes.values(), key=lambda n: n.seq)
        ]

    @classmethod
    def from_list(cls, items: list[dict]) -> "ProvGraph":
        g = cls()
        max_seq = 0
        for it in items:
            node = ProvNode(
                id=it["id"], code=it["code"], deps=list(it.get("deps", [])),
                kind=it.get("kind", ARTIFACT), label=it.get("label"),
                params=dict(it.get("params", {})), stale=bool(it.get("stale", False)),
                seq=int(it.get("seq", 0)),
            )
            g._nodes[node.id] = node
            max_seq = max(max_seq, node.seq)
        g._counter = max_seq
        return g


# ── derivation ────────────────────────────────────────────────────────────────
def graph_to_cells(
    graph: ProvGraph,
    include_terminals: bool = True,
    keep_ids: Optional[Iterable[str]] = None,
) -> list[Cell]:
    """Derive an ordered list of :class:`Cell` from the graph.

    - ``keep_ids`` (optional): keep only these nodes and their transitive
      dependencies, dropping abandoned experiment branches.
    - ``include_terminals``: drop terminal (plot/export) nodes when False.
    """
    order = graph.topo_sort()
    include = set(order)
    if keep_ids is not None:
        include = graph.ancestors_closure(keep_ids)
    cells: list[Cell] = []
    for nid in order:
        if nid not in include:
            continue
        node = graph.get(nid)
        if node is None:
            continue
        if not include_terminals and node.kind == TERMINAL:
            continue
        if node.label:
            cells.append(Cell("markdown", f"## {node.label}", node_id=nid))
        cells.append(Cell("code", node.code, node_id=nid, stale=node.stale))
    return cells


def graph_to_script(graph: ProvGraph, include_terminals: bool = True) -> str:
    """Flat ``.py`` rendering: the code cells joined in topological order."""
    parts = [
        c.source for c in graph_to_cells(graph, include_terminals=include_terminals)
        if c.cell_type == "code"
    ]
    return "\n".join(parts) + "\n"
