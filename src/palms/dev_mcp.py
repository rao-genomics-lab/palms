"""Development-only MCP bridge: expose the running viewer to an AI assistant.

Off by default and never imported unless ``--mcp`` is passed. ``napari-mcp``
lives in the ``mcp`` extra rather than ``environment.yml`` because it is alpha
software and because the bridge is an **unauthenticated HTTP server that
executes arbitrary Python in the viewer process**. That is acceptable as a
localhost dev tool and unacceptable as a default.

What it buys: the GUI proper and the spatial-analysis tabs have no automated
coverage, so exercising them means a human clicking. The bridge lets an
assistant drive the *already loaded* dataset — full-window screenshots that
include the Xenium Controls and Plots docks, plus ``execute_code`` on the Qt
main thread with ``ctx`` in reach.

Two things to know before using it:

- **Bridge mode, not standalone.** ``napari_mcp``'s own ``init_viewer`` would
  create a second, empty viewer. ``NapariBridgeServer`` wraps the viewer we
  already built, which is the only one with a dataset in it.
- **A modal dialog deadlocks it.** Every bridge call marshals onto the Qt main
  thread, so while one of this app's modals is up (stale cache, corrupt cache,
  close-with-running-job) nothing gets serviced until a human clicks it. Long
  main-thread work — a cache build, a pyramid load — has the same effect and
  hits ``NapariBridgeServer``'s 300 s timeout.
"""
from __future__ import annotations

DEFAULT_PORT = 9999


def start_bridge(viewer, port: int = DEFAULT_PORT):
    """Start the napari-mcp bridge over ``viewer``. Returns the server, or None.

    Never raises: a dev aid that cannot start must not stop the viewer from
    opening. The import is deliberately inside the function so a normal launch
    never pays for it and never requires the extra to be installed.
    """
    try:
        from napari_mcp.bridge_server import NapariBridgeServer
    except ImportError:
        print(
            "  --mcp: napari-mcp is not installed. "
            "pip install 'palms[mcp]' to enable it."
        )
        return None

    try:
        server = NapariBridgeServer(viewer, port=port)
        server.start()
    except Exception as exc:
        print(f"  --mcp: bridge failed to start on port {port}: {exc}")
        return None

    print(f"  --mcp: MCP bridge listening on http://127.0.0.1:{port}/mcp")
    print("         Arbitrary code execution, no auth — localhost dev use only.")
    return server


def publish_context(viewer, ctx) -> None:
    """Make ``ctx`` reachable from the bridge's ``execute_code`` namespace.

    That namespace is seeded with ``viewer`` and nothing else, so without this
    an assistant can inspect layers but cannot reach the provenance graph, the
    clusterings dict or any tab's exports — i.e. everything worth checking.

    Called from ``_push_to_console`` because that is the one place a *fresh*
    ``ctx`` is published to interactive surfaces, and ``reload_dataset`` builds
    a new one; hanging it anywhere else would leave a stale object behind after
    a dataset switch, which is worse than having none.
    """
    try:
        viewer._xv_ctx = ctx
    except Exception:
        pass
