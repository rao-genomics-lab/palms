"""What the analysis was run *with* — captured as a recorded cell.

A notebook that reproduces a result only reproduces it against the same
software. The recorded code said which functions were called but never which
versions answered the call, so a replay that disagreed gave no way to tell a
real difference from a scanpy upgrade.

This module supplies the ``environment`` node: a comment block pinning the
versions present when the analysis was recorded, plus two executable lines that
seed the global RNGs and print the versions the *replay* is actually running
with. Reading the two side by side is the whole point — the cell is deliberately
not an assertion, because a version mismatch is information, not a failure.

Pure Python (no Qt/napari), so it can be tested and called from scripts.
"""

from __future__ import annotations

import sys
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version

# Distributions worth pinning: the ones whose behaviour a recorded step can
# actually depend on. Distribution names (``scikit-learn``), not import names.
RECORDED_PACKAGES = (
    "scanpy",
    "anndata",
    "squidpy",
    "spatialdata",
    "spatialdata-io",
    "zarr",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "leidenalg",
    "igraph",
    "infercnvpy",
    "xenium-viewer",
)

# Seeds the recorded cell sets. Steps that take an explicit ``random_state``
# still carry it in their own source; this covers everything that does not.
DEFAULT_SEED = 0


def package_versions(names=RECORDED_PACKAGES) -> dict:
    """Map distribution name → version string, ``None`` when not installed."""
    out = {"python": sys.version.split()[0]}
    for name in names:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = None
    return out


def pins(code: str) -> list[str]:
    """The version-pin comment lines of an ``environment`` cell."""
    return [line for line in code.splitlines() if line.startswith("#   ")]


def same_environment(a: str, b: str) -> bool:
    """True when two ``environment`` cells pin the same versions.

    Re-opening a dataset re-records the node. Comparing the *pins* rather than
    the whole cell keeps a session that changed nothing from rewriting the
    timestamp — which would churn the exported notebook and the sidecar on
    every launch for no information gained.
    """
    return pins(a) == pins(b)


def environment_code(versions: dict | None = None,
                     timestamp: str | None = None,
                     seed: int = DEFAULT_SEED) -> str:
    """Return the source of the ``environment`` node.

    *versions* defaults to the live interpreter's; *timestamp* to now. Both are
    injectable so the cell is testable and so a re-record does not churn.
    """
    if versions is None:
        versions = package_versions()
    if timestamp is None:
        timestamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")

    width = max((len(k) for k in versions), default=0)
    pins = "\n".join(
        f"#   {name:<{width}}  {ver}"
        for name, ver in versions.items() if ver is not None
    )
    missing = [name for name, ver in versions.items() if ver is None]
    if missing:
        pins += f"\n#   (not installed: {', '.join(missing)})"

    return (
        f"# Recorded {timestamp} with:\n"
        f"{pins}\n"
        "#\n"
        "# The versions below are this run's. Compare them with the block above\n"
        "# before treating any difference in the results as a real one.\n"
        "import random\n"
        "import numpy as np\n"
        "import scanpy as sc\n"
        "\n"
        f"random.seed({seed})\n"
        f"np.random.seed({seed})\n"
        "sc.logging.print_header()"
    )
