"""Install the CopyKAT R package (GitHub-only) into the active env's R.

CopyKAT's *dependencies* come from conda (see environment.yml), but copykat
itself is only on GitHub, so it can't live in environment.yml. This module
installs it via ``remotes::install_github`` through rpy2 — idempotently, so it's
safe to call on every CopyKAT run. Called lazily from the CopyKAT pipeline
branch; ``python -m palms.install_copykat`` runs it standalone. There is no
console script on purpose: this needs rpy2 and R, which only the
``palms_copykat`` env has, so an entry point installed by the main env would
be a command that can never work.
"""

from __future__ import annotations

import sys

COPYKAT_REPO = "navinlabcode/copykat"


def _copykat_available() -> bool:
    import rpy2.robjects as ro
    return bool(ro.r('requireNamespace("copykat", quietly=TRUE)')[0])


def ensure_copykat_installed(quiet: bool = False) -> bool:
    """Ensure the copykat R package is importable; install it from GitHub if not.

    Returns True if copykat is available afterwards. Raises RuntimeError with an
    actionable message if rpy2/R is missing or the install fails.
    """
    try:
        import rpy2.robjects as ro  # noqa: F401
    except ImportError as e:
        raise RuntimeError(
            "CopyKAT needs rpy2 + R. Install the env from environment.yml "
            "(conda/mamba env create -f environment.yml)."
        ) from e

    import rpy2.robjects as ro

    if _copykat_available():
        return True

    if not quiet:
        print(f"Installing the copykat R package from GitHub ({COPYKAT_REPO}) — one-time...",
              flush=True)
    try:
        # dependencies=FALSE: copykat's CRAN deps are already provided by conda.
        ro.r(f'remotes::install_github("{COPYKAT_REPO}", dependencies=FALSE, upgrade="never")')
    except Exception as e:  # rpy2 RRuntimeError etc.
        raise RuntimeError(
            f"Failed to install the copykat R package ({COPYKAT_REPO}). Check network "
            f"access and R, then re-run the CopyKAT analysis — the install is retried "
            f"automatically. Underlying error: {e}"
        ) from e

    if not _copykat_available():
        raise RuntimeError(
            f"copykat R package still not importable after install_github({COPYKAT_REPO})."
        )
    if not quiet:
        print("copykat R package installed.", flush=True)
    return True


def main() -> int:
    try:
        ensure_copykat_installed()
        print("copykat is available.")
        return 0
    except Exception as e:  # noqa: BLE001
        print(f"error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
