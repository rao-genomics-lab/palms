#!/usr/bin/env bash
#
# Create the xenium_viewer conda env, applying the Linux-only GL overlay when the
# host needs it. One command on Linux, macOS and WSL:
#
#     ./scripts/install.sh
#     conda activate xenium_viewer
#
# Extra arguments are passed through to `env create`, e.g.:
#
#     ./scripts/install.sh --name xv-test
#
# The OS branch exists because conda env files have no platform selectors, and
# `libglx-devel` — which fixes a Qt6 startup abort on remote X displays — is a
# linux-only package. See environment-linux.yml for the root cause.
#
# Deliberately POSIX-ish bash: macOS still ships bash 3.2, so no arrays, no
# `mapfile`, no `${var,,}`.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"

# mamba is much faster than conda and solves the same env; prefer it if present.
# Override with CONDA_EXE_OVERRIDE=/path/to/micromamba for CI or an unusual install.
if [ -n "${CONDA_EXE_OVERRIDE:-}" ]; then
    conda_exe="$CONDA_EXE_OVERRIDE"
elif command -v mamba >/dev/null 2>&1; then
    conda_exe="mamba"
elif command -v conda >/dev/null 2>&1; then
    conda_exe="conda"
else
    echo "install.sh: neither mamba nor conda is on PATH." >&2
    echo "Install miniforge/mambaforge first: https://github.com/conda-forge/miniforge" >&2
    exit 1
fi

echo "==> Creating the environment with $conda_exe"
"$conda_exe" env create -f "$repo/environment.yml" "$@"

# Determine the env name we just created, so a passed-through --name/-n still
# gets the overlay. `env create` has no way to report it, so re-read the file.
env_name="$(sed -n 's/^name:[[:space:]]*//p' "$repo/environment.yml" | head -1)"
prev=""
for arg in "$@"; do
    case "$prev" in
        --name|-n) env_name="$arg" ;;
    esac
    prev="$arg"
done

case "$(uname -s)" in
    Linux)
        # Covers WSL too: `uname -s` is Linux there and conda installs linux-64
        # packages, so WSL needs the overlay exactly as a native Linux box does.
        echo "==> Linux detected: applying the GL overlay (environment-linux.yml)"
        "$conda_exe" env update -n "$env_name" -f "$repo/environment-linux.yml"
        ;;
    Darwin)
        echo "==> macOS detected: skipping the Linux GL overlay (no GLX on macOS)"
        ;;
    *)
        echo "==> $(uname -s) is not a tested platform; skipping the Linux GL overlay" >&2
        ;;
esac

echo
echo "Done. Activate it with:"
echo "    conda activate $env_name"
