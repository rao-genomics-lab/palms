# Installation

## Prerequisites

- **Operating system**: Linux, macOS (Apple Silicon or Intel), or WSL2. Linux is the primary development and CI platform; macOS is supported but less exercised. Native Windows (outside WSL) is not supported.
- **Memory**: 8 GB RAM minimum; 32 GB or more is recommended for large datasets (500k+ cells)
- **Package manager**: `pip` (Python 3.10+) for the released package, or
  [mamba](https://mamba.readthedocs.io/) / [conda](https://docs.conda.io/) for the
  development environment. Only the CopyKAT CNV backend requires conda — see below.

## Install from PyPI

```bash
pip install palms
palms /path/to/xenium/output/
```

That is the whole install for using the viewer. PyQt6 brings Qt inside its wheel, so
there is no system Qt or GL package to add, and nothing to clone. Verified on Linux and
macOS against the published 1.0.0.

Add the CNV stack, or everything, with an extra:

```bash
pip install "palms[cnv]"     # inferCNV backend
pip install "palms[full]"    # every optional feature
```

One thing a wheel install does not get: the **CopyKAT** CNV backend. It needs a second
conda environment carrying rpy2, R 4.3 and the `copykat` R package, so it is reachable
only from a conda install (and, on Apple Silicon, not at all — see below). The
**inferCNV** backend runs in-process and works everywhere, so CNV inference itself is not
lost.

## Install from source (conda)

Use this to develop PALMS, to run the test suite, or to get the CopyKAT backend.

```bash
git clone https://github.com/rao-genomics-lab/palms.git
cd palms
./scripts/install.sh
conda activate palms
```

This installs the `palms` package in editable mode along with all required dependencies. It also installs the CNV inference stack (`infercnvpy` and `insitucnv-copykat`), so the CNV tab's **inferCNV** backend works immediately — no extra step required.

### What `install.sh` does, and the manual equivalent

It is `conda env create` plus one OS-dependent step. On **Linux and WSL** it additionally applies `environment-linux.yml`, which adds `libglx-devel` — the package that fixes the `Could not initialize GLX` abort described below. That package has no macOS build, so keeping it in `environment.yml` made the environment impossible to solve on a Mac; it lives in a separate overlay file instead. conda environment files have no platform selectors (`# [linux]` is a conda-*build* feature that `conda env create` ignores), which is why this needs a script rather than a conditional line.

To run conda yourself:

```bash
conda env create -f environment.yml                          # all platforms
conda env update -n palms -f environment-linux.yml   # Linux and WSL only
```

Forgetting the second command on Linux is not silent — the viewer detects it at startup and prints the command to run, instead of aborting with no traceback.

`install.sh` passes any extra arguments through to `env create`, so `./scripts/install.sh --name xv-test` works and the overlay follows the name you chose. Set `CONDA_EXE_OVERRIDE` to force a particular solver binary; otherwise it prefers `mamba` and falls back to `conda`.

## Optional Extras

These are the extras declared in `pyproject.toml`. From PyPI they are
`pip install "palms[<extra>]"`; in a checkout, with the environment activated, the same
extra is `pip install -e ".[<extra>]"`.

| Extra | Install command | What it adds |
|-------|----------------|--------------|
| CellTypist | `pip install "palms[celltypist]"` | Automated cell type annotation using pre-trained models |
| R integration | `pip install "palms[r]"` | rpy2 + anndata2ri, for the R-based reference fetcher (needs a system R) |
| GPU support | `pip install "palms[gpu]"` | torch, torch-geometric, pytorch-lightning, xgboost |
| Reference datasets | `pip install "palms[references]"` | readfcs + rasterio, for reading reference panel formats |
| CNV inference | `pip install "palms[cnv]"` | infercnvpy + insitucnv-copykat — **already included** in `environment.yml` |
| Everything | `pip install "palms[full]"` | All of the above |
| Test suite | `pip install "palms[test]"` | pytest + nbformat, nbclient and ipykernel, needed to run the notebook-replay test |
| Development | `pip install "palms[dev]"` | The `test` extra plus ruff, matching what CI runs |

`full` deliberately does *not* include `test` or `dev` — install those explicitly if you intend to run the suite. A plain `pip` install needs Python 3.10 or newer; the conda environment pins 3.12.

Each extra is independent — install only what you need. A tab whose optional dependency is missing still appears in the UI; it reports a clear "not installed" error naming the command to run, rather than failing at startup.

The `cnv` extra is the exception to "optional": `environment.yml` already installs it, so you only need it for a plain `pip` install. It pulls `insitucnv-copykat` — the [`insituCNV-copykat`](https://github.com/sraorao/insituCNV-copykat) fork, published under that name because `insitucnv` on PyPI is upstream's. It imports as `insitucnv` either way, so don't install both. The fork exists because upstream's release pins `anndata<0.12` and `pandas<3`, which cannot resolve against this app.

To download pre-built reference panels, use the `palms-fetch-references` command (see [Console Scripts](#console-scripts) below) — it ships with the package and needs no extra install.

## CopyKAT Backend (Optional Second Environment)

The CNV tab's inferCNV backend runs in the main environment. Its **CopyKAT** backend does not: CopyKAT needs rpy2 with R 4.3, a stack that only builds on python 3.11, which cannot coexist with the main environment's python 3.12. It therefore lives in a separate environment:

```bash
conda env create -f environment-copykat.yml   # creates 'palms_copykat'
```

**Linux only.** This environment does not solve on Apple Silicon: `r-dlm`, a CopyKAT dependency, is published only on the Anaconda `r` channel, which has no `osx-arm64` builds. The inferCNV backend runs in the main environment and works on every supported platform.

The viewer locates that environment by name and launches CopyKAT there as a detached background job, so the analysis (typically a couple of hours) survives closing the GUI. The `copykat` R package is GitHub-only and installs automatically on first run.

To use a differently-named environment, set one of:

| Variable | Meaning |
|----------|---------|
| `PALMS_COPYKAT_ENV` | Name of the conda environment to use |
| `PALMS_COPYKAT_PYTHON` | Full path to that environment's `python` |

Skip this section entirely if you only need inferCNV.

## Reinstalling

If the conda environment becomes broken or corrupted, remove it and start fresh:

```bash
conda env remove -n palms
git clone https://github.com/rao-genomics-lab/palms.git
cd palms
./scripts/install.sh
conda activate palms
```

Your analyses stored in the dataset directory (`sdata_cached.zarr/viewer_session/`) are not affected by reinstalling the environment.

## Launching the Viewer

```bash
# Optional one-time transcript preprocessing per dataset (~30-60 min)
# Speeds up per-gene transcript loading from ~5 s to <100 ms
palms-preprocess /path/to/xenium/output/

# Launch viewer (a file dialog opens if no path is given)
palms [/path/to/xenium/output/]

# Launch without building or reading the SpatialData zarr cache.
# Expect a slow start and high memory use — see "Memory" under Troubleshooting.
palms /path/to/xenium/output/ --no-cache

# Launch ignoring any customised analysis templates, running the shipped ones
palms /path/to/xenium/output/ --no-user-templates
```

`python -m palms` is equivalent to `palms` and accepts the same arguments.

### Building the zarr cache ahead of time

The first launch on a dataset reads the raw Xenium output and writes `sdata_cached.zarr/` beside it, which takes tens of minutes and tens of GB. That happens automatically, so this step is optional — but it holds a napari window open for the duration, which is awkward over ssh. `palms-build-cache` does the same work with no GUI:

```bash
# Build it (or refresh it) — safe to run detached or overnight
palms-build-cache /path/to/xenium/output/

# Report on an existing cache and exit; builds nothing, writes nothing.
# Exits non-zero if the cache is missing, stale, or does not verify.
palms-build-cache /path/to/xenium/output/ --check
```

With no GUI attached there is nobody to answer the "this cache looks stale — rebuild it?" question, so by default the command **keeps the existing cache** rather than discarding work you cannot get back. `--on-stale` answers in advance:

| `--on-stale` | Meaning |
|--------------|---------|
| `ask` (default) | Prompt if a GUI is available, otherwise keep the cache untouched |
| `restore` | Rebuild and carry your ROIs, registrations, clusterings and CNV results over |
| `rebuild` | Rebuild and discard them (the old cache is still moved aside, never deleted) |
| `keep` | Load the cache as it is, and stop asking |

Run `--check` first if you are not sure what a rebuild would cost you — it lists the user-generated data the cache holds. This is separate from `palms-preprocess`, which builds the per-gene transcript cache; the two are independent and can be run in either order.

### Environment variables

| Variable | Meaning |
|----------|---------|
| `PALMS_COPYKAT_ENV` | Name of the conda environment holding the CopyKAT stack |
| `PALMS_COPYKAT_PYTHON` | Full path to that environment's `python` |
| `PALMS_TEMPLATE_PATH` | Colon-separated search path replacing the default location for customised analysis templates. Set it to an empty value to disable overrides entirely. See [Templates](Tab-Templates). |

## Console Scripts

The following commands are available after activating the environment:

| Command | Description |
|---------|-------------|
| `palms` | Launch the viewer |
| `palms-preprocess` | Preprocess transcripts for fast per-gene loading |
| `palms-build-cache` | Build or inspect the SpatialData zarr cache without the GUI |
| `palms-fetch-references` | Download pre-built label transfer reference datasets |
| `palms-build-custom-segmentation` | Run the custom cell segmentation pipeline |

All commands accept `--help` for usage details.

The repository also ships `scripts/verify_notebook.py`, which exports a dataset's recorded analysis as a notebook, replays it against the raw Xenium output, and compares the replayed results against what the viewer saved. Use `--dry-run` for a summary in seconds without replaying.

## Troubleshooting

**ICE/X11 disconnect on startup**
Handled automatically at startup. No action required.

**Permission errors on read-only filesystems**
If your Xenium output directory is on a read-only filesystem, copy the dataset to a writable location before launching. `--no-cache` will also start, but on a full slide it is expensive rather than merely slower — see **Memory** below — so copying the dataset is the better answer where there is room for it.

**Memory**
Building the cache for a full slide peaks at around 9 GB, and reports its usage per element as it goes, to both the terminal and the session log. Once a cache exists, launching against it is a fraction of that.

On **macOS** those per-element figures read as unknown. The probe uses `/proc/self/statm` and `/proc/self/status`, which macOS does not have; it fails soft, so this costs you the numbers and nothing else. The peak memory itself is unchanged, so the guidance in this section still applies.

`--no-cache` is the exception, and the reason is worth understanding. The image pyramid the viewer displays exists on disk in a cache; without one it has to be *computed*, and napari draws the smallest, most zoomed-out level first — which stands on every chunk of the full-resolution image below it. On a 57887×51217 slide that is tens of GB, and it is not something a smaller machine can page its way through: the computation is not streamable, so a whole pyramid level must be resident at once. The viewer warns at startup when it is about to do this. Use `--no-cache` for a quick look at a small dataset, not as a way to avoid writing to disk on a large one.

**The viewer vanished with no error, and took the terminal with it** (Linux)
That is almost certainly not a crash — a crashing process leaves a traceback, and closing every tab in a terminal window is beyond what one process can do to itself. On systemd-based Linux desktops it is `systemd-oomd`, which kills an entire cgroup scope when a slice is under sustained memory *pressure*, not when memory is exhausted. It is why you can see this at 60% of RAM and never catch the process anywhere near the limit.

```bash
journalctl -u systemd-oomd --since "1 hour ago" | grep Killed
```

A line naming a `vte-spawn-*.scope` and a pressure percentage confirms it. Nothing appears in `dmesg`, because this is not the kernel OOM killer. The fix is to reduce the memory rather than to disable oomd, which is what keeps a runaway allocation from taking the desktop down with it: build the cache once and launch against it rather than using `--no-cache`, and close other large applications during a first load.

**Corrupt zarr cache**
The viewer tries to repair the cache itself before asking you anything — replaying any interrupted write, clearing debris, rebuilding the metadata index, and restoring elements from its own backup copies. Most damage is resolved there and you never see a dialog.

If the repair does not succeed, a dialog offers three choices:

| Choice | What happens |
|--------|--------------|
| Rebuild and restore my data | The cache is set aside as `sdata_cached_backup_<timestamp>.zarr`, a fresh one is built from the raw files, and your clusterings, ROIs, registrations and session state are copied across. The backup is removed once the rebuild succeeds. |
| Rebuild without restoring | The cache is set aside as `sdata_cached_prev_<timestamp>.zarr` and **kept**, and a fresh empty cache is built. Use this when you suspect the saved analyses themselves are the problem. |
| Quit | Nothing is changed, so you can back the directory up first. |

The viewer never deletes a cache — every path here is a rename. With no display available it raises an error rather than choosing for you. To inspect or repair a cache deliberately rather than at startup, use the [Cache](Tab-Cache) tab; for the whole procedure, see [Recovering a Cache](Tutorial-Recovering-a-Cache).

**Stale zarr cache**
Freshness is decided by a content hash of `experiment.xenium` recorded when the cache was built, so copying a dataset, re-downloading it or restoring it from backup does not make the cache look stale — only a genuine change to the file does. If a change is detected, a dialog at startup offers to rebuild; your existing analyses are preserved.

Caches built before content hashing existed fall back to comparing modification times. That is treated as an uncertain signal: it prompts rather than rebuilding, and choosing to keep the cache records a hash so it stops asking.

**Where the viewer writes**
Besides the zarr cache, a dataset directory accumulates `viewer_cache/` (normalised expression, CNV profiles, cached DEG tables, the analysis provenance graph), `plots/` (saved figures), `transcript_cache/` (per-gene transcript files), `analysis.py`, `analysis_notebook.ipynb`, and a rotating `palms.log`. The [Dataset](Tab-Dataset) tab lists all of it with sizes and can delete the regenerable parts.


## Troubleshooting: `Could not initialize GLX` / `Aborted (core dumped)`

**Linux and WSL only** — macOS has no GLX at all; Qt uses the `cocoa` platform plugin there.

If the viewer dies immediately with

```
WARNING: Could not initialize GLX
Aborted (core dumped)
```

your environment is missing **`libglx-devel`**. That happens if the environment was created
with a bare `conda env create -f environment.yml` instead of `./scripts/install.sh`, or
predates the package being added at all. Apply the Linux overlay:

```bash
conda env update -n palms -f environment-linux.yml
```

The viewer checks for this at startup and prints a fix naming *your* environment, so you
should not have to diagnose the abort from the message above. A wheel-only install
(`pip install palms`) cannot hit this and is not warned: PyQt6 brings Qt inside its wheel and
nothing pulls conda's `libglx`, so there is only one copy of the library in the process.

Why: conda ships `libGLX.so.0` but not the unversioned `libGLX.so`. PyOpenGL's loader looks
for the unversioned name first, misses the environment, and loads your *system's* copy
instead — leaving two different builds of the same library in one process. Qt6 then resolves
GL calls across both and aborts. There is no Python traceback because the abort comes from
Qt itself, and no environment variable can work around it: the fix has to put the missing
name inside the environment.

Most visible on remote desktops (ThinLinc/VNC/x2go/xrdp) on a machine that also has the
system `libglx-dev` package, which is a common combination.
