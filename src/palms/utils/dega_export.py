"""Export a dataset as Celldega DegaFiles — a shareable, web-native view.

A PALMS session already ends as a replayable notebook. This makes it end as a
*viewable* artifact too: DegaFiles are WebP Deep Zoom image pyramids plus Apache
Parquet vector data, which Celldega's ``Landscape`` widget renders in a browser
with no server and no install. Celldega's own Discussion asks for exactly this
interoperability, and the conversion is theirs — ``celldega.pre.main`` — so this
module is a *wrapper*, not a second implementation of their format.

Export only. Reading DegaFiles is deliberately out of scope: that is a second
format reader, and the viewer already reads the raw 10x output directly.

Three things this wrapper exists to do, none of which the bare API does:

**The raw 10x output is never written to.** ``celldega.pre.main`` calls
``_xenium_unzipper``, which ``os.chdir``s into the dataset directory and runs
``gzip -dk cells.csv.gz`` and ``unzip cells.zarr.zip`` *in place*, leaving several
GB of decompressed duplicates beside the user's data — and failing outright on a
read-only mount, which a shared NAS dataset usually is. Nothing else in this
package writes to the raw output; ``store_inventory.deletable_roots`` will not
even let the user delete a file there. So the export runs against a directory of
**symlinks** to the raw entries, in ``viewer_cache/``: the reads follow the links
to the real files, and every write lands in a directory the viewer owns and the
Dataset tab can clear. The archives are then unpacked into that directory by
:func:`extract_archives`, which is not an optimisation: ``gzip -dk`` refuses to
read through a symlink, so the staged run does not work without it.

**It says which optional dependency is missing, and how to get it.** See
``dega_available`` and :func:`install_hint`; the constraint is sharper than
"pip install celldega" — sharper, it turned out, than ``celldega[pre]`` too,
since ``--no-deps`` makes pip ignore extras — and a bare ``ImportError`` in a
napari worker would not say any of it.

**It refuses a dataset that has nothing to convert.** A Crop Dataset export
carries ``experiment.xenium``, a zarr cache and derived transcripts, and that is
the whole dataset: its zarr *is* the data. celldega reads the raw 10x bundle, so
the export cannot work, and without :func:`require_exportable` it fails minutes
in with ``CalledProcessError: Command '['gzip', '-dk', 'cells.csv.gz']' returned
non-zero exit status 1`` — measured on the ``crop_6`` demo dataset — which reads
as a broken archive rather than as the real answer. This is the rule issue #17
established for every *rebuild* path, applied to the one path that reads the raw
output instead of the cache.

Measured 2026-09-03 on ``Xenium_V1_human_Pancreas_FFPE_outs`` with
celldega 0.24.2, twice on the same machine: a complete export in **322 s and
then 637 s** — budget five to eleven minutes rather than a number — to 220 MB
(``pyramid_images/`` as WebP Deep Zoom, ``transcript_tiles/``,
``cell_segmentation/``, ``cbg/``, ``cell_clusters/``, ``cell_metadata.parquet``,
``meta_gene.parquet``, ``df_sig.parquet``, ``micron_to_image_transform.csv``,
``landscape_parameters.json``), with the raw output directory byte-for-byte
unchanged, verified on disk afterwards: ``cells.csv``, ``cells.zarr``,
``analysis`` and ``cell_feature_matrix`` all appeared in the staging farm and
none of them beside the data. Both runs were against anndata 0.13.2,
spatialdata 0.8.0 and pandas 3.0.5 — over every upper cap celldega declares —
which is the evidence for the ``--no-deps`` advice below.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Optional

#: Default Celldega tile size, in microns. Theirs, from the tutorials.
TILE_SIZE_UM = 250

#: Image layers to tile, matching celldega's own default —
#: ``run_pre_processing.main(image_tile_layer="all")``. For Xenium that is the
#: four channels ``get_image_info('Xenium', 'all')`` names — dapi, bound, rna,
#: prot — which are exactly the four ``morphology_focus_000{0..3}.ome.tif``
#: files a bundle ships, and what their Landscape viewer expects to offer.
#:
#: This was ``"dapi"`` until 2026-09-05, on the reading that the nuclear channel
#: alone "is what their Xenium tutorial uses". That is only half true — their
#: docstring example does pass ``--image_tile_layer 'dapi'`` and
#: ``get_image_info`` defaults to it — but their pipeline entry point does not,
#: and a published export that silently drops three of four channels is the
#: wrong default for a tool whose point is that the artifact is the record.
#: The cost is real (four pyramids, so roughly 4x the runtime and output size),
#: which is why the tab still offers 'dapi' and now says what the choice buys.
IMAGE_TILE_LAYER = "all"

#: The pin. A floating version here would be the ``insitucnv`` mistake twice —
#: see the note in CLAUDE.md about why that git URL became a PyPI pin.
#:
#: **Not** ``celldega[pre]``. ``--no-deps`` makes pip ignore extras entirely, so
#: that spelling installs exactly the same thing while *reading* as though it
#: brought pyvips along — measured with ``pip install --dry-run --no-deps
#: 'celldega[pre]==0.24.2'``: "Would install celldega-0.24.2", and nothing else.
#: The remedy below therefore names every dependency by hand.
CELLDEGA_PIN = "celldega==0.24.2"

#: celldega's runtime imports that ``--no-deps`` skips, found by importing
#: ``celldega.pre`` and installing whatever it asked for until it stopped asking.
#: None of them is version-capped by celldega, so they resolve normally — it is
#: only anndata, spatialdata, pandas and ome-zarr that must be kept away from
#: pip's resolver, which is what ``--no-deps`` on the celldega line is for.
CELLDEGA_RUNTIME_DEPS = (
    "mudata", "ipywidgets", "anywidget", "libpysal", "polars",
)

#: pyvips, and the wheel that carries its own libvips. The second is a
#: safeguard against an **order-dependent** collision, and the ordering half of
#: the story is what :func:`dega_available` has to keep honouring.
#:
#: Plain ``pyvips`` runs in API mode against whatever ``libvips.so.42`` the host
#: provides. On a box with Ubuntu's libvips that build pulls the system HDF5
#: 1.10 into the process, while the conda env's h5py is built for 1.14 — two
#: builds of one library in one process, the ``libglx`` collision in
#: ``utils/gl_check.py`` wearing a different hat. Measured here, without
#: ``pyvips-binary``: pyvips loads ``/usr/lib/x86_64-linux-gnu/libvips.so.42``
#: and ``/usr/lib/x86_64-linux-gnu/libhdf5_serial.so.103``.
#:
#: What it does depends entirely on **who got there first**, and this was
#: measured both ways:
#:
#: * ``import pyvips`` *then* ``import h5py`` → ``ValueError: Not a datatype``.
#:   h5py is unusable, so anndata is, for the rest of the process.
#: * ``import h5py`` *then* ``import pyvips`` → fine. Both libraries are mapped,
#:   but h5py already bound to the env's HDF5 and keeps working; an h5ad
#:   round-trip afterwards succeeds.
#:
#: The viewer is always in the second case (anndata is imported at startup), and
#: so is :func:`dega_available`, which imports ``celldega.pre`` — and therefore
#: anndata, and therefore h5py — *before* pyvips. That ordering is load bearing
#: rather than stylistic, and ``tests/test_dega_export.py`` pins it.
#:
#: ``pyvips-binary`` removes the hazard rather than sequencing around it: it
#: ships its own libvips wheel, which pyvips prefers, so the host build is never
#: mapped at all. Recommended for that reason, not because the GUI needs it.
#: conda-forge's ``libvips`` would also fix it and is **not** the route —
#: solving it into the live env downgrades napari 0.8→0.7, spatialdata
#: 0.8→0.6.1 and squidpy 1.8.2→1.7.
PYVIPS_PINS = ("pyvips~=2.2.2", "pyvips-binary")


def install_hint() -> str:
    """The install command that actually works, as one block of shell."""
    return (f"    pip install --no-deps {CELLDEGA_PIN}\n"
            f"    pip install {' '.join(PYVIPS_PINS + CELLDEGA_RUNTIME_DEPS)}")


class DegaUnavailable(RuntimeError):
    """Celldega is not importable, with the reason and the remedy attached."""


class NotExportable(RuntimeError):
    """The dataset has no raw 10x output, so there is nothing to convert."""


def dega_available() -> tuple[bool, str]:
    """``(importable, message)`` — the message is empty when it is.

    Not a bare ``importlib.util.find_spec``: celldega caps ``anndata<0.13`` and
    ``spatialdata<0.8`` while this package runs 0.13 and 0.8, so
    ``pip install celldega`` into the viewer's environment *succeeds* and
    silently downgrades anndata, spatialdata, pandas and ome-zarr underneath it.
    Measured 2026-09-03 against celldega 0.24.2: pip would install
    anndata 0.12.19, spatialdata 0.7.3, pandas 2.3.3, ome-zarr 0.15.0. Those
    caps are conservative rather than real — celldega imports and runs against
    the viewer's own versions — so the remedy is ``--no-deps``, and saying that
    is most of this function's value.
    """
    try:
        # ORDER IS LOAD BEARING: celldega.pre first, pyvips second. Importing
        # celldega pulls in anndata and therefore h5py, which binds to this
        # env's HDF5 before pyvips can map the host libvips and its own. The
        # reverse order raises `ValueError: Not a datatype` out of h5py and
        # leaves anndata broken for the rest of the process -- measured, and see
        # PYVIPS_PINS. Do not sort these two imports.
        import celldega.pre  # noqa: F401

        # celldega imports pyvips in a try/except that leaves the module bound
        # to None, so an install without it reaches "generating dapi image
        # tiles" and dies there with an AttributeError. Checking here turns a
        # failure ten minutes into an export into one before it starts.
        import pyvips  # noqa: F401
    except ImportError as exc:
        return False, (
            f"Celldega is not installed ({exc}). Two commands, and the split "
            f"matters:\n"
            f"{install_hint()}\n"
            f"The first is --no-deps because celldega caps anndata<0.13 and "
            f"spatialdata<0.8 while the viewer runs 0.13 and 0.8, so a plain "
            f"install succeeds and downgrades the stack underneath it; those "
            f"caps are conservative and it runs against the versions already "
            f"here. The second puts back the dependencies --no-deps skipped — "
            f"note that pip ignores extras under --no-deps, so 'celldega[pre]' "
            f"would NOT have brought pyvips. pyvips-binary is strongly "
            f"advised: without it pyvips maps the host libvips, whose HDF5 can "
            f"collide with this env's h5py depending on import order."
        )
    return True, ""


def require_dega() -> None:
    """Raise :class:`DegaUnavailable` with the remedy, or return."""
    ok, message = dega_available()
    if not ok:
        raise DegaUnavailable(message)


#: What ``celldega.pre.main`` reads that only raw 10x output has. A Crop Dataset
#: export has none of them — it ships ``experiment.xenium`` plus the zarr cache
#: and derived transcripts, and that *is* the whole dataset.
_REQUIRED_RAW = ("cells.csv.gz", "cells.zarr.zip", "cell_feature_matrix.tar.gz",
                 "morphology_focus")


def missing_raw_inputs(data_path) -> list[str]:
    """Which of celldega's raw inputs *data_path* does not have.

    Empty means the export can run. Named individually rather than answered with
    a bool because the two cases need different advice: a Crop Dataset export is
    missing all of them and can never be published, while a bundle missing one is
    an incomplete download.
    """
    data_path = Path(data_path)
    return [name for name in _REQUIRED_RAW if not (data_path / name).exists()]


def is_cache_only(data_path) -> bool:
    """True when *none* of celldega's raw inputs is present.

    The distinction :func:`missing_raw_inputs` draws, as the predicate a caller
    usually wants: this dataset can never be published, as against one that is
    merely an incomplete copy and could be. Mirrors
    ``loader.has_raw_xenium_source``'s "all or nothing" reading, and for the same
    reason — a bundle missing *some* files is broken raw output, not a crop.
    """
    return len(missing_raw_inputs(data_path)) == len(_REQUIRED_RAW)


def require_exportable(data_path) -> None:
    """Raise :class:`NotExportable` unless the raw 10x output is here.

    Checked *before* anything is staged, because the failure without it is both
    late and misleading: ``celldega.pre.main`` gets as far as its own unzipper
    and dies with ``CalledProcessError: Command '['gzip', '-dk',
    'cells.csv.gz']' returned non-zero exit status 1`` — measured on the
    ``crop_6`` demo dataset. That reads as a broken gzip or a corrupt archive
    rather than as the real answer, which is that this dataset has no raw output
    and never will.

    This is the rule issue #17 established for every rebuild path, applied to the
    one path that reads the raw output *instead of* the cache: a crop export's
    zarr **is** the data. ``loader.has_raw_xenium_source`` is the general form;
    the check here is narrower on purpose, since celldega needs specific files
    that a merely *partial* bundle could still be missing.
    """
    missing = missing_raw_inputs(data_path)
    if not missing:
        return
    data_path = Path(data_path)
    if is_cache_only(data_path):
        raise NotExportable(
            f"{data_path.name} has no raw 10x output, so there is nothing for "
            f"celldega to convert: it is a Crop Dataset export (or another "
            f"cache-only dataset), whose SpatialData zarr *is* the data.\n"
            f"DegaFile export reads the original 10x bundle — "
            f"{', '.join(_REQUIRED_RAW)} — not the cache, so it cannot run "
            f"here. Publish the dataset this one was cropped from instead."
        )
    raise NotExportable(
        f"{data_path.name} is missing {', '.join(missing)}, which celldega "
        f"needs. This looks like an incomplete copy of a 10x bundle rather "
        f"than a Crop Dataset export — the rest of the raw output is here."
    )


def degafiles_dir(data_path) -> Path:
    """Where an export goes: ``<data_path>/degafiles``.

    Beside ``plots/`` rather than inside ``viewer_cache/``: this is a deliverable
    the user shares, not a cache the viewer may clear.
    """
    return Path(data_path) / "degafiles"


def staging_dir(data_path) -> Path:
    """The symlink farm celldega is pointed at, under ``viewer_cache/``.

    Named after the dataset, one level down, because celldega's contract is a
    sample *name* plus its parent directory rather than one path -- and that
    name is what it stamps into ``landscape_parameters.json``. A flat
    ``dega_staging`` would publish DegaFiles whose sample is called
    "dega_staging".
    """
    data_path = Path(data_path)
    return data_path / "viewer_cache" / "dega_staging" / data_path.name


def stage_raw_output(data_path, staging: Optional[Path] = None) -> Path:
    """Link every entry of the dataset directory into a directory we own.

    Returns the staged directory. Existing links are refreshed, and anything
    celldega decompressed on a previous run is left in place — its own
    "skip if it already exists" checks then make a second export much cheaper.

    Symlinks and not copies: ``transcripts.parquet`` alone is 1.4 GB on a full
    slide, and celldega only reads it.
    """
    data_path = Path(data_path)
    staging = Path(staging) if staging is not None else staging_dir(data_path)
    staging.mkdir(parents=True, exist_ok=True)

    for entry in sorted(data_path.iterdir()):
        # Not the viewer's own directories: linking viewer_cache into a child of
        # itself is a loop, and the rest is ours rather than 10x's.
        if entry.name in {"viewer_cache", "degafiles", "plots", "transcript_cache"}:
            continue
        if entry.name.startswith("sdata_cached"):
            continue
        link = staging / entry.name
        if link.is_symlink():
            if link.readlink() == entry:
                continue
            link.unlink()
        elif link.exists():
            # A real file here is something celldega extracted on a previous
            # run. It is ours, it is correct, and re-extracting it is minutes.
            continue
        link.symlink_to(entry)
    return staging


#: What ``celldega.pre._xenium_unzipper`` would extract, and from what. Kept as
#: data because it is *their* list, copied from celldega 0.24.2 -- re-check it
#: when the pin moves. Entries are (produced, archive, kind).
_ARCHIVES = (
    ("cells.csv", "cells.csv.gz", "gzip"),
    ("cells.zarr", "cells.zarr.zip", "zip"),
    ("analysis", "analysis.tar.gz", "tar"),
    ("cell_feature_matrix", "cell_feature_matrix.tar.gz", "tar"),
)


def extract_archives(staging: Path) -> list[str]:
    """Decompress into *staging* what celldega would otherwise decompress in place.

    Returns the names produced; an already-extracted one is skipped and not
    listed, which is what makes a second export cheap.

    This is not merely a nicer place to put the output. ``gzip -dk`` **refuses a
    symlink** -- "not a regular file", exit 1 -- so linking the archives and
    letting celldega extract them does not work at all; that is how this
    function came to exist. Doing it in Python also drops celldega's unstated
    dependency on the ``gzip``, ``unzip`` and ``tar`` binaries, and leaves
    ``_xenium_unzipper`` a no-op, since it skips any target that already exists.

    ``tarfile`` is given ``filter='data'`` -- the default from Python 3.14, a
    warning before it -- so an archive member cannot write outside *staging*.
    """
    import gzip
    import tarfile
    import zipfile

    produced = []
    for target, archive, kind in _ARCHIVES:
        out = staging / target
        source = staging / archive
        if out.exists() or not source.exists():
            continue
        if kind == "gzip":
            with gzip.open(source, "rb") as src, open(out, "wb") as dst:
                shutil.copyfileobj(src, dst)
        elif kind == "zip":
            with zipfile.ZipFile(source) as zf:
                zf.extractall(out)
        else:
            with tarfile.open(source) as tf:
                tf.extractall(staging, filter="data")
        produced.append(target)
    return produced


def clear_staging(data_path) -> None:
    """Remove the staging directory, decompressed copies and all."""
    staging = staging_dir(data_path).parent
    if staging.exists():
        shutil.rmtree(staging)


def export_degafiles(data_path, out_dir=None, tile_size: int = TILE_SIZE_UM,
                     image_tile_layer: str = IMAGE_TILE_LAYER,
                     use_int_index: bool = True, max_workers: int = 1) -> Path:
    """Write DegaFiles for *data_path* and return the directory they went to.

    Thin over ``celldega.pre.main``, whose contract is ``(sample,
    data_root_dir)`` — a directory *name* and its parent — rather than one path.
    Here that pair names the staged symlink farm, which is why
    :func:`staging_dir` puts it one level down under the dataset's own name: the
    sample Celldega stamps into ``landscape_parameters.json`` is then the name a
    reader expects rather than an internal directory of ours.

    Note the process-global side effect this cannot remove: celldega's unzipper
    ``os.chdir``s and restores in a ``finally``, so for the duration of that step
    every relative path in this process resolves somewhere else. It is why the
    cwd is asserted afterwards rather than assumed, and it is the strongest
    argument for eventually running the export in a subprocess.
    """
    require_dega()
    data_path = Path(data_path).resolve()
    # Before the output directory is made and before anything is staged: a
    # dataset that cannot be exported should leave no trace of having tried.
    require_exportable(data_path)
    import celldega.pre as dega_pre

    out = Path(out_dir) if out_dir is not None else degafiles_dir(data_path)
    out.mkdir(parents=True, exist_ok=True)

    staging = stage_raw_output(data_path)
    extract_archives(staging)
    before = Path.cwd()
    try:
        dega_pre.main(
            sample=staging.name,
            data_root_dir=str(staging.parent),
            tile_size=tile_size,
            image_tile_layer=image_tile_layer,
            path_dega_files=str(out),
            use_int_index=use_int_index,
            max_workers=max_workers,
        )
    finally:
        if Path.cwd() != before:
            os.chdir(before)
    return out
