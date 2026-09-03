"""Export a dataset as Celldega DegaFiles — a shareable, web-native view.

A PALMS session already ends as a replayable notebook. This makes it end as a
*viewable* artifact too: DegaFiles are WebP Deep Zoom image pyramids plus Apache
Parquet vector data, which Celldega's ``Landscape`` widget renders in a browser
with no server and no install. Celldega's own Discussion asks for exactly this
interoperability, and the conversion is theirs — ``celldega.pre.main`` — so this
module is a *wrapper*, not a second implementation of their format.

Export only. Reading DegaFiles is deliberately out of scope: that is a second
format reader, and the viewer already reads the raw 10x output directly.

Two things this wrapper exists to do, both of which the bare API does not:

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
``dega_available``; the constraint is sharper than "pip install celldega", and a
bare ``ImportError`` in a napari worker would not say so.

Measured 2026-09-03 on ``Xenium_V1_human_Pancreas_FFPE_outs`` with
celldega 0.24.2: a complete export in **322 s** to 220 MB
(``pyramid_images/`` as WebP Deep Zoom, ``transcript_tiles/``,
``cell_segmentation/``, ``cbg/``, ``cell_clusters/``, ``cell_metadata.parquet``,
``meta_gene.parquet``, ``df_sig.parquet``, ``micron_to_image_transform.csv``,
``landscape_parameters.json``), with the raw output directory byte-for-byte
unchanged. That run was against anndata 0.13.2, spatialdata 0.8.0 and
pandas 3.0.5 — over every upper cap celldega declares — which is the evidence
for the ``--no-deps`` advice below.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Optional

#: Default Celldega tile size, in microns. Theirs, from the tutorials.
TILE_SIZE_UM = 250

#: Image layers to tile. 'dapi' is the nuclear channel alone and is what their
#: Xenium tutorial uses; 'all' tiles every morphology channel and costs
#: proportionally more.
IMAGE_TILE_LAYER = "dapi"

#: The pin. A floating version here would be the ``insitucnv`` mistake twice —
#: see the note in CLAUDE.md about why that git URL became a PyPI pin.
CELLDEGA_PIN = "celldega[pre]==0.24.2"


class DegaUnavailable(RuntimeError):
    """Celldega is not importable, with the reason and the remedy attached."""


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
        import celldega.pre  # noqa: F401

        # celldega imports pyvips in a try/except that leaves the module bound
        # to None, so a bare install reaches "generating dapi image tiles" and
        # dies there with an AttributeError. It is declared under celldega's own
        # `pre` extra; checking here turns a failure twenty minutes into an
        # export into one before it starts.
        import pyvips  # noqa: F401
    except ImportError as exc:
        return False, (
            f"Celldega is not installed ({exc}). Install it *without* its "
            f"dependency pins:\n"
            f"    pip install --no-deps {CELLDEGA_PIN}\n"
            f"A plain install downgrades anndata, spatialdata and pandas "
            f"underneath the viewer; celldega's upper caps are conservative and "
            f"it runs against the versions already here."
        )
    return True, ""


def require_dega() -> None:
    """Raise :class:`DegaUnavailable` with the remedy, or return."""
    ok, message = dega_available()
    if not ok:
        raise DegaUnavailable(message)


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
    import celldega.pre as dega_pre

    data_path = Path(data_path).resolve()
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
