"""
H&E image loading, landmark-based affine registration, and landmark I/O.

Provides utilities for aligning an H&E image to Xenium spatial data using
manually placed landmark points. The affine is estimated as a similarity
transform (rotation + uniform scale + translation) via scikit-image.

Coordinate convention:
  - napari uses (row, col) = (y, x)
  - skimage estimate_transform uses (x, y)
  - We convert between them via column-swap and matrix permutation.
"""

import json
from dataclasses import dataclass, field

import numpy as np
import dask.array as da
import tifffile
import cv2


#: Extra levels synthesised for a TIFF that carries no internal pyramid.
_SYNTHETIC_LEVELS = 4


def _aszarr_dask(store):
    """Wrap a tifffile ZarrTiffStore as a dask array, working around the
    fact that dask.array.from_zarr uses RegularChunkGrid, which was
    removed in zarr 3. Open the store with zarr first, then wrap with
    da.from_array — same end result, no internals access."""
    import zarr
    z = zarr.open(store, mode="r")
    chunks = getattr(z, "chunks", None) or "auto"
    return da.from_array(z, chunks=chunks)


def describe_pyramid(pyramid, label: str = "image") -> str:
    """One-line summary of what a loaded pyramid will cost to draw.

    ``load_he_pyramid`` can hand back either the file's own stored levels or
    levels it synthesised, and the memory behaviour is completely different
    between the two. Nothing said so before, which is why "the viewer died while
    loading an H&E" carried no information about the file that killed it.
    """
    from palms.utils.raster_io import level_is_computed

    base = pyramid[0]
    chained = any(level_is_computed(level) for level in pyramid)
    chunk = getattr(base, "chunksize", None)
    note = ""
    if chunk is not None and base.dtype.itemsize * int(np.prod(chunk)) > 512 << 20:
        # tifffile gives one chunk per plane for a file that is not tiled, so
        # nothing downstream can stream — the same shape of problem raster_io.py
        # documents for morphology_focus.
        note = " — NOT TILED, one chunk is the whole plane"
    return (
        f"{label}: {len(pyramid)} level(s), base {tuple(base.shape)} {base.dtype}, "
        f"chunks {chunk}, "
        f"levels {'CHAINED (materialise on draw)' if chained else 'stored in the file'}"
        f"{note}"
    )


def load_he_pyramid(path):
    """Load an H&E OME-TIFF/SVS as a list of dask arrays (one per pyramid level).

    Parameters
    ----------
    path : str or Path
        Path to the H&E image file.

    Returns
    -------
    pyramid : list of dask.array
        Pyramid levels sorted from highest to lowest resolution.
        Each array has shape (Y, X, C) for RGB images.
    tif : tifffile.TiffFile
        The open TiffFile object — must be kept alive to keep the zarr
        store valid for lazy reading.
    """
    tif = tifffile.TiffFile(str(path))
    store = tif.aszarr()

    # Determine number of pyramid levels
    n_levels = len(tif.series[0].levels)

    if n_levels > 1:
        # Image has an internal pyramid — read each level as dask
        pyramid = []
        for level_idx in range(n_levels):
            level_store = tif.aszarr(level=level_idx)
            arr = _aszarr_dask(level_store)
            pyramid.append(arr)
    else:
        # No internal pyramid — build one via 2x mean-pooling. Deliberately left
        # lazy: measured on a 16384x12288 RGB TIFF, materialising each level as
        # it is built cost *more* (tiled 0.83 -> 2.14 GB peak, strip-encoded
        # 2.63 -> 2.88 GB) because a coarsen over a tile-chunked base already
        # streams. Unlike morphology_focus's whole-page chunks, this chain is
        # not the expensive part — see describe_pyramid for what is.
        base = _aszarr_dask(store)
        pyramid = [base]
        current = base
        for _ in range(_SYNTHETIC_LEVELS):
            # Trim to even dimensions
            h, w = current.shape[0], current.shape[1]
            trimmed = current[:h - h % 2, :w - w % 2]
            # (Y, X, C) or (Y, X) — either way, coarsen the leading spatial dims.
            current = da.coarsen(
                np.mean, trimmed, {0: 2, 1: 2}, trim_excess=True,
            ).astype(base.dtype)
            pyramid.append(current)

    return pyramid, tif


def load_multichannel_pyramid(path):
    """Load a possibly-multichannel OME-TIFF as a dask pyramid with explicit channel axis.

    Unlike ``load_he_pyramid`` (which lets napari infer RGB), this probes
    ``tif.series[0].axes`` and OME-XML to locate the channel axis. Works for
    single-channel, RGB, and N-channel IF images.

    Parameters
    ----------
    path : str or Path
        Path to the OME-TIFF (or plain TIFF / SVS).

    Returns
    -------
    pyramid : list of dask.array
        Levels sorted from highest to lowest resolution. Dimensions match the
        source axes layout (e.g. (C, Y, X) or (Y, X, C) or (Y, X)).
    tif : tifffile.TiffFile
        Open file handle — keep alive for lazy reads.
    channel_axis : int | None
        Index of the channel axis within each pyramid level, or None if the
        image has no channel dimension (single-channel 2-D).
    channel_names : list[str]
        Human-readable channel labels. Parsed from OME-XML when possible,
        otherwise ``["C0", "C1", ...]``.
    """
    import re
    import xml.etree.ElementTree as ET

    tif = tifffile.TiffFile(str(path))
    series = tif.series[0]
    axes = (series.axes or "").upper()

    # Build pyramid (same as load_he_pyramid) — no rgb= flag so layout is preserved.
    n_levels = len(series.levels)
    if n_levels > 1:
        pyramid = [_aszarr_dask(tif.aszarr(level=i)) for i in range(n_levels)]
    else:
        base = _aszarr_dask(tif.aszarr())
        pyramid = [base]

    base_shape = pyramid[0].shape
    ndim = len(base_shape)

    # Locate channel axis from axes string ('C' for channels, 'S' for samples/RGB)
    channel_axis = None
    for marker in ("C", "S"):
        idx = axes.find(marker)
        if idx != -1 and idx < ndim:
            channel_axis = idx
            break

    # Fallback heuristics if axes string missing or ambiguous
    if channel_axis is None and ndim >= 3:
        # Trailing small dim ∈ {3,4} → treat as YXC RGB(A)
        if base_shape[-1] in (3, 4):
            channel_axis = ndim - 1
        else:
            # Assume leading axis is channels (CYX layout, typical for IF)
            channel_axis = 0

    n_channels = base_shape[channel_axis] if channel_axis is not None else 1

    # Parse channel names from OME-XML if present
    channel_names = [f"C{i}" for i in range(n_channels)]
    ome_xml = getattr(tif, "ome_metadata", None)
    if ome_xml:
        try:
            root = ET.fromstring(ome_xml)
            ns = re.match(r"\{.*\}", root.tag)
            ns = ns.group(0) if ns else ""
            pixels = root.find(f".//{ns}Pixels")
            if pixels is not None:
                parsed = []
                for ch in pixels.findall(f"{ns}Channel"):
                    nm = ch.get("Name") or ch.get("Fluor") or ch.get("ID")
                    if nm:
                        parsed.append(nm)
                if parsed and len(parsed) == n_channels:
                    channel_names = parsed
        except Exception:
            pass

    # If single-channel, drop channel_axis to None (napari handles as mono image)
    if n_channels == 1:
        channel_axis = None

    # Build extra pyramid levels when none are present (multichannel-safe coarsen)
    if n_levels == 1 and pyramid[0].ndim >= 2:
        current = pyramid[0]
        for _ in range(_SYNTHETIC_LEVELS):
            # Determine which axes are spatial (all non-channel axes of size > 1)
            axes_to_coarsen = {
                ax: 2 for ax in range(current.ndim)
                if ax != channel_axis and current.shape[ax] > 1
            }
            if not axes_to_coarsen:
                break
            # Trim to even along coarsened axes
            slicer = [slice(None)] * current.ndim
            for ax in axes_to_coarsen:
                s = current.shape[ax]
                slicer[ax] = slice(0, s - s % 2)
            trimmed = current[tuple(slicer)]
            try:
                current = da.coarsen(
                    np.mean, trimmed, axes_to_coarsen, trim_excess=True,
                ).astype(pyramid[0].dtype)
            except Exception:
                break
            pyramid.append(current)

    return pyramid, tif, channel_axis, channel_names


def parse_rgb_image_for_store(level):
    """Parse an H&E/ARMS base level into an ``Image2DModel`` element, lazily.

    Returns ``(parsed, (height, width))``.

    The eager version of this — ``np.asarray(pyramid[0])``, then a numpy
    ``.astype(np.uint8)`` copy — pulled the full-resolution slide into RAM twice
    on the Qt main thread before a byte reached disk, and then handed
    ``Image2DModel.parse`` a dense array whose four ``scale_factors`` levels
    spatialdata computes in a single ``da.compute``. Keeping it dask means the
    write streams chunk by chunk; the dims, scale factors and chunking below are
    unchanged, so what lands in the store is byte-identical to before.

    ``da.map_blocks`` detaches the graph from tifffile's ``ZarrTiffStore``, which
    has no ``.root`` for zarr v3 / spatialdata to introspect — the same guard
    ``adata_persistence.save_external_image_to_sdata`` already carries. It hides
    the store; it does not detach the *file*, so the ``TiffFile`` returned by
    ``load_he_pyramid`` must stay alive until the write completes. The tabs keep
    it in ``he_state["he_tif"]`` / ``arms_state["he_tif"]`` for exactly that.
    """
    from spatialdata.models import Image2DModel

    base = level
    if not isinstance(base, da.Array):
        base = da.asarray(base)
    base = da.map_blocks(lambda x: x, base, dtype=base.dtype)

    if base.ndim == 3 and base.shape[-1] in (3, 4):
        shape_yx = (base.shape[0], base.shape[1])
        base_cyx = da.transpose(base, (2, 0, 1))
    else:
        shape_yx = (base.shape[-2], base.shape[-1])
        base_cyx = base

    parsed = Image2DModel.parse(
        base_cyx.astype(np.uint8), dims=("c", "y", "x"),
        scale_factors=[2, 2, 2, 2], chunks=(3, 1024, 1024),
    )
    return parsed, shape_yx


def compute_landmark_affine(xenium_pts_yx, he_pts_yx):
    """Estimate a similarity affine from paired landmark points.

    Parameters
    ----------
    xenium_pts_yx : ndarray, shape (N, 2)
        Landmark positions in Xenium/napari space, (y, x) convention.
    he_pts_yx : ndarray, shape (N, 2)
        Corresponding landmark positions in H&E pixel space, (y, x) convention.

    Returns
    -------
    affine_3x3 : ndarray, shape (3, 3)
        Affine matrix in napari (y, x) convention, suitable for
        ``layer.affine = affine_3x3``.
    residuals : ndarray, shape (N,)
        Per-landmark Euclidean residual (in napari pixel units) after transform.
    """
    from skimage.transform import estimate_transform

    # Convert (y, x) → (x, y) for skimage
    src_xy = he_pts_yx[:, ::-1].astype(np.float64)
    dst_xy = xenium_pts_yx[:, ::-1].astype(np.float64)

    tform = estimate_transform('similarity', src=src_xy, dst=dst_xy)
    M_xy = tform.params  # 3x3 in (x, y) space

    # Permute to (y, x) convention: M_yx = P @ M_xy @ P
    P = np.array([[0, 1, 0],
                  [1, 0, 0],
                  [0, 0, 1]], dtype=np.float64)
    affine_3x3 = P @ M_xy @ P

    # Compute residuals in napari (y, x) space
    n = len(he_pts_yx)
    he_homo = np.hstack([he_pts_yx, np.ones((n, 1))])  # (N, 3)
    transformed = (affine_3x3 @ he_homo.T).T[:, :2]    # (N, 2) in (y, x)
    residuals = np.sqrt(np.sum((transformed - xenium_pts_yx) ** 2, axis=1))

    return affine_3x3, residuals


def save_landmarks(path, xenium_yx, he_yx, affine=None, he_filename=None,
                   flip_v=None, flip_h=None):
    """Save landmark points and optional affine to a JSON file.

    Parameters
    ----------
    path : str or Path
        Output JSON file path.
    xenium_yx : ndarray, shape (N, 2)
        Xenium landmark positions in (y, x) convention.
    he_yx : ndarray, shape (N, 2)
        H&E landmark positions in (y, x) convention.
    affine : ndarray or None
        3x3 affine matrix (napari convention).
    he_filename : str or None
        Name of the H&E file for reference.
    flip_v, flip_h : bool or None
        The orientation ``he_yx`` and ``affine`` are expressed in. The file stays
        self-consistent — ``affine`` maps ``he_yx`` onto ``xenium_yx`` — which is
        what any reader assumes, and these say how to put the points back on an
        unflipped layer. Omitted entirely when there is no flip, so a file
        written for an unflipped H&E is unchanged.
    """
    data = {
        "xenium_landmarks_yx": xenium_yx.tolist(),
        "he_landmarks_yx": he_yx.tolist(),
    }
    if flip_v:
        data["flip_v"] = bool(flip_v)
    if flip_h:
        data["flip_h"] = bool(flip_h)
    if affine is not None:
        data["affine_3x3_yx"] = affine.tolist()
    if he_filename is not None:
        data["he_filename"] = str(he_filename)

    with open(str(path), "w") as f:
        json.dump(data, f, indent=2)


def load_landmarks(path):
    """Load landmark points and optional affine from a JSON file.

    Parameters
    ----------
    path : str or Path
        Input JSON file path.

    Returns
    -------
    data : dict
        Keys: ``xenium_landmarks_yx``, ``he_landmarks_yx`` (as ndarrays),
        optionally ``affine_3x3_yx`` (as ndarray) and ``he_filename`` (str).
    """
    with open(str(path)) as f:
        raw = json.load(f)

    result = {
        "xenium_landmarks_yx": np.array(raw["xenium_landmarks_yx"], dtype=np.float64),
        "he_landmarks_yx": np.array(raw["he_landmarks_yx"], dtype=np.float64),
    }
    if "affine_3x3_yx" in raw:
        result["affine_3x3_yx"] = np.array(raw["affine_3x3_yx"], dtype=np.float64)
    if "he_filename" in raw:
        result["he_filename"] = raw["he_filename"]
    for key in ("flip_v", "flip_h"):
        if key in raw:
            result[key] = bool(raw[key])

    return result


# ─── Tissue masks ─────────────────────────────────────────────────────────────

def pick_level(pyramid, min_long_side=384):
    """Return the smallest pyramid level whose longest spatial side is >= *min_long_side*.

    Coarse alignment used to take ``pyramid[-1]``, which is not a size — it is
    however many levels the file happens to carry. For one and the same pancreas
    H&E that is 860x466 read from the OME-TIFF and 1718x931 read back from the
    zarr cache, so the same dataset was aligned at two different resolutions
    depending on whether the session had been restored.
    """
    best = pyramid[0]
    for level in pyramid:
        shape = tuple(level.shape)
        # (Y, X, C) or (C, Y, X) or (Y, X) — the two largest dims are the spatial ones.
        long_side = max(sorted(shape)[-2:])
        if long_side >= min_long_side and np.prod(shape) <= np.prod(best.shape):
            best = level
    return best


def extract_tissue_mask(image_gray, blur_ksize=5, open_ksize=5,
                        close_ksize=5, min_area_ratio=0.01,
                        method="otsu", outline=False, outline_sigma=6):
    """Extract a binary tissue mask via thresholding + morphological cleanup.

    Parameters
    ----------
    image_gray : ndarray, uint8, shape (H, W)
        Grayscale image (0–255).
    blur_ksize : int
        Median-blur kernel size (must be odd, 0 to skip).
    open_ksize : int
        Morphological-opening disk radius (removes noise).
    close_ksize : int
        Morphological-closing disk radius (fills holes).
    min_area_ratio : float
        Connected components smaller than this fraction of total area are removed.
    method : {"otsu", "triangle"}
        Which automatic threshold to use. **This choice is not cosmetic.** Otsu
        splits a histogram at the valley between its two largest modes, which is
        the right question for an H&E saturation image — unstained glass against
        stained tissue — and the wrong one for a fluorescence max-projection,
        where the two modes are *both tissue*. On the pancreas reference dataset
        Otsu cut between the bright acinar region and the dim fibrotic one and
        kept 17.4% of the frame; the triangle threshold, which anchors on the
        single dominant background peak, keeps 88.0% and is the actual section.
    outline : bool
        Reduce the mask to one filled outline: smooth, keep the largest connected
        component, fill its contour. A tissue *outline* is what coarse alignment
        compares; the speckle inside it is texture, and it makes the area — and
        therefore any scale estimated from the area — depend on the threshold.
    outline_sigma : float
        Gaussian sigma (mask pixels) used to smooth before tracing the outline.

    Returns
    -------
    mask : ndarray, uint8, shape (H, W)
        Binary mask with tissue=255, background=0.
    """
    # Median blur
    if blur_ksize > 0:
        ksize = blur_ksize if blur_ksize % 2 == 1 else blur_ksize + 1
        blurred = cv2.medianBlur(image_gray, ksize)
    else:
        blurred = image_gray

    if method == "triangle":
        flag = cv2.THRESH_TRIANGLE
    elif method == "otsu":
        flag = cv2.THRESH_OTSU
    else:
        raise ValueError(f"unknown threshold method {method!r}")
    _, mask = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY + flag)

    # Morphological opening (remove small noise)
    if open_ksize >= 1:
        kernel = cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE, (2 * open_ksize + 1, 2 * open_ksize + 1))
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)

    # Morphological closing (fill holes)
    if close_ksize >= 1:
        kernel = cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE, (2 * close_ksize + 1, 2 * close_ksize + 1))
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)

    # Remove small connected components
    if min_area_ratio > 0:
        total_area = mask.shape[0] * mask.shape[1]
        min_area = int(total_area * min_area_ratio)
        n_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask)
        for i in range(1, n_labels):  # skip background label 0
            if stats[i, cv2.CC_STAT_AREA] < min_area:
                mask[labels == i] = 0

    if outline:
        mask = _filled_outline(mask, outline_sigma)

    return mask


def _filled_outline(mask, sigma=6):
    """Smooth a binary mask, keep its largest component, and fill it."""
    if sigma > 0:
        mask = cv2.GaussianBlur(mask, (0, 0), sigma)
        mask = ((mask > 127).astype(np.uint8)) * 255
    n_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask)
    if n_labels > 1:
        keep = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
        mask = np.where(labels == keep, 255, 0).astype(np.uint8)
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    filled = np.zeros_like(mask)
    cv2.drawContours(filled, contours, -1, 255, -1)
    return filled


def extract_tissue_mask_fluorescence(image_cyx):
    """Build a tissue mask from a multi-channel fluorescence image.

    Uses the triangle threshold, not Otsu — see ``extract_tissue_mask``.

    Parameters
    ----------
    image_cyx : ndarray, shape (C, Y, X), uint16
        Multi-channel fluorescence (e.g. morphology_focus lowest pyramid level).

    Returns
    -------
    mask : ndarray, uint8, shape (Y, X)
    """
    proj = _normalise_u8(np.max(np.asarray(image_cyx), axis=0))
    return extract_tissue_mask(proj, method="triangle", outline=True)


def extract_tissue_mask_he(image_rgb):
    """Build a tissue mask from an H&E RGB image using HSV saturation.

    Tissue regions have high saturation (purple/pink stain) while background
    (white slide) has near-zero saturation — a genuinely bimodal histogram, so
    Otsu is the right threshold here.

    Parameters
    ----------
    image_rgb : ndarray, shape (Y, X, 3), uint8

    Returns
    -------
    mask : ndarray, uint8, shape (Y, X)
    """
    rgb = np.asarray(image_rgb)
    if rgb.ndim == 3 and rgb.shape[-1] == 4:
        rgb = rgb[..., :3]
    hsv = cv2.cvtColor(np.ascontiguousarray(rgb), cv2.COLOR_RGB2HSV)
    saturation = hsv[:, :, 1]  # uint8, 0–255
    return extract_tissue_mask(saturation, method="otsu", outline=True)


def _normalise_u8(image, lo_pct=1.0, hi_pct=99.5):
    """Percentile-stretch any 2-D array to uint8."""
    a = np.asarray(image).astype(np.float32)
    lo, hi = np.percentile(a, (lo_pct, hi_pct))
    return np.clip((a - lo) / (hi - lo + 1e-8) * 255, 0, 255).astype(np.uint8)


def mask_area_scale(target_mask, source_mask, target_downsample=1.0,
                    source_downsample=1.0):
    """Scale (source px -> target px) implied by two tissue-outline areas.

    Only meaningful once both outlines are correct — with the pre-2026-09 masks
    this was out by 2.4x on the pancreas dataset. With the corrected masks it is
    within 0.04% on a crop and 3.5% on a whole section, which makes it a usable
    prior when the H&E declares no pixel size of its own.
    """
    src_area = float((np.asarray(source_mask) > 0).sum())
    tgt_area = float((np.asarray(target_mask) > 0).sum())
    if src_area <= 0 or tgt_area <= 0:
        raise ValueError("cannot estimate scale — one tissue mask is empty")
    return np.sqrt(tgt_area / src_area) * target_downsample / source_downsample


# ─── Nuclear-density fields ───────────────────────────────────────────────────
#
# Coarse alignment scores *these*, not the binary masks. A tissue outline says
# nothing at all when the tissue fills both frames, which is exactly what a crop
# export looks like: on demo_data/crop_6 the corrected masks cover 100.0% and
# 99.6% of their images, and an outline search scores a perfect 1.0 at a scale
# 36% wrong. Blurred nuclear density still carries the internal structure, and
# recovers that dataset's 180 deg rotation.

def nuclear_density_fluorescence(image_cyx, channel=0):
    """Nuclear density from a morphology_focus thumbnail, as a 0-1 float field.

    Channel 0 of ``morphology_focus`` is DAPI. No blurring here: the smoothing
    that matters is applied by ``compute_coarse_affine`` at its own working
    resolution, which is the only place that knows how big a working pixel is.
    """
    arr = np.asarray(image_cyx)
    if arr.ndim == 3:
        arr = arr[min(channel, arr.shape[0] - 1)]
    return _normalise_u8(arr).astype(np.float32) / 255.0


def nuclear_density_he(image_rgb):
    """Haematoxylin density from an H&E thumbnail, as a 0-1 float field.

    ``1 - green`` is a serviceable haematoxylin proxy: haematoxylin absorbs most
    strongly in green, so the green channel is the one that darkens with nuclear
    density. Full colour deconvolution buys nothing here — the field is about to
    be reduced to a density at a few hundred pixels across anyway.
    """
    rgb = np.asarray(image_rgb)
    if rgb.ndim == 3 and rgb.shape[-1] >= 3:
        green = rgb[..., 1].astype(np.float32)
    else:
        green = np.squeeze(rgb).astype(np.float32)
    return 1.0 - green / 255.0


def build_alignment_fields(morph_cyx, he_rgb):
    """The pair of density fields coarse alignment compares."""
    return nuclear_density_fluorescence(morph_cyx), nuclear_density_he(he_rgb)


def he_pixel_size_um(tif):
    """Physical pixel size of an H&E TIFF in microns, or None if it declares none.

    Worth reaching for first: on the pancreas reference the OME-XML says
    0.27377 um, which against the dataset's 0.2125 um gives a scale of 1.28833
    where the truth is 1.28912 — 0.06%. No search over tissue outlines gets
    close to that, and it costs one XML attribute.

    Never raises; a file with no such metadata is the ordinary case.
    """
    import re
    import xml.etree.ElementTree as ET

    _TO_UM = {"µm": 1.0, "um": 1.0, "micron": 1.0, "microns": 1.0,
              "nm": 1e-3, "mm": 1e3, "cm": 1e4, "m": 1e6}
    try:
        ome = getattr(tif, "ome_metadata", None)
        if ome:
            root = ET.fromstring(ome)
            ns = re.match(r"\{.*\}", root.tag)
            ns = ns.group(0) if ns else ""
            pixels = root.find(f".//{ns}Pixels")
            if pixels is not None and pixels.get("PhysicalSizeX"):
                value = float(pixels.get("PhysicalSizeX"))
                unit = pixels.get("PhysicalSizeXUnit", "µm")
                factor = _TO_UM.get(unit, _TO_UM.get(unit.replace("μ", "µ"), None))
                if factor and value > 0:
                    return value * factor
    except Exception:
        pass
    try:
        tags = tif.pages[0].tags
        xres = tags["XResolution"].value
        unit = int(tags["ResolutionUnit"].value)
        pixels_per_unit = xres[0] / xres[1] if isinstance(xres, tuple) else float(xres)
        if pixels_per_unit <= 0:
            return None
        if unit == 2:      # inch
            return 25400.0 / pixels_per_unit
        if unit == 3:      # centimetre
            return 10000.0 / pixels_per_unit
    except Exception:
        pass
    return None


# ─── Coarse alignment ─────────────────────────────────────────────────────────

#: How far to search around a scale prior, by where the prior came from. A
#: declared pixel size is a measurement (0.06% on the pancreas reference); a
#: tissue-area ratio is an estimate (1.0% on a crop, 3.5% on a whole section).
#:
#: The scale is searched even when the prior is a measurement, because the two
#: errors trade off: holding the pancreas scale at its declared 1.2883 leaves the
#: rotation 0.87 deg out and 78 um of disagreement, where letting the search move
#: it to 1.2749 costs 1.10% of scale and gives 45 um. A band of 0 does hold it
#: fixed, for a caller that wants that.
SCALE_BAND_FOR = {"metadata": 0.05, "tissue-area": 0.25, "unknown": 0.75}

#: Below this NCC the match is not worth presenting as an answer, and it must
#: beat the best *materially different* orientation by MIN_COARSE_MARGIN, or the
#: search has found a tie rather than a match.
#:
#: Both floors are set from measurement, not taste. Images with nothing to match
#: score 0.24-0.34 with margins of 0.04-0.08; the two datasets with known
#: transforms score 0.49 and 0.78 with margins of 0.12 and 0.20. The floors sit
#: just under the tighter real case (crop_6, 0.494 / 0.121) rather than midway,
#: because the asymmetry runs one way: a low-confidence label on a good result
#: costs the user a glance at the overlay, while a confident label on a wrong
#: one sends them off to place landmarks against a transform that is not close.
#: Neither refuses the result — the transform is applied either way.
MIN_COARSE_SCORE = 0.40
MIN_COARSE_MARGIN = 0.10
#: Candidates whose warped source covers less of the target than this are skipped
#: — a transform that pushes the image off the canvas can otherwise score well on
#: the sliver that remains.
_MIN_VALID_FRACTION = 0.15
#: Two orientations closer than this are the same answer, not a rival one.
_DISTINCT_ANGLE_DEG = 20.0
#: Ceiling on how many hypotheses the global pass scores. A wider scale band
#: spends this on scales and coarsens the angle, rather than taking longer.
_GLOBAL_BUDGET = 2600
#: Longest side of the fields each pass runs on.
_SEARCH_LONG_SIDE = 256
_REFINE_LONG_SIDE = 512

#: Smoothing applied at the working resolution, in working pixels — which is the
#: only frame in which it can be stated once and be right everywhere.
#:
#: DAPI density and haematoxylin density agree about tissue *architecture*, not
#: about individual nuclei, so the comparison has to be made at a coarse spatial
#: scale. Blurring in thumbnail pixels instead made the amount of smoothing an
#: accident of which pyramid level the file happened to carry: at sigma 1 working
#: pixel the objective peaks within 0.12 deg of the pancreas truth at every
#: resolution tried (256, 512, 768), while at sigma 4 the peak sits 4 deg away.
_WORK_SIGMA = 1.0

_P_SWAP = np.array([[0, 1, 0],
                    [1, 0, 0],
                    [0, 0, 1]], dtype=np.float64)


@dataclass
class CoarseAlignResult:
    """What the coarse search found, and how much to believe it."""

    affine_3x3_yx: np.ndarray
    score: float
    runner_up_score: float
    scale: float
    rotation_deg: float
    mirrored: bool
    margin: float = float("inf")
    scale_source: str = "search"
    confident: bool = True
    notes: list = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"Match score: {self.score:.3f} "
            f"(next distinct orientation {self.runner_up_score:.3f}, "
            f"margin {self.margin:.3f})",
            f"Scale: {self.scale:.4f}  (from {self.scale_source})",
            f"Rotation: {self.rotation_deg:.2f} deg",
            f"Mirrored: {'yes' if self.mirrored else 'no'}",
        ]
        if not self.confident:
            lines.append(
                "LOW CONFIDENCE — the tissue may not be distinctive enough to "
                "align automatically. Check the overlay before placing landmarks.")
        lines.extend(self.notes)
        return "\n".join(lines)


def _rescale_field(field, long_side):
    """Downscale a float field so its longest side is ~*long_side*.

    Returns ``(small, fx, fy)`` where ``x_original = x_small * fx``. The two
    factors differ by a fraction of a percent, because the output size is
    rounded to whole pixels; that is why the final matrix is projected back onto
    an exact similarity rather than assumed to be one.
    """
    h, w = field.shape[:2]
    factor = max(h, w) / float(long_side)
    if factor <= 1.0:
        return np.ascontiguousarray(field, dtype=np.float32), 1.0, 1.0
    # Round to a size the FFT likes: the translation for every hypothesis comes
    # from cv2.phaseCorrelate, and a prime dimension is its worst case — 336x601
    # costs 9.0 ms a call against 6.3 ms for 336x600.
    size = (cv2.getOptimalDFTSize(max(int(round(w / factor)), 1)),
            cv2.getOptimalDFTSize(max(int(round(h / factor)), 1)))
    small = cv2.resize(field, size, interpolation=cv2.INTER_AREA)
    return (np.ascontiguousarray(small, dtype=np.float32),
            w / float(small.shape[1]), h / float(small.shape[0]))


def _nearest_similarity(M_xy, about_xy):
    """The closest scaled rotation (or reflection) to *M_xy*, fixing *about_xy*.

    The rescale factors above are anisotropic by up to ~0.3%, and a similarity
    is what every consumer of this matrix assumes: ``compute_landmark_affine``
    fits one, ``decompose`` reads a single scale off one, and the H&E layer is
    a rigid overlay. Projecting once at the end is cheaper than pretending.
    """
    A = np.asarray(M_xy, dtype=np.float64)[:2, :2]
    U, S, Vt = np.linalg.svd(A)
    R = U @ Vt
    scale = float(S.mean())
    A2 = scale * R
    c = np.asarray(about_xy, dtype=np.float64)
    t = np.asarray(M_xy, dtype=np.float64)[:2, 2]
    out = np.eye(3)
    out[:2, :2] = A2
    out[:2, 2] = A @ c + t - A2 @ c
    return out


def _ncc(a, b, valid):
    a = a[valid]
    b = b[valid]
    if a.size < 100:
        return -1.0
    a = a - a.mean()
    b = b - b.mean()
    denom = np.sqrt(float((a * a).sum()) * float((b * b).sum()))
    return float((a * b).sum() / denom) if denom > 0 else -1.0


class _Candidate:
    """One (scale, rotation, mirror) hypothesis, scored on a pair of fields."""

    __slots__ = ("target", "source", "ones", "shape", "src_centre", "tgt_centre", "mirror_xy")

    def __init__(self, target, source):
        self.target = target
        self.source = source
        self.ones = np.ones_like(source, dtype=np.float32)
        self.shape = target.shape
        sh, sw = source.shape
        h, w = target.shape
        self.src_centre = (sw / 2.0, sh / 2.0)
        self.tgt_centre = (w / 2.0, h / 2.0)
        self.mirror_xy = np.array([[-1.0, 0.0, sw - 1.0],
                                   [0.0, 1.0, 0.0],
                                   [0.0, 0.0, 1.0]])

    def matrix(self, scale, deg, mirror):
        t = np.radians(deg)
        c, s = scale * np.cos(t), scale * np.sin(t)
        M = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        if mirror:
            M = M @ self.mirror_xy
        scx, scy = self.src_centre
        cx = M[0, 0] * scx + M[0, 1] * scy + M[0, 2]
        cy = M[1, 0] * scx + M[1, 1] * scy + M[1, 2]
        M[0, 2] += self.tgt_centre[0] - cx
        M[1, 2] += self.tgt_centre[1] - cy
        return M

    def _warp(self, M):
        h, w = self.shape
        warped = cv2.warpAffine(self.source, M[:2, :], (w, h), flags=cv2.INTER_LINEAR)
        valid = cv2.warpAffine(self.ones, M[:2, :], (w, h), flags=cv2.INTER_NEAREST) > 0.5
        return warped, valid

    def score(self, M, align_translation=True):
        """Score *M*, optionally letting phase correlation choose the translation."""
        warped, valid = self._warp(M)
        if valid.mean() < _MIN_VALID_FRACTION:
            return -1.0, M
        if align_translation:
            shift, _ = cv2.phaseCorrelate(
                np.ascontiguousarray(self.target * valid),
                np.ascontiguousarray(warped * valid))
            M = M.copy()
            M[0, 2] -= shift[0]
            M[1, 2] -= shift[1]
            warped, valid = self._warp(M)
            if valid.mean() < _MIN_VALID_FRACTION:
                return -1.0, M
        return _ncc(self.target, warped, valid), M


def _search(cand, scales, degrees, mirrors):
    """Score a (scale x rotation x mirror) grid. Returns a list of (score, scale, deg, mirror, M)."""
    out = []
    for mirror in mirrors:
        for scale in scales:
            for deg in degrees:
                score, M = cand.score(cand.matrix(scale, deg, mirror))
                if score > -1.0:
                    out.append((score, float(scale), float(deg), bool(mirror), M))
    out.sort(key=lambda r: -r[0])
    return out


def _polish(cand, scale0, deg0, mirror, step, spacing, rounds=3):
    """Coordinate descent on (rotation, scale) from one grid hypothesis.

    Fixed refinement windows are the trap here: with a 25% band over 7 scales the
    grid steps by 7.9%, so a +/-3% window around a grid point cannot reach the
    truth *even when the grid straddles it* — the pancreas section sat at exactly
    the prior, 3.5% out, with the answer permanently out of reach. Alternating
    the two axes and halving the window each round has no such blind spot.
    """
    best_score, best_M = cand.score(cand.matrix(scale0, deg0, mirror))
    best_scale, best_deg = scale0, deg0
    for _ in range(rounds):
        for deg in np.arange(best_deg - 1.5 * step, best_deg + 1.5 * step + 1e-9, step / 8.0):
            score, M = cand.score(cand.matrix(best_scale, deg, mirror))
            if score > best_score:
                best_score, best_deg, best_M = score, float(deg), M
        for scale in np.geomspace(best_scale / spacing, best_scale * spacing, 9):
            score, M = cand.score(cand.matrix(scale, best_deg, mirror))
            if score > best_score:
                best_score, best_scale, best_M = score, float(scale), M
        step, spacing = step / 2.0, spacing ** 0.5
    return best_score, best_scale, best_deg, mirror, best_M


def _angle_gap(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def _distinct_top(ranked, k, gap=_DISTINCT_ANGLE_DEG):
    """The best *k* hypotheses that are not each other's neighbours."""
    kept = []
    for row in ranked:
        if all(row[3] != o[3] or _angle_gap(row[2], o[2]) >= gap for o in kept):
            kept.append(row)
        if len(kept) >= k:
            break
    return kept


def compute_coarse_affine(target_field, source_field,
                          target_downsample=1.0, source_downsample=1.0,
                          scale_prior=None, scale_band=0.25, n_scales=7,
                          rotation_step=2.0, mirror=True, refine=True,
                          scale_source="search", source_shape_yx=None):
    """Find the similarity that best overlays an H&E on the Xenium morphology.

    Scores the *normalised cross-correlation of blurred nuclear density*, over a
    global search in rotation (and, optionally, reflection), with the translation
    for each hypothesis taken from phase correlation. It replaces a moment-based
    estimate — centroid, area ratio, principal axis, four 90-degree hypotheses —
    which on the pancreas reference dataset returned scale 0.5515 against a truth
    of 1.2891 and -55.0 deg against -89.88 deg. Two things were wrong with it and
    only one was the search: the fluorescence tissue mask kept a third of the
    section, so the *true* transform scored 0.216 where the wrong answer scored
    0.493. A better search over a broken objective finds a wrong answer faster.

    Density rather than a tissue outline, because an outline says nothing when
    tissue fills both images — which is what every crop export looks like. On
    ``demo_data/crop_6`` the corrected masks cover 100.0% and 99.6% of their
    frames and an outline search scores a perfect 1.0 at a scale 36% wrong. The
    outlines are still used, for the scale prior.

    Measured 2026-09-02 against the two datasets that carry an independent
    transform: scale within 0.7% and rotation within 0.04 deg on both, which is
    15 um and 17 um of mean disagreement across the whole slide. See
    ``scripts/score_coarse_align.py``, which is how those numbers are produced.

    Parameters
    ----------
    target_field, source_field : ndarray, float32, 2-D
        Nuclear-density fields, from ``nuclear_density_fluorescence`` and
        ``nuclear_density_he``, at thumbnail resolution.
    target_downsample, source_downsample : float
        full_res / thumbnail ratio for each, so the result applies to full-res
        pixels.
    scale_prior : float or None
        Expected full-resolution scale (H&E px -> Xenium px). ``he_pixel_size_um``
        divided by the dataset pixel size is the best source; ``mask_area_scale``
        is the fallback. ``None`` anchors the search at 1.0 in *thumbnail* space
        and widens the band, which is a guess and is reported as one.
    scale_band : float
        Fractional half-width of the geometric scale search around the prior.
        **Zero holds the scale fixed at the prior.** ``SCALE_BAND_FOR`` maps a
        ``scale_source`` to the band that was measured to work for it.
    n_scales : int
        Scales in the global grid. With ``rotation_step`` this sets the size of
        the joint grid, which ``_GLOBAL_BUDGET`` then caps by coarsening the
        angle — a wider scale band costs angular resolution, not wall-clock.
    rotation_step : float
        Degrees between hypotheses in the global pass. It is a floor: the budget
        may enlarge it.
    mirror : bool
        Search reflections too. A mirrored H&E cannot be corrected by
        ``compute_landmark_affine`` (a similarity has no reflection), so finding
        it here is the only way it is ever found.
    refine : bool
        Run the coordinate-descent polish. Without it the answer is only as good
        as the global grid, which is coarse by design.
    source_shape_yx : tuple or None
        Full-resolution H&E shape. Used to factor a detected reflection out of
        the returned matrix, so the caller can carry it as an explicit flip; see
        the note on the return value.

    Returns
    -------
    CoarseAlignResult
        ``affine_3x3_yx`` maps full-resolution H&E pixels to full-resolution
        Xenium pixels, in napari (y, x). **When ``mirrored`` is true it maps the
        horizontally flipped H&E**, so the caller must compose it with the same
        flip — which is what ``tab_he_registration`` does, and what lets the
        landmark refinement inherit the reflection instead of discarding it.
    """
    target = np.ascontiguousarray(np.asarray(target_field, dtype=np.float32))
    source = np.ascontiguousarray(np.asarray(source_field, dtype=np.float32))
    if target.ndim != 2 or source.ndim != 2:
        raise ValueError("coarse alignment needs two 2-D density fields")

    notes = []
    # The prior is stated in full-resolution terms; the search works in thumbnails.
    if scale_prior is None:
        thumb_prior = 1.0
        scale_band = max(scale_band, 0.75)
        scale_source = "unknown"
        notes.append("No scale prior was available; the scale was searched blind.")
    else:
        thumb_prior = float(scale_prior) * source_downsample / target_downsample

    mirrors = (False, True) if mirror else (False,)
    # A band of zero holds the scale at the prior. See SCALE_BAND_FOR for why
    # that is not the default even when the prior is a measured pixel size.
    fixed_scale = scale_band <= 0

    def _tier(long_side):
        """Both fields at a working resolution, smoothed by the same amount.

        Both are reduced to the same *longest side*, not to the same physical
        pixel size. Matching physical pixel size is the more principled-sounding
        choice and it measures worse — on the pancreas section it scores 0.552
        and puts the optimum 3.5 deg off the truth, against 0.811 and 0.1 deg for
        equal long sides. Equal long sides is defensible on its own terms: an H&E
        and the morphology image of the same section cover comparable ground, so
        the same pixel count over it is the same physical sampling either way.
        """
        t, ftx, fty = _rescale_field(target, long_side)
        sfield, fsx, fsy = _rescale_field(source, long_side)
        return (_Candidate(cv2.GaussianBlur(t, (0, 0), _WORK_SIGMA),
                           cv2.GaussianBlur(sfield, (0, 0), _WORK_SIGMA)),
                (ftx, fty), (fsx, fsy), np.sqrt(fsx * fsy) / np.sqrt(ftx * fty))

    thumb_scales = (np.array([thumb_prior]) if fixed_scale else
                    np.geomspace(thumb_prior / (1 + scale_band),
                                 thumb_prior * (1 + scale_band), n_scales))
    # Keep the joint grid inside a fixed budget by coarsening the angle, so a
    # wider scale band costs angular resolution rather than wall-clock.
    step = max(rotation_step, 360.0 * len(thumb_scales) * len(mirrors) / _GLOBAL_BUDGET)

    # --- global pass over the *joint* (scale, rotation) grid ---
    # Not rotation alone at a bracketed scale: the two are coupled, and a scale
    # prior that is 3.4% out moves the best rotation by more than 3 degrees. A
    # rotation-first search then hands the polish a seed it cannot walk back
    # from — measured at 214 um of disagreement on the pancreas section, against
    # 45 um for the same code given a prior that happened to be good.
    cand, tfac, sfac, to_work = _tier(_SEARCH_LONG_SIDE)
    ranked = _search(cand, thumb_scales * to_work, np.arange(-180.0, 180.0, step), mirrors)
    if not ranked:
        raise ValueError("coarse alignment found no usable overlap between the images")

    # How far the winner beats the best *materially different* orientation, read
    # in the global pass — the only pass that looked at every orientation. It is
    # the ambiguity measure: on a pair of images with nothing to match, every
    # orientation scores alike and this collapses.
    best = ranked[0]
    runner_up = next((r[0] for r in ranked[1:]
                      if r[3] != best[3] or _angle_gap(r[2], best[2]) >= _DISTINCT_ANGLE_DEG),
                     -1.0)
    margin = best[0] - runner_up if runner_up > -1.0 else float("inf")
    seeds = _distinct_top(ranked, 5)
    score, scale_work, deg_work, mir_work, M_work = best

    # --- polish each surviving orientation, then the winner at finer resolution ---
    if refine:
        # How far apart the global grid's scales are; the polish must be able to
        # cross that gap, not merely wobble inside a cell.
        spacing = 1.0 if fixed_scale or len(thumb_scales) < 2 else float(
            (thumb_scales[-1] / thumb_scales[0]) ** (1.0 / (len(thumb_scales) - 1)))
        polished = [_polish(cand, sc, deg, mir, step, spacing)
                    for _, sc, deg, mir, _ in seeds]
        polished.sort(key=lambda r: -r[0])
        score, scale_work, deg_work, mir_work, M_work = polished[0]

        cand2, tfac2, sfac2, to_work2 = _tier(_REFINE_LONG_SIDE)
        if cand2.target.shape != cand.target.shape:
            score, scale_work, deg_work, mir_work, M_work = _polish(
                cand2, scale_work / to_work * to_work2, deg_work, mir_work,
                step / 4.0, spacing ** 0.25, rounds=2)
            tfac, sfac, to_work = tfac2, sfac2, to_work2

    # --- working space -> thumbnail space -> full resolution -> napari (y, x) ---
    M_thumb_xy = (np.diag([tfac[0], tfac[1], 1.0]) @ M_work
                  @ np.diag([1.0 / sfac[0], 1.0 / sfac[1], 1.0]))
    M_full_xy = (np.diag([target_downsample, target_downsample, 1.0]) @ M_thumb_xy
                 @ np.diag([1.0 / source_downsample, 1.0 / source_downsample, 1.0]))
    M_full_xy = _nearest_similarity(
        M_full_xy,
        (source.shape[1] * source_downsample / 2.0,
         source.shape[0] * source_downsample / 2.0))

    if mir_work:
        # Hand back the transform of the *flipped* image, exactly: the caller
        # composes it with the same full-resolution flip. Leaving the reflection
        # inside the matrix would work until the user pressed Compute
        # Registration, which fits a similarity and would silently drop it.
        if source_shape_yx is None:
            source_shape_yx = (source.shape[0] * source_downsample,
                               source.shape[1] * source_downsample)
        full_w = float(source_shape_yx[1])
        flip_full_xy = np.array([[-1.0, 0.0, full_w - 1.0],
                                 [0.0, 1.0, 0.0],
                                 [0.0, 0.0, 1.0]])
        M_full_xy = M_full_xy @ np.linalg.inv(flip_full_xy)

    affine_3x3_yx = _P_SWAP @ M_full_xy @ _P_SWAP

    full_scale = float(np.hypot(affine_3x3_yx[0, 0], affine_3x3_yx[0, 1]))
    rotation = float(np.degrees(np.arctan2(affine_3x3_yx[0, 1], affine_3x3_yx[0, 0])))
    confident = bool(score >= MIN_COARSE_SCORE and margin >= MIN_COARSE_MARGIN)

    return CoarseAlignResult(
        affine_3x3_yx=affine_3x3_yx,
        score=float(score),
        runner_up_score=float(runner_up),
        scale=full_scale,
        rotation_deg=rotation,
        mirrored=bool(mir_work),
        margin=float(margin),
        scale_source=scale_source,
        confident=confident,
        notes=notes,
    )
