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
import numpy as np
import dask.array as da
import tifffile
import cv2


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
            arr = da.from_zarr(level_store)
            pyramid.append(arr)
    else:
        # No internal pyramid — build one via 2x mean-pooling
        base = da.from_zarr(store)
        pyramid = [base]
        current = base
        for _ in range(4):  # build 4 additional levels
            # Trim to even dimensions
            h, w = current.shape[0], current.shape[1]
            trimmed = current[:h - h % 2, :w - w % 2]
            if current.ndim == 3:
                # (Y, X, C) — coarsen spatial dims only
                current = da.coarsen(
                    np.mean, trimmed, {0: 2, 1: 2}, trim_excess=True,
                ).astype(current.dtype)
            else:
                current = da.coarsen(
                    np.mean, trimmed, {0: 2, 1: 2}, trim_excess=True,
                ).astype(current.dtype)
            pyramid.append(current)

    return pyramid, tif


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


def save_landmarks(path, xenium_yx, he_yx, affine=None, he_filename=None):
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
    """
    data = {
        "xenium_landmarks_yx": xenium_yx.tolist(),
        "he_landmarks_yx": he_yx.tolist(),
    }
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

    return result


# ─── Coarse tissue-outline alignment (OpenCV moments) ─────────────────────────

def extract_tissue_mask(image_gray, blur_ksize=5, open_ksize=5,
                        close_ksize=5, min_area_ratio=0.01):
    """Extract a binary tissue mask via Otsu thresholding + morphological cleanup.

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

    # Otsu threshold
    _, mask = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

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

    return mask


def extract_tissue_mask_fluorescence(image_cyx):
    """Build a tissue mask from a multi-channel fluorescence image.

    Parameters
    ----------
    image_cyx : ndarray, shape (C, Y, X), uint16
        Multi-channel fluorescence (e.g. morphology_focus lowest pyramid level).

    Returns
    -------
    mask : ndarray, uint8, shape (Y, X)
    """
    # Max-project across channels → single 2-D image
    proj = np.max(image_cyx, axis=0)  # (Y, X), uint16

    # Normalize to uint8
    pmin, pmax = np.percentile(proj, (1, 99.5))
    proj = np.clip((proj.astype(np.float32) - pmin) / (pmax - pmin + 1e-8) * 255,
                   0, 255).astype(np.uint8)

    return extract_tissue_mask(proj)


def extract_tissue_mask_he(image_rgb):
    """Build a tissue mask from an H&E RGB image using HSV saturation.

    Tissue regions have high saturation (purple/pink stain) while background
    (white slide) has near-zero saturation.

    Parameters
    ----------
    image_rgb : ndarray, shape (Y, X, 3), uint8

    Returns
    -------
    mask : ndarray, uint8, shape (Y, X)
    """
    hsv = cv2.cvtColor(image_rgb, cv2.COLOR_RGB2HSV)
    saturation = hsv[:, :, 1]  # uint8, 0–255
    return extract_tissue_mask(saturation)


def compute_coarse_affine(target_mask, source_mask,
                          target_downsample=1.0, source_downsample=1.0):
    """Compute a coarse similarity transform aligning tissue outlines.

    Uses image moments (centroid, area, principal-axis angle) with
    multi-rotation hypothesis testing (0/90/180/270 deg offsets) scored by IoU.

    The returned 3×3 matrix is in napari (y, x) convention and accounts for the
    downsample factors so it can be applied directly to the full-resolution layer.

    Parameters
    ----------
    target_mask : ndarray, uint8, (H, W)
        Binary tissue mask of the target (morphology_focus thumbnail).
    source_mask : ndarray, uint8, (H, W)
        Binary tissue mask of the source (H&E thumbnail).
    target_downsample : float
        full_res / thumbnail ratio for target (e.g. 32 if thumbnail is 32× smaller).
    source_downsample : float
        full_res / thumbnail ratio for source.

    Returns
    -------
    affine_3x3_yx : ndarray, shape (3, 3)
        Affine in napari (row, col) = (y, x) convention, mapping full-res H&E
        pixel coords → full-res Xenium pixel coords.
    """
    # --- moment extraction helper ---
    def _moments_info(mask):
        binary = (mask > 0).astype(np.uint8)
        M = cv2.moments(binary)
        if M['m00'] == 0:
            return None
        cx = M['m10'] / M['m00']
        cy = M['m01'] / M['m00']
        mu20 = M['mu20'] / M['m00']
        mu02 = M['mu02'] / M['m00']
        mu11 = M['mu11'] / M['m00']
        theta = 0.5 * np.arctan2(2 * mu11, mu20 - mu02)
        area = M['m00']
        return {'cx': cx, 'cy': cy, 'theta': theta, 'area': area}

    tgt = _moments_info(target_mask)
    src = _moments_info(source_mask)
    if tgt is None or src is None:
        raise ValueError("Could not compute moments — one mask is empty")

    # Base similarity parameters (in thumbnail pixel space, OpenCV x/y convention)
    scale = np.sqrt(tgt['area'] / (src['area'] + 1e-8))
    base_theta = tgt['theta'] - src['theta']

    target_h, target_w = target_mask.shape[:2]
    target_binary = (target_mask > 0).astype(np.uint8)

    # Test 4 rotation hypotheses
    best_affine_2x3 = None
    best_iou = -1.0
    source_binary = (source_mask > 0).astype(np.uint8)

    for offset in [0, np.pi / 2, np.pi, 3 * np.pi / 2]:
        theta = base_theta + offset
        cos_t, sin_t = np.cos(theta), np.sin(theta)

        tx = tgt['cx'] - scale * (cos_t * src['cx'] - sin_t * src['cy'])
        ty = tgt['cy'] - scale * (sin_t * src['cx'] + cos_t * src['cy'])

        A = np.array([
            [scale * cos_t, -scale * sin_t, tx],
            [scale * sin_t,  scale * cos_t, ty],
        ], dtype=np.float64)

        warped = cv2.warpAffine(source_binary * 255, A, (target_w, target_h),
                                flags=cv2.INTER_NEAREST)

        intersection = np.sum((warped > 0) & (target_binary > 0))
        union = np.sum((warped > 0) | (target_binary > 0))
        iou = intersection / (union + 1e-8)

        if iou > best_iou:
            best_iou = iou
            best_affine_2x3 = A

    print(f"  Coarse alignment: best IoU = {best_iou:.4f}, "
          f"scale = {scale:.4f}, "
          f"rotation = {np.degrees(np.arctan2(best_affine_2x3[1, 0], best_affine_2x3[0, 0])):.1f} deg")

    # --- Convert thumbnail-space 2×3 (x,y) → full-res 3×3 (y,x) for napari ---
    # 1. Promote to 3×3 in (x,y) space
    M_xy = np.eye(3, dtype=np.float64)
    M_xy[:2, :] = best_affine_2x3

    # 2. Scale from thumbnail to full-res:
    #    A_full_xy = S_target @ M_xy @ inv(S_source)
    #    where S = diag(ds, ds, 1)
    S_tgt = np.diag([target_downsample, target_downsample, 1.0])
    S_src_inv = np.diag([1.0 / source_downsample, 1.0 / source_downsample, 1.0])
    M_full_xy = S_tgt @ M_xy @ S_src_inv

    # 3. Permute (x,y) → (y,x) for napari: M_yx = P @ M_xy @ P
    P = np.array([[0, 1, 0],
                  [1, 0, 0],
                  [0, 0, 1]], dtype=np.float64)
    affine_3x3_yx = P @ M_full_xy @ P

    return affine_3x3_yx
