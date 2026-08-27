"""Utilities for manual tissue annotation analysis.

Provides:
  - sample_annotation_centroids: grid-sample virtual cell positions inside annotation polygons
  - compute_distance_to_annotation: minimum distance from each cell to an annotation boundary
  - get_annotation_types: list unique annotation type names in the annotation layer
"""

from __future__ import annotations

import numpy as np


def get_annotation_types(annotation_layer) -> list[str]:
    """Return sorted list of unique non-empty annotation type strings in the layer."""
    if annotation_layer is None:
        return []
    types = annotation_layer.properties.get("annotation_type", [])
    return sorted({t for t in types if t and str(t).strip()})


def sample_annotation_centroids(
    annotation_layer,
    annotation_type: str,
    pixel_size: float,
    density_um2: float = 100.0,
) -> np.ndarray:
    """Sample virtual cell centroids inside polygons of the given annotation type.

    Parameters
    ----------
    annotation_layer : napari Shapes layer
        The annotation layer (shapes in yx-pixel coordinates).
    annotation_type : str
        The type to sample (must match entries in layer.properties["annotation_type"]).
    pixel_size : float
        Microns per pixel — used to convert polygon coordinates to microns.
    density_um2 : float
        One virtual cell per this many µm². Default 100 µm²/cell.

    Returns
    -------
    np.ndarray, shape (N, 2), columns = (x, y) in microns.
        Empty array if no polygons of the requested type exist.
    """
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    from shapely import contains_xy

    shapes = annotation_layer.data
    types = list(annotation_layer.properties.get("annotation_type", []))
    while len(types) < len(shapes):
        types.append("")

    # Collect polygons of this type, converting yx-pixels → xy-microns
    polys = []
    for arr, t in zip(shapes, types):
        if str(t) != annotation_type:
            continue
        xy_um = arr[:, ::-1] * pixel_size  # yx-px → xy-µm
        poly = Polygon(xy_um)
        if not poly.is_valid:
            poly = poly.buffer(0)
        if not poly.is_empty:
            polys.append(poly)

    if not polys:
        return np.empty((0, 2), dtype=np.float64)

    merged = unary_union(polys)
    minx, miny, maxx, maxy = merged.bounds

    # Spacing between grid points (µm) such that grid density ≈ 1 cell per density_um2
    step = float(np.sqrt(density_um2))

    xs = np.arange(minx + step / 2, maxx, step)
    ys = np.arange(miny + step / 2, maxy, step)
    if len(xs) == 0 or len(ys) == 0:
        return np.empty((0, 2), dtype=np.float64)

    gx, gy = np.meshgrid(xs, ys)
    grid = np.column_stack([gx.ravel(), gy.ravel()])  # (M, 2) xy-µm

    mask = contains_xy(merged, grid[:, 0], grid[:, 1])
    return grid[mask]


def compute_distance_to_annotation(
    centroids_um_xy: np.ndarray,
    annotation_layer,
    annotation_type: str,
    pixel_size: float,
) -> np.ndarray:
    """Compute minimum distance (µm) from each cell centroid to the annotation boundary.

    Uses the shapely vectorised API (shapely >= 2.0) for performance.

    Parameters
    ----------
    centroids_um_xy : np.ndarray, shape (N, 2)
        Cell centroids in xy-micron coordinates (from adata.obsm['spatial']).
    annotation_layer : napari Shapes layer
    annotation_type : str
        The annotation type whose boundary to measure against.
    pixel_size : float
        Microns per pixel.

    Returns
    -------
    np.ndarray of float64, length N.
        Minimum distance in microns to the nearest point on the annotation boundary.
        NaN for any cell if no polygons of the given type exist.
    """
    import shapely
    from shapely.geometry import Polygon
    from shapely.ops import unary_union

    shapes = annotation_layer.data
    types = list(annotation_layer.properties.get("annotation_type", []))
    while len(types) < len(shapes):
        types.append("")

    polys = []
    for arr, t in zip(shapes, types):
        if str(t) != annotation_type:
            continue
        xy_um = arr[:, ::-1] * pixel_size
        poly = Polygon(xy_um)
        if not poly.is_valid:
            poly = poly.buffer(0)
        if not poly.is_empty:
            polys.append(poly)

    n = len(centroids_um_xy)
    if not polys:
        return np.full(n, np.nan, dtype=np.float64)

    merged = unary_union(polys)
    boundary = merged.boundary  # MultiLineString / LineString

    # Vectorised shapely distance (shapely >= 2.0)
    points = shapely.points(centroids_um_xy[:, 0], centroids_um_xy[:, 1])
    distances = shapely.distance(points, boundary)
    return np.asarray(distances, dtype=np.float64)
