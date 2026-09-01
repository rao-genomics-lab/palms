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


def _merged_polygon_of_type(annotation_layer, annotation_type: str, pixel_size: float):
    """Union of the polygons of one annotation type, in xy-micron coordinates.

    Returns ``None`` if the layer holds no polygon of that type.

    ``make_valid``, not ``buffer(0)``: on a self-intersecting polygon — a bowtie
    is easy to draw by hand in napari — ``buffer(0)`` silently *deletes* a lobe
    rather than repairing it, so the caller would sample virtual cells from, or
    measure distances to, half the region the user drew. Measured on shapely
    2.1.2: a 10x10 bowtie comes back as one 25 µm² triangle from ``buffer(0)``
    and as the full 50 µm² MultiPolygon from ``make_valid``. Same rule, and the
    same reason, as ``roi.polygons.tmpl``.

    Only the *polygonal* parts of the repair are kept. ``make_valid`` answers a
    degenerate ring — a shape whose vertices are collinear, i.e. a line the user
    drew with the polygon tool — with a LineString, which has bounds and a
    boundary and would therefore sail on through both callers as if it enclosed
    a region. ``buffer(0)`` returned an empty polygon there, and that shape being
    ignored is the behaviour worth keeping.
    """
    from shapely import make_valid
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
        xy_um = np.asarray(arr)[:, ::-1] * pixel_size  # yx-px → xy-µm
        poly = Polygon(xy_um)
        if not poly.is_valid:
            poly = make_valid(poly)
        parts = getattr(poly, "geoms", [poly]) if poly.geom_type == "GeometryCollection" else [poly]
        polys += [g for g in parts if g.geom_type in ("Polygon", "MultiPolygon") and not g.is_empty]

    return unary_union(polys) if polys else None


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
    from shapely import contains_xy

    merged = _merged_polygon_of_type(annotation_layer, annotation_type, pixel_size)
    if merged is None:
        return np.empty((0, 2), dtype=np.float64)

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

    n = len(centroids_um_xy)
    merged = _merged_polygon_of_type(annotation_layer, annotation_type, pixel_size)
    if merged is None:
        return np.full(n, np.nan, dtype=np.float64)

    boundary = merged.boundary  # MultiLineString / LineString

    # Vectorised shapely distance (shapely >= 2.0)
    points = shapely.points(centroids_um_xy[:, 0], centroids_um_xy[:, 1])
    distances = shapely.distance(points, boundary)
    return np.asarray(distances, dtype=np.float64)
