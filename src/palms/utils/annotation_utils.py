"""Utilities for manual tissue annotation analysis.

Provides:
  - get_annotation_types: list unique annotation type names in the annotation layer

The geometry used to live here too — ``sample_annotation_centroids`` and
``compute_distance_to_annotation``, which the two annotation tabs called
directly. Both are now *templates* (``annot.virtual_cells``, ``annot.distance``)
run through ``ctx.run_step``, so the code that computes a result is the code the
notebook records. Keeping a second python implementation of the same geometry
beside them would be the drift the Step system exists to remove: nothing would
have called it, and nothing would have noticed it disagreeing.
"""

from __future__ import annotations


def get_annotation_types(annotation_layer) -> list[str]:
    """Return sorted list of unique non-empty annotation type strings in the layer."""
    if annotation_layer is None:
        return []
    types = annotation_layer.properties.get("annotation_type", [])
    return sorted({t for t in types if t and str(t).strip()})
