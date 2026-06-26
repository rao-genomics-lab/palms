# Annotations

Draw and manage labelled tissue annotation polygons on the Annotations layer, assign annotation types, customise colours per type, and import or export annotations as GeoJSON. This tab is in the "Tools" control panel group.

<!-- SCREENSHOT: docs/screenshots/tab-annotations.png -->

## Controls

### Drawing

Annotations are drawn directly on the napari canvas, not through this tab's buttons:

1. Select the "Annotations" Shapes layer in the napari layer list.
2. Choose the polygon tool in the napari toolbar.
3. Draw a polygon outline around a tissue region of interest.
4. Press Enter to close each polygon.

### Management

| Control | Description |
|---|---|
| Annotation type: | Text field. Type the label to assign to selected shapes (for example, "bone", "adipocyte", or "vessel"). |
| Assign to selected shapes | Applies the entered type label to all shapes currently selected on the Annotations layer. |
| Annotation type table | Table with columns Type, Count, and Colour. Shows all annotation types with their shape counts. Click a cell in the Colour column to open a colour picker for that type. |
| Delete selected shapes | Removes the currently selected shapes from the layer. |
| Clear all annotations | Removes all shapes from the layer. Shows a confirmation dialog before proceeding. |
| Import GeoJSON... | Loads a GeoJSON FeatureCollection (polygon or multipolygon geometries with an `annotation_type` or `name` property) and appends the shapes to existing annotations. |
| Export GeoJSON... | Saves all current annotations as a GeoJSON FeatureCollection. |

## Workflow

1. Select the "Annotations" Shapes layer in the napari layer list.
2. Activate the polygon tool in the napari toolbar and draw outlines around regions of interest.
3. Select one or more completed shapes.
4. Type a label in the "Annotation type:" field and click "Assign to selected shapes".
5. Repeat steps 2–4 for each distinct region type.
6. Click a colour cell in the annotation type table to customise per-type colours.
7. Click "Export GeoJSON..." to save annotations for sharing or downstream analysis.

## Notes

- Annotations are persisted to `sdata_cached.zarr` automatically and restored on the next launch.
- Importing a GeoJSON file appends shapes to existing annotations; it does not replace them. Merge conflicts must be resolved manually.
- Annotation types defined here are available to the Annot Nhood and Annot Dist analysis tabs.
