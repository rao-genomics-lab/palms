# Crop Dataset

Draw one or more polygons in the "Crop Regions" layer and export each as its own standalone, independently-openable xenium-viewer data directory — a cropped morphology image, cell/nucleus labels, transcripts, and AnnData table, ready to `xenium-viewer <output_dir>`. This tab is in the "Tools" control panel group.

![Crop Dataset](screenshots/tab-crop-dataset.png)

## Controls

| Control | Description |
|---|---|
| Activate Draw Polygon Tool | Selects the "Crop Regions" Shapes layer and switches it into polygon-drawing mode, so the next clicks on the canvas start a new polygon. |
| Clear All Regions | Removes all drawn polygons from the "Crop Regions" layer. |
| Crop & Export | Exports every drawn polygon as a separate dataset. Prompts for an output folder and a dataset name for each region, in order, then runs the crop in the background with a progress dialog. |

## Workflow

1. Select the "Crop Regions" Shapes layer in the napari layer list (or click "Activate Draw Polygon Tool").
2. Draw one or more polygon outlines around the tissue regions you want to export. Press Enter to close each polygon.
3. Click "Crop & Export". For each drawn polygon, in order:
   - A folder picker asks where to save the exported dataset.
   - A name prompt asks for the dataset's folder name (e.g. `core_A`). If that folder already exists, you're asked to confirm overwriting it.
4. A progress dialog tracks each region's export (image/label cropping, cell and transcript filtering, writing the dataset, building the transcript cache).
5. A summary dialog reports success or failure per region, with the output path for each successful export.
6. Open a completed export directly: `xenium-viewer <output_dir>/<name>`.

## Notes

- Images and cell/nucleus label rasters are cropped to each polygon's **pixel bounding box** (a rectangle), not the exact polygon outline. Cells and transcripts, however, are filtered to the **exact drawn polygon** via a true point-in-polygon test, so the exported cell table and transcript points are precise even though the image extends slightly beyond the drawn shape.
- `cell_labels` and `nucleus_labels` are both masked so they only contain cells/nuclei kept in the exported table. `nucleus_labels` pixel values are an independent numbering unrelated to `cell_id`/`cell_labels`, so this isn't done by matching ID numbers — a nucleus is kept if it spatially overlaps a kept cell's footprint in the cropped raster. A nucleus straddling the boundary between a kept and an excluded cell is kept if it overlaps the kept cell at all.
- Only the core Xenium elements are exported: morphology image, cell labels, nucleus labels, transcripts, and the cell-by-gene table. H&E/ARMS overlays, custom segmentation, and tissue annotations are not carried over.
- Requires native Xenium segmentation — if a custom segmentation is active (see the Segmentation tab), revert to Xenium segmentation before cropping.
- Each exported dataset is fully standalone: it includes its own `experiment.xenium`, a prebuilt SpatialData zarr cache, `transcripts.parquet`, and a prebuilt transcript feather cache, so it opens at full speed without re-processing. It does not inherit the source dataset's saved session (ROIs, registration) — clusterings are the exception (see below).
- Exported datasets have no raw Xenium output files (no `cells.zarr.zip`, `cell_feature_matrix`, etc.) — only the zarr cache and derived transcript files. They always load from that cache, even when the viewer is launched with `--no-cache`, since there's nothing else to load them from.
- Exported datasets have no `analysis/` folder, so there's no separate precomputed UMAP/clustering CSVs to load. All of the source dataset's clusterings (built-in ones like `graphclust`/`kmeans_*` and any custom/Leiden ones), subset to the exported cells, are carried over automatically and available immediately in the Clustering tab and every analysis tab's clustering selector — no need to recompute after cropping. Cluster labels (the human-readable names assigned via "Edit Cluster Labels...") carry over too, per clustering. If the source dataset's UMAP was already computed at crop time, those coordinates (for the kept cells) also carry over and the UMAP tab still works; if not, the UMAP tab starts empty until you compute one.
- A region with no cells inside it, or an invalid/self-intersecting polygon, is reported as a failure for that region without affecting the rest of the batch.
