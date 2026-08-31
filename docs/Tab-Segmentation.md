# Segmentation

Replace the native Xenium cell segmentation with a custom segmentation produced by a separate preprocessing pipeline, enabling re-analysis with alternative cell boundaries. This tab is in the "Tools" control panel group.

![Segmentation](screenshots/tab-segmentation.png)

## Controls

| Control | Description |
|---|---|
| Active segmentation (read-only) | Displays the current segmentation — "Xenium (native)" or "Custom (from cache/file)" — along with the corresponding cell and gene counts. |
| Load Custom Segmentation... | Opens a file dialog for a `custom_segmentation.h5ad` file. If a cached version already exists in `sdata_cached.zarr`, prompts you to load from cache (fast) instead of re-importing the original file. |
| Revert to Xenium Segmentation | Restores the original Xenium cell labels and AnnData object. Disabled until a custom segmentation is loaded. |
| Update SpatialData on disk | Saves the current segmentation state (custom or native) back to `sdata_cached.zarr`. |

## Workflow

1. Produce a `custom_segmentation.h5ad` file using the `palms-build-custom-segmentation` pipeline. This step requires R and Seurat for cell boundary extraction; see the project README for details.
2. In the viewer, open the Segmentation tab and click "Load Custom Segmentation...".
3. Select the `custom_segmentation.h5ad` file. If a cached copy exists, choose whether to load from cache or re-import.
4. The viewer swaps to the custom segmentation; previously computed cluster-dependent analyses are cleared.
5. Re-run Leiden clustering and any downstream analyses on the new segmentation.
6. To revert, click "Revert to Xenium Segmentation".

## Notes

- A custom segmentation is cached in `sdata_cached.zarr` on first load; subsequent loads use the cached version automatically.
- Building a custom segmentation is a two-stage pipeline. The boundaries are extracted first (typically in R/Seurat), then `palms-build-custom-segmentation` turns them into the label raster the viewer loads. Its output is a `custom_labels.zarr` that must stay **beside** the `.h5ad` you select here — selecting an `.h5ad` on its own fails with "not found alongside h5ad".
- Swapping segmentation resets active clusterings. You will need to re-run Leiden clustering on the new segmentation before cluster-dependent visualisations and analyses are available.
- The swap is recorded. The exported notebook's load cell binds `adata` from the custom
  table cached in `sdata_cached.zarr` rather than from the Xenium one, so a replay is
  about the same cells you were looking at. Results recorded *before* the swap are marked
  stale in the Notebook tab for the same reason the viewer clears them here: they were
  computed on cells that are no longer loaded.
