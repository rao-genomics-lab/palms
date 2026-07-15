# CNV

The CNV tab infers copy-number variation from expression data using the [InSituCNV](https://github.com/Moldia/InSituCNV) method: pick an existing clustering, mark some of its categories as the "normal" reference population, and run inference to get a colourable CNV-subclone clustering, a continuous per-cell CNV score, and a chromosome heatmap.

Requires the `cnv` optional extra — see [Installation](Installation).

## Controls

### Reference Population

| Control | Description |
|---|---|
| Reference clustering | Existing clustering/annotation to source the "normal" reference population from |
| Cluster checkboxes | Three-column grid — check which categories of the selected clustering count as reference ("normal") cells; unchecked by default |
| Select All | Checks all reference checkboxes |
| Deselect All | Unchecks all reference checkboxes |

### Parameters

| Control | Description |
|---|---|
| Neighbors (expression graph) | Spin box (5–100, default 15) — neighbors for the expression PCA graph used for smoothing |
| Smoothing neighbors | Spin box (5–200, default 30) — neighbors used by InSituCNV's graph-smoothing step |
| Window size (genes) | Spin box (2–200, default 10) — infercnvpy sliding-window size |
| Window step | Spin box (1–50, default 2) — infercnvpy sliding-window step |
| CNV cluster resolution | Float spin box (0.05–2.0, default 0.2) — Leiden clustering resolution for CNV subclones |

### Run and Results

| Control | Description |
|---|---|
| Run CNV Inference | Runs the pipeline: normalize/PCA/neighbors → smoothing → genomic-position mapping → inferCNV → CNV-profile clustering |
| Results area | Reports reference clustering/categories, genes mapped to the genome, windows produced, CNV clusters found, and CNV score range |
| Save Chromosome Heatmap (PDF/PNG) | Saves an infercnvpy chromosome heatmap for the CNV clusters to `plots/cnv_heatmap.png` and `plots/cnv_heatmap.pdf` in the dataset directory; enabled after a run |
| Color Cells by CNV Score | Recolours the labels layer by the continuous per-cell CNV score (viridis); enabled after a run |

## Workflow

1. Select a **Reference clustering** — an existing clustering or cell-type annotation with a clear normal/non-tumor population (e.g. an immune or stromal cluster from Rank Genes annotation, or a built-in `graphclust`).
2. Check the categories that represent the **reference (normal)** population in the cluster checkbox grid.
3. Adjust **Window size**/**Window step** if needed — the defaults are tuned lower than infercnvpy's own bulk-RNA-seq defaults (60/10) because Xenium gene panels are much smaller; the results panel reports how many genes mapped to the genome and how many windows were produced so you can judge whether to widen or narrow further.
4. Click **Run CNV Inference**. On the first run, infercnvpy's human (GRCh38) gene-position reference is downloaded and cached automatically (requires internet access).
5. Once complete, the result is registered as a new clustering (`cnv_leiden_<resolution>`) and automatically applied as the active coloring — it's also usable everywhere clusterings are, e.g. as the groupby in Rank Genes or the cluster source in ROI DEG.
6. Click **Save Chromosome Heatmap (PDF/PNG)** to export gain/loss patterns across chromosomes per CNV cluster as `plots/cnv_heatmap.png` and `plots/cnv_heatmap.pdf` (the heatmap is not displayed in a window — building it can be slow, so it's saved directly instead).
7. Click **Color Cells by CNV Score** to switch the spatial view to a continuous CNV-burden overlay instead of discrete clusters.

## Notes

- Human genome only for now — infercnvpy's default GRCh38 gene-position reference is used; non-human panels will report 0 genes mapped.
- The reference population must be explicitly chosen; no cluster is treated as "normal" by default, and at least one reference category is required to run.
- The full CNV profile (per-window values, gene positions) is cached alongside the zarr store so the chromosome heatmap can be regenerated after reopening the dataset without recomputing the pipeline.
- Results persist across sessions like any other clustering; reopening the dataset restores the results summary and re-enables the heatmap/score-coloring buttons.
