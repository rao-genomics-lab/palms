# Nhood Enrich

The Nhood Enrich tab tests whether cell clusters tend to be spatially co-located or segregated using permutation-based neighbourhood enrichment (squidpy), producing a Z-score matrix that summarises spatial associations between all cluster pairs.

<!-- SCREENSHOT: docs/screenshots/tab-neighborhood-enrichment.png -->

## Controls

| Control | Description |
|---|---|
| Clustering | Clustering to test for spatial neighbourhood enrichment |
| Permutations | Slider (100–1000, default 1000) — number of permutations used to derive Z-scores |
| N neighbors | Slider (3–20, default 6) — neighbourhood size used to define spatial proximity |
| Run Nhood Enrichment | Runs the permutation test and computes the Z-score and count matrices |
| Display mode | Controls which matrix is shown in the heatmap: `zscore` or `count` |
| Show Heatmap | Displays the enrichment heatmap; enabled after running the analysis |
| Export Z-scores CSV... | Saves the Z-score matrix to a CSV file; enabled after running |
| Results area | Shows the top 10 most enriched and most depleted cluster pairs |

## Workflow

1. Select a **Clustering** from the dropdown.
2. Adjust **Permutations** and **N neighbors** as needed.
3. Click **Run Nhood Enrichment** and wait for completion.
4. Choose a **Display mode** (`zscore` or `count`), then click **Show Heatmap**.
5. Review the top enriched/depleted pairs in the results area.
6. Click **Export Z-scores CSV...** to save the full matrix.

## Notes

- Positive Z-scores indicate that two clusters are co-localised more than expected by chance; negative Z-scores indicate spatial segregation.
- Results are persisted to `adata.uns['nhood_enrichment']` and restored automatically on the next load.
