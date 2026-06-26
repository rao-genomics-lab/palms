# Lig-Rec

The Lig-Rec tab tests for significant ligand-receptor interactions between spatially adjacent cell clusters using squidpy and the OmniPath database, producing mean interaction strengths and permutation-based p-values.

<!-- SCREENSHOT: docs/screenshots/tab-ligand-receptor.png -->

## Controls

| Control | Description |
|---|---|
| Clustering | Clustering that defines sender and receiver cell types |
| Permutations | Slider (100–1000, default 1000) — number of permutations for significance testing |
| N neighbors | Slider (3–20, default 6) — neighbourhood size used to define spatial proximity |
| OmniPath | Include OmniPath interactions (default checked) |
| LigRecExtra | Include LigRecExtra interactions (default checked) |
| PathwayExtra | Include PathwayExtra interactions (default checked) |
| KinaseExtra | Include KinaseExtra interactions (default checked) |
| CellPhoneDB only | Restrict the analysis to CellPhoneDB interactions only (default unchecked) |
| Run L-R Analysis | Computes mean interaction strengths and permutation p-values across all selected databases |
| P-value threshold | Filters which interactions are shown in the plot: `0.001`, `0.005`, `0.01`, or `0.05` |
| Show L-R Plot | Displays a dotplot of significant interactions; enabled after running the analysis |
| Export Means CSV... | Saves the mean interaction strength matrix to a CSV file; enabled after running |
| Export P-values CSV... | Saves the full p-value matrix to a CSV file; enabled after running |

## Workflow

1. Select a **Clustering** from the dropdown.
2. Adjust **Permutations** and **N neighbors** as needed.
3. Check or uncheck the interaction database checkboxes to define which resources to query.
4. Click **Run L-R Analysis** and wait for completion.
5. Set a **P-value threshold** and click **Show L-R Plot** to inspect significant interactions.
6. Export results with **Export Means CSV...** or **Export P-values CSV...**.

## Notes

- A clustering must be selected before running the analysis.
- The p-value threshold only affects what is shown in the plot; the full result set is always available for export regardless of the threshold.
- Results are persisted to `adata.uns['ligrec']` and restored automatically on the next load.
