# Rank Genes

The Rank Genes tab lets you compute per-cluster marker genes using differential expression testing, visualise results as dotplots, rank-gene plots, and volcano plots, and annotate clusters automatically using CellTypist, a language model, or label transfer from a reference dataset.

<!-- SCREENSHOT: docs/screenshots/tab-rank-genes.png -->

## Controls

### Rank Genes Computation

| Control | Description |
|---|---|
| Clustering | Clustering whose groups define foreground vs. background in the statistical test |
| Method | Statistical test: `wilcoxon`, `t-test`, or `logreg` |
| Top N genes | Slider (5–50, default 25) — number of top marker genes to retain per cluster |
| Run Rank Genes | Executes marker gene analysis; results are required before plots or annotation are available |

### Visualisation

| Control | Description |
|---|---|
| Genes per cluster | Slider (3–20, default 5) — number of genes shown per cluster in the dotplot |
| Dendrogram | When checked (default), includes a cluster dendrogram in the dotplot |
| Show Dotplot | Generates a scanpy dotplot (mean expression + fraction detected) and saves it |
| Show Rank Genes Plot | Generates a multi-panel rank-genes overview plot and saves it |
| Edit Cluster Labels... | Opens a dialog to rename clusters; labels propagate to all exports and plots |
| Reset Labels | Clears all custom cluster label edits and restores the original cluster names |
| Export Full Results CSV... | Saves the full ranked gene table to a CSV file |
| Generate All Volcano Plots... | Saves one volcano plot per cluster pair to a chosen directory |

### CellTypist Annotation

Requires the `celltypist` package.

| Control | Description |
|---|---|
| CellTypist Model | Select from downloaded models |
| Min confidence | Float slider (0.0–1.0, default 0.5) — minimum per-cell prediction confidence to accept |
| Download Models | Downloads available CellTypist models from the internet |
| Annotate with CellTypist | Runs annotation; predictions are majority-voted per cluster and stored as cluster labels |

### LLM Annotation

| Control | Description |
|---|---|
| LLM Provider | Language model to use: `Claude (claude)`, `Gemini (gemini)`, or `Codex (codex)` |
| Annotate with LLM | Sends top-N genes per cluster to the LLM and stores returned cell-type names as cluster labels; requires the relevant API key in the environment |

### Label Transfer

| Control | Description |
|---|---|
| Reference Dataset / Browse... | Built-in reference dropdown or a custom `.h5ad` file via file dialog |
| Annotation Column | Column in the reference `obs` table to transfer |
| Load Reference | Loads the selected reference h5ad into memory |
| Run Label Transfer | Runs `sc.tl.ingest()` against the reference; predictions are majority-voted per cluster |

## Workflow

1. Select a clustering from the **Clustering** dropdown.
2. Choose a **Method** and set **Top N genes**, then click **Run Rank Genes**.
3. Adjust **Genes per cluster** and toggle **Dendrogram**, then use **Show Dotplot** or **Show Rank Genes Plot** to visualise results.
4. Optionally rename clusters with **Edit Cluster Labels...** before exporting.
5. To annotate clusters automatically, choose one of the three annotation methods:
   - **CellTypist** — download a model, set confidence threshold, click **Annotate with CellTypist**.
   - **LLM** — select a provider and click **Annotate with LLM** (API key required).
   - **Label Transfer** — load a reference dataset, pick an annotation column, click **Run Label Transfer**.
6. Export ranked genes with **Export Full Results CSV...** or save volcano plots with **Generate All Volcano Plots...**.

## Notes

- You must click **Run Rank Genes** before any visualisation or annotation method becomes available.
- CellTypist and LLM annotation both use the ranked gene lists computed in step 2.
- Label Transfer does not require rank genes; it operates directly on the expression matrix.
- Cluster labels set here propagate to exports in all other tabs.
