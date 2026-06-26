# Correlation

The Correlation tab computes Pearson and Spearman correlation between the expression of two selected genes across cells, with optional normalisation and cluster filtering, and displays a scatter plot alongside the computed statistics.

![Gene Correlation](screenshots/tab-gene-correlation.png)

## Controls

| Control | Description |
|---|---|
| Gene A | First gene to correlate |
| Gene B | Second gene to correlate |
| Normalisation | How to prepare expression values before computing correlation: `Raw counts`, `Fraction of total`, or `Log1p(CPM)` |
| Filter by current cluster selection | When checked, restricts the analysis to cells currently visible in the Cell Coloring filter |
| Plot Correlation | Computes Pearson r and Spearman rho, generates a scatter plot with a statistics box, saves the figure, and displays correlation values in the status bar |

## Workflow

1. Select **Gene A** and **Gene B** from their respective dropdowns.
2. Choose a **Normalisation** method appropriate for your comparison.
3. Optionally enable **Filter by current cluster selection** to restrict the analysis to a cell subset defined in the Cell Coloring tab.
4. Click **Plot Correlation** to generate the figure and read off the statistics.

## Notes

- At least two genes must be present in the dataset for this tab to be functional.
- If cluster filtering is enabled but no filter is active in Cell Coloring, all cells are used.
- The saved figure includes both correlation coefficients and the number of cells used.
