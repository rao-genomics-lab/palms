# Preprocess

The Preprocess tab controls how the expression matrix is normalised before any analysis reads it. It sits after [QC](Tab-QC) in the Tools group because that is the order the two steps run in: QC decides *which cells*, Preprocess decides *on what scale*, and both change what every later tab is about.

Normalisation is the one step with no tab of its own. It runs when something needs it — [Clustering](Tab-Clustering), [Rank Genes](Tab-Rank-Genes), [Markers](Tab-Markers), [UMAP](Tab-UMAP), [Correlation](Tab-Gene-Correlation) and the spatial statistics all ask for a normalised copy — so its settings had nowhere to live and its scaling target was fixed in the template text. This tab is that home.

## Controls

Each control's tooltip names the template parameter its value lands in, so a caption here can be traced to the corresponding argument in the exported notebook.

| Control | Description |
|---|---|
| Median counts (scanpy default) | Checkbox (default unticked). Scales every cell to the *median* count across cells, by passing no `target_sum` at all. This is `scanpy`'s own default and what most published scanpy/squidpy pipelines use. Ticking it disables the spinbox and swaps the recorded line for the one-argument form. |
| Counts per cell | Spinbox (100–1,000,000, default 10000). Every cell is scaled to this many counts before `log1p`. 10,000 is the long-standing convention and the viewer's historical behaviour. Sets `target_sum`. |
| Records line | Read-only. The exact `sc.pp.normalize_total(...)` call the next analysis would record, updated as you change the setting. |

## Workflow

1. Set the target *before* running an analysis, or accept the default of 10,000.
2. Run whatever you were going to run. The normalisation step executes as part of it and is recorded once.
3. If you change the setting afterwards, the next analysis re-runs normalisation and the results that were computed on the old scaling show as stale in the [Notebook](Tab-Notebook) tab.

## Notes

- **Nothing recomputes when you change the setting.** There is no Apply button because there is nothing to apply on its own: the step runs when an analysis asks for `adata_norm`. The status bar says so when you change it while a normalised copy already exists.
- **Changing the target marks earlier results stale, and should.** They were computed on a differently scaled matrix. This is the ordinary `upsert` behaviour — the step's code changed, so what depended on it is no longer current.
- **The two conventions are two blocks, not a parameter that can be `None`.** A `target_sum=None` in the recorded cell is legal scanpy and means the median, but it reads as an argument someone forgot to fill in. Picking the median swaps in a block that passes no argument and says why, so the notebook explains itself.
- **ROI DEG and inferCNV are not yet affected.** Both normalise their own copies of the matrix rather than reading `adata_norm`, and both still scale to 10,000 whatever is set here. The tab says so. Making every path honour one setting is a larger question — how many normalisation policies may one session hold? — and is tracked separately.
- **The setting is saved with the session** and restored on reopening, like the QC cutoffs. A dataset that has never been through this tab opens on 10,000, the viewer's historical behaviour, so nothing recorded before this tab existed is disturbed.
- **Upgrading changes the recorded text once.** The target used to be written as `1e4` and is now substituted as a value, so it records as `10000.0`. The first analysis you run after upgrading re-records the normalisation step with that different text and marks its dependents stale — on an unchanged analysis whose results are numerically identical.
