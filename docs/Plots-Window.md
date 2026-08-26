# Plots Window

Every figure the viewer produces — dotplots, UMAPs, neighbourhood-enrichment heatmaps, co-occurrence curves, ligand–receptor dotplots, CNV heatmaps, the provenance DAG — appears in one dockable **Plots** panel and is written to `<dataset>/plots/`.

## Opening it

The dock is hidden until the first figure arrives, then reveals itself at the bottom of the window. **View → Show Plots** (`Ctrl+Shift+P`) toggles it, and closing it with its own close button unticks the menu item.

## The gallery

Figures are listed newest first. Each card shows a thumbnail, the figure's title, and the files it was written to.

| Button | What it does |
|---|---|
| Open | Opens the figure in its own resizable window with matplotlib's pan / zoom / save toolbar. |
| Save as… | Writes an extra copy wherever you choose, in whatever format the extension names. |
| Remove | Drops the figure from the gallery. Files already written are untouched. |
| Clear all | Empties the gallery. |

The gallery keeps the twenty most recent figures; beyond that the oldest is dropped. Matplotlib figures are not small, and a long session produces a lot of them.

## Where the files go

Every plot is written to `<dataset>/plots/`, in each format selected under **Preferences → Plot format**:

| Setting | Files written per plot |
|---|---|
| **PNG + PDF** (default) | a 300 dpi raster to look at, and a vector version to publish |
| PNG | one raster |
| PDF | one vector file |
| SVG | one vector file, for further editing |

Names are keyed by what the figure is about — `dotplot_leiden_r1.0.png`, `nhood_enrichment_graphclust.pdf`, `umap_EPCAM_KRT5.png` — so a second run against a different clustering does not overwrite the first. A run against the *same* clustering does overwrite it, which is usually what you want.

Two kinds of output do not go through the gallery:

- **Pairwise volcano batches** ([Rank Genes](Tab-Rank-Genes), [ROI Analysis](Tab-ROI-Analysis), [ARMS Overlay](Tab-ARMS-Overlay)) produce one figure per cluster pair — dozens at once, which would bury everything else. They still ask where to put them, defaulting to `<dataset>/plots/volcano_<clustering>/`, and honour the format setting.
- **CopyKAT heatmaps** are drawn by a detached background process in a separate conda environment, which cannot reach the GUI. They are written straight to `<dataset>/plots/`.

## In the notebook

Each figure is recorded as a `plot:*` node in the [provenance graph](Tab-Notebook), and the recorded `savefig` call names the file that was actually written, relative to the dataset directory. A replayed notebook therefore reproduces the figure *and* puts it where the session put it.

Two figures are recorded as **notes** rather than code — the [annotation neighbourhood](Tab-Annot-Nhood) heatmap and the [annotation distance](Tab-Annot-Distance) plot. Both are computed over shapes drawn in the viewer, which a notebook has no access to; the note records that the figure exists and where it went rather than pretending there is code to replay.

## Deleting them

`<dataset>/plots/` is viewer-created, so [Tools → Dataset](Tab-Dataset) lists it as deletable. Nothing else depends on it.
