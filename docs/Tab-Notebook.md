# Notebook

An interactive Python code editor embedded in the viewer, with per-cell execution, inline output display (print output, return values, and figures), and synchronisation with the provenance graph that other tabs record their actions into. This tab is in the "Tools" control panel group.

![Notebook](screenshots/tab-notebook.png)

## Controls

### Toolbar

| Control | Description |
|---|---|
| Sync Graph | Rebuilds the recorded cells from the provenance graph, in dependency (topological) order. Steps recorded by other tabs since the last sync appear as cells, and any step whose inputs have since changed is flagged with a ⚠ stale badge. |
| + Cell | Adds a new blank code cell at the bottom of the notebook. |
| Run All | Executes all cells from top to bottom in order. |
| Clear Outputs | Hides all cell output areas. Cell code is not affected. |
| Show DAG | Opens a rendered diagram of the provenance graph, showing each step and the dependencies between them. |
| Export .ipynb | Writes the graph out as `analysis_notebook.ipynb` in the dataset directory. The notebook is code-only and replays from the raw Xenium output. |

### Per-cell controls

| Control | Description |
|---|---|
| Code editor | Monospace, auto-resizing text area. Write or edit Python code for this cell. |
| Run | Executes this cell and displays output inline below it. |
| Delete | Removes this cell from the notebook. |

### Output area

Shown below each cell after execution:

| Output type | Display |
|---|---|
| Return value | Shown as `Out: ...` |
| Standard output | Captured and displayed in monospace. |
| Matplotlib figures | Rendered as inline images. |
| Standard error | Displayed in orange. |
| Exceptions | Displayed in red with a traceback. |

## Workflow

- A welcome cell listing available objects is inserted automatically when the tab is first opened.
- Actions taken in other tabs (clustering, DEG, spatial analysis, and so on) are recorded into the provenance graph and appear in the notebook when you click "Sync Graph".
- Cells share a single Python execution context: variables defined in one cell are available in all subsequent cells.
- Write exploratory code directly in new cells added with "+ Cell", then run them individually or with "Run All".
- To capture a Matplotlib figure inline, assign it to a variable or call `plt.show()` at the end of the cell.

## Notes

- Available objects in the notebook context: `adata`, `sdata`, `viewer`, `ctx`, `clusterings`, `color_manager`, `gene_names`, `data_path`.
- The provenance graph **is** saved with the session, in `sdata_cached.zarr/viewer_session/`, and restored when you reopen the dataset — so recorded steps accumulate across sessions into a single notebook. Free-form cells you type into the tab yourself are scratch space and are not persisted; move anything worth keeping into a file, or rely on the recorded steps.
- Two files are written into the dataset directory automatically: `analysis.py` (a flat script) and `analysis_notebook.ipynb` (the notebook export). Both replay from the raw Xenium output.
