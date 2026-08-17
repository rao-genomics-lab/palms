# Templates

Inspect and edit the analysis code every other tab runs. Each analysis step is a template with a declared contract, made of named blocks; this tab shows the shipped source, lets you override it, and previews the exact string that would be executed with the parameters currently set in the owning tab. The preview is not a reconstruction of what might run — it calls the same render method the executor calls, so the code shown is the code that would run. This is the **Templates** tab in the "Tools" control panel group.

<!-- SCREENSHOT: docs/screenshots/tab-templates.png -->

## Controls

### Template list

| Column | Description |
|---|---|
| Template | Every template, grouped by what it belongs to: **Setup** (`normalize`, `spatial_neighbors`), **Clustering**, **Genes**, **Spatial** and **ROI**. |
| Blocks | The number of blocks in the template, or a status badge: **● customised (n)** for one you have edited, **⚠ review** for one whose shipped version has changed since you forked it, and **✕ not used** for a customisation that failed validation and is therefore being ignored. |

### Panes

| Pane | Description |
|---|---|
| Contract | Read-only. The template's documentation line, then its required and optional `params`, the namespace names it `requires`, the results it `outputs`, its block names, how many assemblies it declares, and which blocks are frozen. |
| Default (read-only) | The shipped template body, exactly as installed. |
| Yours (editable) | Your version, and the only editable pane. Blocks are separated by `#--- block <name>` marker lines. Empty until you type something — an unmodified template shows the shipped body here. |
| Preview — what would run | The rendered result: blocks assembled and parameters substituted. When the owning tab has registered a preview provider the header reads `# preview — current widget values` and the numbers are the ones in that tab right now; otherwise it reads `# preview — sample values` and uses synthesised literals of the right type. |
| Diff — yours vs the new default | Replaces the preview pane while a template is flagged **⚠ review**, showing your version against the shipped one that has moved on. |
| Problems | Errors and warnings from the last validation. Empty when the template is sound. |

### Buttons

| Control | Description |
|---|---|
| Validate | Parses your edit, merges it onto the shipped blocks and checks it against the contract. Runs in milliseconds and needs no data, so it can be used freely while editing. |
| Save & Activate | Writes the override to your config directory. **It saves even when validation fails** — what is gated is activation, not writing, so an invalid edit is preserved for you to fix rather than forcing you into an external editor. |
| Revert to default | Removes the override entirely and goes back to the shipped template. |
| Take new default for changed blocks | Only visible while a template is flagged **⚠ review**. Replaces just the blocks whose shipped source has changed, keeping your other customisations. |

## Workflow

- **To see what a button will actually run**, select its template and read the **Preview — what would run** pane. Set the parameters in the owning tab first and they will be reflected here.
- **To change what it runs**, edit the **Yours (editable)** pane, click **Validate**, and then **Save & Activate**. The status line tells you whether the edit is active or whether the shipped template is still running in its place.
- **To edit one part of a template only**, delete the blocks you do not care about from your version. Resolution is per block, so an omitted block keeps tracking the shipped template and continues to receive upstream fixes.
- **To undo**, click **Revert to default**. Saving a body that is identical to the shipped one has the same effect.
- **After an upgrade**, look for **⚠ review** badges. Your edit is still active; the badge means the shipped block it was forked from has since changed, and the update might be a fix your version is now shadowing. Read the diff, then either click **Take new default for changed blocks** or save again to confirm you have reviewed it.
- **When a result looks wrong**, relaunch with `xenium-viewer <path> --no-user-templates` to run the shipped templates for one session. That is the fastest way to tell a customisation apart from a genuine problem.

## Notes

- Overrides are written to the user config directory — on Linux typically `~/.config/xenium-viewer/templates/<template.id>.tmpl`. The directory is not created until your first save. A manifest alongside them records which shipped version each edited block was forked from, which is what makes the **⚠ review** flag possible.
- **An override is resolved per block, not per file.** This is what lets you change one step of an analysis without freezing the rest of it at the version you forked.
- **An invalid override is skipped, never fatal.** The shipped template runs, the tab badges it **✕ not used**, a napari warning appears once per session, and the [Cache](Tab-Cache) tab's header counts it among the session's failures.
- The strictest check is that a **required parameter the template no longer mentions is a hard error** — a template like that runs, succeeds, and silently ignores the setting you chose in the GUI, which is worse than failing.
- Blocks are delimited by `#--- block <name>` lines. Nothing may be added to a marker line: everything after the block name is part of the name, so a note there creates a differently-named block. Frozen blocks are listed in the **Contract** pane for the same reason.
- **The diff for a stale template is two-way, not a three-way merge.** Conflict markers in Python source would be a syntax error rather than a visible annotation, so none are inserted.
- A missing or unreadable manifest means nothing is flagged for review. That is deliberate: warning about every override after every upgrade, with nothing specific to point at, teaches people to dismiss the warning.
- `XENIUM_VIEWER_TEMPLATE_PATH` (colon-separated) replaces the user scope for both reading and writing; setting it to an empty value disables overrides entirely. Saving derives its destination from the same search path that reading uses, so a write cannot land somewhere the reader does not look.
- Twelve of the fourteen templates supply live parameters to the preview. `normalize` takes no parameters, and `spatial_neighbors` takes its neighbour count from whichever tab requested the graph, so both fall back to declared sample values.
- What a step actually ran, as opposed to what it would run, is recorded in the provenance graph and readable in the [Notebook](Tab-Notebook) tab. An exported notebook containing any non-shipped template opens with a banner saying so.
