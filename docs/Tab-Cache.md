# Cache

Check the health of the dataset's zarr cache, rebuild its metadata index, recover user-generated data out of a backup, and force a full rebuild from the raw Xenium files. Every action here is one the loader already attempts automatically at startup; this tab makes them available deliberately, and shows you what it found. This is the **Cache** tab in the "Tools" control panel group.

<!-- SCREENSHOT: docs/screenshots/tab-cache.png -->

## Controls

| Control | Description |
|---|---|
| Header block | The cache path, its size and the free space on that disk, when it was built and with which spatialdata version, an overall **Status: healthy** / **needs attention** line, a count of any write failures this session, and the path to the session log. |
| Verify (read-only) | Reads the store and reports on it. Changes nothing, and is safe to run at any time — including on a store too broken for spatialdata to open, because it parses the root metadata as plain JSON rather than opening the store. |
| Report area | The result of the last action. Lists interrupted writes still to finish, elements listed in the metadata but missing from disk (annotated when a backup exists), elements present on disk but missing from the metadata, invalid entries safe to drop, leftover files, and then an inventory of elements, sidecars and available backups. |
| Re-consolidate Metadata | Finishes or unwinds any interrupted write, clears debris, and rebuilds the root metadata index. Reversible, and does not touch your data. This fixes the most common form of corruption: an element that is on disk but absent from the index, so the viewer cannot see it. |
| Backup | Dropdown of recovery sources, filled by **Verify**. Two kinds of entry: `<element>  (previous version)` restores one element from `.xv_trash`; `<store name>  (whole cache)` salvages from a sibling backup store. |
| Recover from Backup... | Disabled when the dropdown is empty. Copies data from the selected source into the live cache, leaving anything already present alone. Confirms first. |
| Force Rebuild + Restore... | Moves the current cache aside so the next launch rebuilds it from the raw Xenium files. Confirms first, and warns if the disk has too little free space. **Requires restarting the viewer** — see the notes. |
| Open Log File | Opens this session's rotating log in the system's default viewer. |
| Copy Report | Copies the report area to the clipboard, for pasting into a bug report. |

## Workflow

1. Click **Verify (read-only)** and read the report. If it ends `✓ Cache is healthy.`, there is nothing to do here.
2. If elements are reported present on disk but missing from the metadata, click **Re-consolidate Metadata**. This is the cheap fix and resolves most problems.
3. If elements are reported missing from disk, choose the matching entry in **Backup** and click **Recover from Backup...**. Entries marked `(previous version)` restore a single element and are reversible; `(whole cache)` entries salvage everything recoverable out of an older store.
4. Accept the reload the tab offers afterwards. Recovered data is on disk but not yet in memory, and the viewer does not re-read the store on its own.
5. Use **Force Rebuild + Restore...** only when the first three steps have not worked. Then quit and relaunch the viewer on this dataset — the rebuild happens at startup, and it will ask whether to restore your data.
6. When reporting a problem, click **Copy Report** and **Open Log File** — together they describe both the store's state and what the session did to reach it.

## Notes

- **Force Rebuild does not rebuild in the session.** It renames the store to `sdata_cached_prev_<timestamp>.zarr` and stops, because every napari layer and manager in a running viewer points into the live store; rebuilding underneath them would leave them pointing at freed memory. The rebuild happens on the next launch. Nothing is deleted at any point.
- The rebuild refuses to start when free space is under roughly 1.1× the store size — it keeps the old copy alongside the new one, so it would rather stop than fill the disk.
- **Recovery never opens the backup as a SpatialData object**, deliberately: a cache worth recovering from is usually one that will not open. It works at the filesystem and zarr-array level instead.
- A whole-cache recovery pulls across, in order: user shapes and images (ROIs, annotations, landmarks, tiles); clustering and CNV `obs` columns, reindexed onto the live cells; the session attributes and the H&E/ARMS affines; and the sidecar files from `viewer_cache/`. Anything already present in the live cache is left untouched, so recovering twice is harmless.
- Because the H&E and ARMS transforms are rewritten from memory when the session is saved, a recovery also loads them back into the running viewer — otherwise the images would return unaligned at the next save.
- **Freshness is a content hash, not a timestamp.** The cache carries a sha256 of `experiment.xenium` in its manifest, so copying or re-downloading a dataset no longer condemns a good cache. A cache built before manifests existed falls back to comparing modification times, which is treated as uncertain: it asks rather than rebuilding.
- Every action here takes the store lock, so nothing in this tab can race a background write from another tab.
- All four actions are disabled under `--no-cache`, and when no cache exists for the dataset.
- To reclaim the space left behind by a rebuild — the moved-aside `sdata_cached_prev_*.zarr` — use the [Dataset](Tab-Dataset) tab. Be aware that doing so removes what a future recovery would draw on.
- See [Installation](Installation) for the startup dialog that appears when a cache is found to be unreadable before the viewer opens, and [Recovering a Cache](Tutorial-Recovering-a-Cache) for the whole procedure end to end.
