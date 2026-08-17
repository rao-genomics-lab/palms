# Recovering a Cache

**Prerequisites:** Xenium Viewer installed; a dataset whose zarr cache the viewer is complaining about, or that you simply want to check

**Time required:** ~10 minutes for a check and repair; considerably longer if a full rebuild turns out to be necessary

---

## Overview

Your clusterings, ROIs, registrations and cluster names live in the dataset's zarr cache. If a write is interrupted — a crash, a full disk, a machine going to sleep mid-save — that cache can end up in a state the viewer will not open, and it looks as though the work is gone.

It usually is not. The viewer never deletes a cache: every recovery path is a rename, previous copies of elements are kept aside rather than removed, and repair is attempted automatically before you are asked anything. This tutorial covers what to do when the automatic attempt is not enough, in the order you should try things.

## Steps

### 1. Launch, and read the dialog if one appears

Open the dataset as usual:

```bash
xenium-viewer /path/to/xenium/output/
```

If the cache is fine — or was repairable on its own — the viewer opens normally and you can skip to step 3 to confirm.

If it is not, a dialog appears offering three choices. **Which one you pick matters, so do not dismiss it quickly:**

| Choice | Use it when |
|--------|-------------|
| Rebuild and restore my data | The default, and the right answer almost always. The damaged cache is set aside, a fresh one is built from the raw Xenium files, and your analyses are copied across. |
| Rebuild without restoring | You suspect the saved analyses themselves are what is broken, and you would rather start clean. The old cache is kept as `sdata_cached_prev_<timestamp>.zarr`, so this is not a destructive choice. |
| Quit | You want to back up the dataset directory before anything is touched. Nothing is changed. |

The dialog's details section carries the health report explaining what was actually wrong.

### 2. If you quit, back up first

Copy the cache aside before trying again. It is the only copy of your session state:

```bash
cp -r /path/to/xenium/output/sdata_cached.zarr /path/to/backup/
```

Then relaunch and choose **Rebuild and restore my data**.

### 3. Check the cache from inside the viewer

Once the viewer is open, go to **Tools → Cache** and click **Verify (read-only)**. This changes nothing and is safe to run at any time.

Read the report. The line you want is `✓ Cache is healthy.` — if you see it, you are done. Otherwise the report names the problem, and the next two steps follow from which one it is.

### 4. Re-consolidate, if elements are present but unlisted

If the report says elements are **present on disk but missing from metadata**, the data is intact and only the index is wrong. This is the most common form of damage and the cheapest to fix.

Click **Re-consolidate Metadata**. It finishes any interrupted write, clears leftover files and rebuilds the index. It does not touch your data, and it is reversible.

Click **Verify (read-only)** again to confirm.

### 5. Recover, if elements are missing from disk

If the report says elements are **listed in metadata but missing from disk**, look for `(backup available)` beside them. That means a previous copy was kept.

1. Open the **Backup** dropdown. Entries ending `(previous version)` restore a single element; entries ending `(whole cache)` salvage everything recoverable out of an older store.
2. Choose one and click **Recover from Backup...**, then confirm.
3. When the recovery finishes, **accept the reload it offers.** Recovered data is on disk but not yet in memory, and the viewer does not re-read the store on its own — declining means the recovery looks as though it did nothing until you reopen the dataset.

A whole-cache recovery brings back ROIs and annotations, clustering and CNV columns, session state including the H&E and ARMS transforms, and the analysis sidecars. Anything already present in the live cache is left alone, so recovering twice is harmless.

### 6. Force a rebuild, as a last resort

If none of the above worked, click **Force Rebuild + Restore...**.

This does **not** rebuild there and then. It renames the cache to `sdata_cached_prev_<timestamp>.zarr` and stops — a live rebuild would pull the store out from under every layer in the running viewer. Nothing is deleted.

Quit and relaunch on the same dataset. The rebuild happens at startup, and you are asked whether to restore your data from the cache that was moved aside.

If the disk has less free space than roughly 1.1× the cache size, the rebuild refuses to start rather than filling the disk. Free some space and try again.

### 7. Clean up

A successful rebuild leaves the old cache behind on purpose. Once you have confirmed your analyses are back — check the clustering dropdowns, the ROI layer and the H&E alignment — you can reclaim the space.

Go to **Tools → Dataset**, click **Scan Dataset**, and look under **Backups & trash**. Tick the `sdata_cached_prev_*.zarr` you no longer need and click **Delete Selected...**.

Do this only once you are satisfied, because it removes exactly what a future recovery would draw on.

---

## Notes

- **Copying a dataset does not damage its cache.** Freshness is decided by a content hash of `experiment.xenium`, not by file timestamps, so `cp`, `rsync` and re-downloads no longer make a good cache look stale.
- If you are reporting a problem, click **Copy Report** and **Open Log File** in the Cache tab first. Together they describe both the state of the store and what the session did to reach it.
- The Cache tab's header counts write failures for the current session. A non-zero count there is worth investigating even when the store still verifies clean.
- Analysis outputs that are expensive to recompute — normalised expression, CNV profiles, cached DEG tables, the analysis provenance graph — live in `viewer_cache/` beside the store, not inside it. A cache rebuild does not touch them.
- Running with `--no-cache` disables every action in the Cache tab, since there is no cache to act on.

---

## Next steps

- [Tab-Cache](Tab-Cache) — full reference for every control used here
- [Tab-Dataset](Tab-Dataset) — what else is on disk, and what is safe to remove
- [Installation](Installation) — the startup dialog and other troubleshooting
