# `sq.gr.spatial_neighbors` → `spatial_neighbors_knn` — done

Status: **done (2026-08-17)**, issue #19. Kept as the record of what was measured and
what turned out to be wrong in the original assessment, because the equivalence claim is
the reason no saved result had to be recomputed.

## What changed

Two call sites, changed together — the template is what the notebook replays, the
function is what the GUI runs:

```python
# was
sq.gr.spatial_neighbors(adata_norm, coord_type='generic', n_neighs=k)
# now
sq.gr.spatial_neighbors_knn(adata_norm, n_neighs=k)
```

1. `src/palms/utils/step_templates/builtin/spatial_neighbors.tmpl` — the recorded
   code, reached through `ctx.ensure_spatial_neighbors()` in `tabs/_helpers.py`.
2. `src/palms/utils/spatial_analysis.py::compute_spatial_neighbors` — used
   directly by `tabs/tab_annot_nhood.py`, whose call is unrecorded. Its signature
   (`adata_norm`, `n_neighs=6`) is unchanged, so that tab needed no edit.

Two more that the original note missed: `scripts/generate_docs.py` hard-codes the
`coord_type='generic'` prose, and `docs/Analysis-Templates.md` carries the template body
verbatim — `tests/test_generated_docs.py` byte-compares the regenerated page, so the
generator must be edited and `python scripts/generate_docs.py` re-run.

`adata_norm` is passed **positionally**: the replacement names its first parameter
`data`, not `adata`, and everything after it is keyword-only.

## Why it was safe

squidpy 1.8.2 warned that `spatial_neighbors` would be removed in 1.9.0, naming five
replacements. The graph is a *dependency*, not a leaf feature, so the concern worth
checking was whether the swap changes results — which would invalidate saved
neighbourhood-enrichment and ligand-receptor output.

It does not. Measured on 500 random points under squidpy 1.8.2:

| | `spatial_neighbors(coord_type='generic', n_neighs=6)` | `spatial_neighbors_knn(n_neighs=6)` |
|---|---|---|
| `obsp` keys | `spatial_connectivities`, `spatial_distances` | same |
| `uns` keys | `spatial_neighbors` | same, **including** `params: {'coord_type': 'generic', …}` |
| `spatial_connectivities` | 3000 nnz | 3000 nnz, **identical matrix** |
| `spatial_distances` | 3000 nnz | 3000 nnz, **identical matrix** |

`coord_type='generic'` with `n_neighs=k` *is* k-nearest-neighbours. Re-run the check
below against whatever squidpy version is current before trusting this table again:

```python
import numpy as np, anndata as ad, squidpy as sq
rng = np.random.default_rng(0)
coords = rng.random((500, 2)) * 1000

def build(fn, **kw):
    a = ad.AnnData(rng.random((500, 20)).astype("float32"))
    a.obsm["spatial"] = coords.copy()
    fn(a, **kw)
    return a

old = build(sq.gr.spatial_neighbors, coord_type="generic", n_neighs=6)
new = build(sq.gr.spatial_neighbors_knn, n_neighs=6)
for k in ("spatial_connectivities", "spatial_distances"):
    assert (old.obsp[k] != new.obsp[k]).nnz == 0, k
```

## Four things the original assessment got wrong

Recorded because each was stated confidently and each was checkable:

1. **Co-occurrence was never affected.** The note and `CHANGELOG.md` both said three tabs
   break together. `tabs/tab_co_occurrence.py` declares only
   `deps=[f"clustering:{key}"]` and `spatial.cooccur.tmpl` calls `sq.gr.co_occurrence`,
   which computes its own radii and needs no `obsp` graph. Two tabs, not three:
   Neighbourhood Enrichment and Ligand-Receptor.
2. **The version floor did need raising — and the doc's own evidence was a trap.**
   `spatial_neighbors`'s docstring says `.. deprecated:: 1.7.0`, which reads as "the
   replacements exist from 1.7.0". They do not. Checked against the actual wheels,
   `spatial_neighbors_knn` is **absent from 1.7.0, 1.8.0 and 1.8.1, and first appears in
   1.8.2**:

   ```
   squidpy 1.7.0  no spatial_neighbors_knn
   squidpy 1.8.0  no spatial_neighbors_knn
   squidpy 1.8.1  no spatial_neighbors_knn
   squidpy 1.8.2  spatial_neighbors_knn present
   ```

   The old `squidpy>=1.8` pin therefore permitted two versions that would raise
   `AttributeError` on every Nhood or Ligand-Receptor run. Both `environment.yml` and
   `pyproject.toml` now pin `squidpy>=1.8.2`, and
   `test_the_squidpy_floor_is_a_version_that_has_spatial_neighbors_knn` checks both the
   installed module and the declared floors so this cannot drift again. Verify a
   replacement against the wheel that is supposed to contain it, not against a
   `deprecated::` directive — the directive dates the *deprecation*, not the replacement.
3. **The `cnv` extra contains no squidpy**, contrary to the old checklist. Nothing to
   change there.
4. **`tests/test_notebook_replay.py` does cover this path** — it runs a
   `spatial_neighbors` step and an `nhood:` step and replays both in a clean kernel. The
   note claimed it did not, and sent the reader to `scripts/verify_notebook.py` as the
   only real gate.

## Known limitation

A notebook **already exported to disk** keeps the removed call: a provenance node stores
the code verbatim from when it ran, so re-exporting from the stored graph reproduces the
old cell. This self-heals rather than needing a migration — the next session that opens
Neighbourhood Enrichment or Ligand-Receptor re-runs `ensure_spatial_neighbors`, `upsert`
revises the node with the new code, and its descendants are flagged stale, which is
exactly what a changed upstream step is supposed to do.

Users with a customised `spatial_neighbors` template are flagged automatically:
`overrides.json` records the hash of the shipped block they forked, so the Templates tab
shows **⚠ review** and offers "Take new default".

Nothing on disk needed migrating — the graph lives in `obsp` and is rebuilt per session.

## Related

- `docs/pyqt6-migration.md` — the other tracked upstream deprecation (napari drops PyQt5
  in fall 2026). Note that its own inventory is understated; see the issue.
