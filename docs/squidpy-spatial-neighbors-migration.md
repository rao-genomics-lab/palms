# TODO: `sq.gr.spatial_neighbors` is removed in squidpy 1.9 — migrate before upgrading

Status: **not started / tracked** (as of 2026-08-17). Priority: **do this before any
squidpy upgrade past 1.8.x.** No code changes made yet.

## Why this matters more than a usual deprecation

squidpy 1.8.2 prints, on every neighbour-graph build:

> `FutureWarning: Calling 'spatial_neighbors' is deprecated and will be removed in squidpy
> v1.9.0. Use 'spatial_neighbors_knn', 'spatial_neighbors_radius',
> 'spatial_neighbors_delaunay', 'spatial_neighbors_grid', or
> 'spatial_neighbors_from_builder' instead.`

The spatial neighbour graph is a **dependency of three analyses**, not a leaf feature.
`ctx.ensure_spatial_neighbors()` builds it, and Neighbourhood Enrichment, Co-occurrence and
Ligand-Receptor all declare `deps=["spatial_neighbors"]`. When squidpy 1.9 lands, all three
tabs stop working at once — and so does **every exported notebook**, including ones already
written, because the recorded `spatial_neighbors` cell calls the removed function. That is
the part that makes this worth doing early: the provenance graph's promise is that an old
notebook still replays.

`environment.yml` and `pyproject.toml` both pin `squidpy>=1.8` with no upper bound, so a
routine `conda env update` is enough to break it. Nothing here is failing today.

## The good news: measured, the swap is result-preserving

The concern that would make this risky — a different graph, hence different enrichment
z-scores and invalidated saved results — **does not apply to how this codebase calls it.**
Measured on 500 random points, squidpy 1.8.2:

| | `spatial_neighbors(coord_type='generic', n_neighs=6)` | `spatial_neighbors_knn(n_neighs=6)` |
|---|---|---|
| `obsp` keys | `spatial_connectivities`, `spatial_distances` | same |
| `uns` keys | `spatial_neighbors` | same |
| `spatial_connectivities` | 3000 nnz | 3000 nnz, **identical matrix** |
| `spatial_distances` | 3000 nnz | 3000 nnz, **identical matrix** |

`sq.gr.nhood_enrichment` then runs unchanged against the new graph. So for our one call
shape this is a rename, not a behaviour change: `coord_type='generic'` with `n_neighs=k`
*is* k-nearest-neighbours, and that is exactly what `spatial_neighbors_knn` does.

Reproduce with the comparison in the "Verification" section below before trusting this —
it was measured against 1.8.2 and the replacement's behaviour could still shift before 1.9.

## Call sites

Only two, and they must change together — the template is what the notebook replays, the
function is what the GUI runs:

1. **`src/xenium_viewer/utils/spatial_analysis.py:35`** — `compute_spatial_neighbors()`,
   used directly by `tabs/tab_annot_nhood.py:175`.
2. **`src/xenium_viewer/utils/step_templates/builtin/spatial_neighbors.tmpl:15`** — the
   recorded code, reached through `ctx.ensure_spatial_neighbors()` in `tabs/_helpers.py`.

**Not a call site:** `tabs/tab_novae.py` calls `novae.spatial_neighbors()`, which is the
Novae library's own function and has nothing to do with squidpy. Leave it alone.

## Checklist

1. **Confirm the version floor.** `spatial_neighbors_knn` exists in 1.8.2; check whether it
   exists in **1.8.0**, the current floor in `environment.yml` and `pyproject.toml`. If not,
   raise the pin to the first version that has it in *both* files and in the `cnv` extra —
   otherwise the swap breaks anyone on an older 1.8.x.
2. **Re-run the equivalence check** below against whatever squidpy version is current.
3. **Change the template first**, since it is the harder half:
   ```python
   # was
   sq.gr.spatial_neighbors(adata_norm, coord_type='generic', n_neighs=$n_neighs)
   # becomes
   sq.gr.spatial_neighbors_knn(adata_norm, n_neighs=$n_neighs)
   ```
   Note the replacements are **keyword-only** after the first argument, and `coord_type` is
   gone entirely. The `params` contract (`n_neighs:int`) is unchanged, so no
   `check_step` breakage — but `tests/test_template_registry.py` renders every template
   against every declared assembly and will catch a typo.
4. **Change `compute_spatial_neighbors()`** to match, keeping its signature
   (`adata_norm`, `n_neighs=6`) so `tab_annot_nhood.py` needs no edit.
5. **Users with a customised `spatial_neighbors` template will be flagged for review**
   automatically — `overrides.json` records the hash of the shipped block they forked, so
   the Templates tab shows **⚠ review** and offers "Take new default". That is the
   mechanism working as designed; no extra migration step is needed for them.
6. **Old caches keep working.** The graph lives in `obsp`, is rebuilt per session by
   `ensure_spatial_neighbors`, and is not persisted — so there is nothing on disk to migrate.

## Verification

- `pytest tests/test_notebook_replay.py` — the real gate. It executes the exported notebook
  in a clean kernel, so a removed-or-renamed squidpy call fails there rather than in a user's
  session. It does not currently cover the spatial tabs, so **also** run
  `scripts/verify_notebook.py <dataset>` against a dataset that has a recorded
  `spatial_neighbors` node.
- Run Neighbourhood Enrichment, Co-occurrence and Ligand-Receptor in the GUI and confirm
  the z-scores match what the dataset already has saved — identical, not merely similar,
  given the equivalence measured above.
- The equivalence check itself:
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

## Related

- `docs/pyqt6-migration.md` — the other tracked upstream deprecation (napari drops PyQt5 in
  fall 2026).
