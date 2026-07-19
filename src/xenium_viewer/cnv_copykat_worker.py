"""Detached subprocess worker for the CopyKAT CNV backend.

CopyKAT is slow (~2 h/sample), so the viewer runs it here as a *separate*,
session-detached OS process (``subprocess.Popen(..., start_new_session=True)``)
that keeps running even if the GUI is closed. Everything this worker produces is
self-sufficient on disk so the viewer can pick the result up either live (by
polling the done-file) or on its next launch:

* ``adata_cnv_cache_copykat.h5ad`` beside the zarr cache — the CNV-profile
  AnnData (``obsm["X_cnv"]`` + ``copykat_leiden_res*`` cluster columns +
  ``cnv_score``), the same shape the inferCNV path caches.
* ``cnv_copykat_result.json`` — a sidecar with the metadata the viewer needs to
  rebuild its registry entry without re-reading the zarr table.
* ``plots/cnv_heatmap_copykat_<key>.png`` / ``.pdf`` — the chromosome heatmap.
* ``plots/copykat_DONE.txt`` — a status/timestamp marker; its presence signals
  completion (``status: ok`` or ``status: error``).

Top-level imports are kept light (no napari/Qt); the heavy work happens inside
``main`` with local imports, mirroring ``utils/leiden_worker.py``.

Usage::

    python -m xenium_viewer.cnv_copykat_worker <input_h5ad> <params_json> <output_dir>
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path


def _write_done(done_file: Path, status: str, **extra) -> None:
    payload = {"status": status, "timestamp": datetime.now().isoformat(), **extra}
    done_file.parent.mkdir(parents=True, exist_ok=True)
    done_file.write_text(json.dumps(payload, indent=2))


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: python -m xenium_viewer.cnv_copykat_worker "
              "<input_h5ad> <params_json> <output_dir>", file=sys.stderr)
        return 2

    input_h5ad, params_json, _output_dir = argv
    params = json.loads(Path(params_json).read_text())
    done_file = Path(params["done_file"])
    running_marker = Path(params.get("running_marker", "")) if params.get("running_marker") else None

    try:
        if running_marker is not None:
            running_marker.parent.mkdir(parents=True, exist_ok=True)
            running_marker.write_text(datetime.now().isoformat())

        import re
        import anndata as ad
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
        from xenium_viewer.utils.cnv_analysis import run_cnv_pipeline, make_cnv_heatmap

        adata = ad.read_h5ad(input_h5ad)
        ref_col = params["reference_obs_col"]
        idx = adata.obs["cell_id"].values if "cell_id" in adata.obs.columns else adata.obs_names
        reference_series = pd.Series(adata.obs[ref_col].astype(str).values, index=idx)

        result = run_cnv_pipeline(
            adata,
            reference_series,
            params["reference_categories"],
            reference_clustering_name=params.get("reference_clustering_name", ""),
            n_neighbors=params.get("n_neighbors", 15),
            smoothing_neighbors=params.get("smoothing_neighbors", 20),
            window_size=params.get("window_size", 60),
            step=params.get("step", 10),
            lfc_clip=params.get("lfc_clip", 4.0),
            resolution=params.get("resolution", 0.2),
            analyze_categories=None,  # cell-type restriction already applied by the parent
            backend="copykat",
            copykat_output_dir=params.get("copykat_workdir"),
        )

        adata_cnv = result["adata_cnv"]
        cluster_key = result["cluster_key"]
        # Stash the per-cell CNV score into obs so it round-trips inside the h5ad.
        score = result["cnv_score"]
        if "cell_id" in adata_cnv.obs.columns:
            adata_cnv.obs["cnv_score"] = score.reindex(adata_cnv.obs["cell_id"].values).to_numpy()
        else:
            adata_cnv.obs["cnv_score"] = score.reindex(adata_cnv.obs_names).to_numpy()

        out_h5ad = Path(params["output_h5ad"])
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        try:
            from xenium_viewer.utils.adata_persistence import _convert_adata_arrow_strings
            _convert_adata_arrow_strings(adata_cnv)
        except Exception:
            pass
        adata_cnv.write_h5ad(out_h5ad)

        # Chromosome heatmap (fork settings, via make_cnv_heatmap).
        plots_dir = Path(params["plots_dir"])
        plots_dir.mkdir(parents=True, exist_ok=True)
        safe = re.sub(r"[^0-9A-Za-z._-]", "_", cluster_key)
        png = plots_dir / f"cnv_heatmap_copykat_{safe}.png"
        pdf = plots_dir / f"cnv_heatmap_copykat_{safe}.pdf"
        fig = make_cnv_heatmap(adata_cnv, cluster_key)
        fig.savefig(png, dpi=200, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)

        # run_cnv_pipeline returns a single cluster_key (one resolution per run);
        # the tab derives the accumulated cluster_keys list, so build it here too.
        cluster_keys = list(result.get("cluster_keys") or [cluster_key])
        resolutions = sorted({
            float(m.group(1)) for k in cluster_keys
            for m in [re.search(r"res([0-9.]+)$", str(k))] if m
        })
        sidecar = {
            "backend": "copykat",
            "cluster_key": cluster_key,
            "cluster_keys": cluster_keys,
            "resolutions": resolutions,
            "analyze_categories": list(params.get("analyze_categories", [])),
            "reference_categories": list(result["reference_categories"]),
            "reference_clustering_name": result.get("reference_clustering_name", ""),
            "reference_obs_key": result.get("reference_obs_key"),
            "params": dict(result["params"]),
            "n_cells": int(result.get("n_cells", adata_cnv.n_obs)),
            "n_genes_total": int(result.get("n_genes_total", adata_cnv.n_vars)),
            "n_genes_mapped": int(result.get("n_genes_mapped", adata_cnv.n_vars)),
            "n_windows": int(result.get("n_windows", 0)),
            "max_cells": params.get("max_cells"),
            "adata_cnv_path": str(out_h5ad),
            "heatmap_png": str(png),
            "heatmap_pdf": str(pdf),
        }
        Path(params["result_json"]).write_text(json.dumps(sidecar, indent=2))

        _write_done(done_file, "ok", cluster_keys=cluster_keys,
                    n_cells=int(result.get("n_cells", adata_cnv.n_obs)))
        return 0
    except Exception as e:  # noqa: BLE001 — the worker must always leave a done-file
        import traceback
        _write_done(done_file, "error", error=f"{type(e).__name__}: {e}",
                    traceback=traceback.format_exc())
        print(f"CopyKAT worker failed: {e}", file=sys.stderr)
        return 1
    finally:
        if running_marker is not None:
            try:
                running_marker.unlink(missing_ok=True)
            except Exception:
                pass


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
