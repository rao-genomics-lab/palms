"""
Spatial analysis utilities for the Xenium viewer.

Provides:
  - Spatial neighbor graph construction via squidpy
  - Ligand-receptor interaction analysis
  - L-R result plotting
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import warnings


def compute_spatial_neighbors(
    adata_norm: sc.AnnData,
    n_neighs: int = 6,
) -> None:
    """Compute spatial neighbor graph (modifies adata_norm in-place)."""
    sq.gr.spatial_neighbors(
        adata_norm, coord_type='generic', n_neighs=n_neighs,
    )


def run_ligrec(
    adata_norm: sc.AnnData,
    cluster_key: str,
    n_perms: int = 1000,
    threshold: float = 0.01,
    seed: int = 42,
) -> dict:
    """Run ligand-receptor interaction analysis.

    Returns dict with keys 'means', 'pvalues' (both DataFrames),
    and 'warning' (str or None).
    """
    warning_msg = None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = sq.gr.ligrec(
                adata_norm,
                cluster_key=cluster_key,
                n_perms=n_perms,
                threshold=threshold,
                seed=seed,
                copy=True,
            )
        means = result['means']
        pvalues = result['pvalues']

        if means.shape[0] == 0:
            warning_msg = (
                "No ligand-receptor interactions found. "
                "The 480-gene Xenium panel may not contain enough L-R pairs."
            )
    except Exception as e:
        means = pd.DataFrame()
        pvalues = pd.DataFrame()
        warning_msg = f"L-R analysis failed: {e}"

    return {
        'means': means,
        'pvalues': pvalues,
        'warning': warning_msg,
    }


def make_ligrec_plot(
    result: dict,
    pvalue_threshold: float = 0.05,
) -> plt.Figure:
    """Create a ligand-receptor interaction dot plot. Returns matplotlib Figure."""
    sq.pl.ligrec(
        result,
        pvalue_threshold=pvalue_threshold,
        show=False,
    )
    return plt.gcf()
