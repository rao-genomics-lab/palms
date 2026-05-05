"""Per-tab UI modules for the Xenium Viewer control panel."""

from xenium_viewer.tabs.tab_clustering import build_tab as build_clustering_tab
from xenium_viewer.tabs.tab_cell_coloring import build_tab as build_cell_coloring_tab
from xenium_viewer.tabs.tab_transcripts import build_tab as build_transcripts_tab
from xenium_viewer.tabs.tab_umap import build_tab as build_umap_tab
from xenium_viewer.tabs.tab_roi import build_tab as build_roi_tab
from xenium_viewer.tabs.tab_he_registration import build_tab as build_he_registration_tab
from xenium_viewer.tabs.tab_gene_analysis import build_tab as build_gene_analysis_tab
from xenium_viewer.tabs.tab_ligrec import build_tab as build_ligrec_tab
from xenium_viewer.tabs.tab_nhood import build_tab as build_nhood_tab
from xenium_viewer.tabs.tab_co_occurrence import build_tab as build_co_occurrence_tab
from xenium_viewer.tabs.tab_arms import build_tab as build_arms_tab

__all__ = [
    "build_clustering_tab",
    "build_cell_coloring_tab",
    "build_transcripts_tab",
    "build_umap_tab",
    "build_roi_tab",
    "build_he_registration_tab",
    "build_gene_analysis_tab",
    "build_ligrec_tab",
    "build_nhood_tab",
    "build_co_occurrence_tab",
    "build_arms_tab",
]
