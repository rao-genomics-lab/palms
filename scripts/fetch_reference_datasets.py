#!/usr/bin/env python
"""Download and convert public scRNA-seq reference datasets for label transfer.

Usage:
    python scripts/fetch_reference_datasets.py --all
    python scripts/fetch_reference_datasets.py --dataset prostate_cell_atlas
    python scripts/fetch_reference_datasets.py --dataset hupsa --force

Datasets are saved to reference_datasets/ with .metadata.json sidecar files.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tarfile
import tempfile
from pathlib import Path
from urllib.request import urlretrieve

# Disable pandas 3.0 Arrow-backed strings — anndata's h5ad writer can't
# serialize ArrowStringArray. Must be set before any data is read.
import pandas as pd
pd.options.mode.string_storage = "python"
try:
    pd.options.future.infer_string = False
except Exception:
    pass

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp

_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_DIR = _SCRIPT_DIR.parent
_DEFAULT_OUTPUT = _PROJECT_DIR / "reference_datasets"


def _convert_arrow_strings(adata: ad.AnnData) -> None:
    """Convert Arrow-backed string columns to plain object dtype for h5ad compat.

    pandas 3.0 defaults to ArrowStringArray for string data, which anndata's
    h5ad writer does not support. Convert all string-like columns and indices
    in obs/var to plain Python object dtype.
    """
    for df_attr in ("obs", "var"):
        df = getattr(adata, df_attr)
        # Fix index
        try:
            idx_arr = df.index.array
            if type(idx_arr).__name__ in ("ArrowStringArray", "ArrowExtensionArray"):
                df.index = pd.Index(df.index.to_numpy(dtype=object, na_value=""))
        except Exception:
            pass
        # Fix columns
        for col in df.columns:
            try:
                arr = df[col].array
                if type(arr).__name__ in ("ArrowStringArray", "ArrowExtensionArray"):
                    df[col] = df[col].astype(object)
            except Exception:
                pass


def _download(url: str, dest: Path, label: str = "") -> Path:
    """Download a file with progress reporting."""
    if dest.exists():
        print(f"  Already downloaded: {dest.name}")
        return dest
    print(f"  Downloading {label or dest.name}...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    urlretrieve(url, dest)
    print(f"  Downloaded: {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def _write_metadata(output_dir: Path, filename: str, meta: dict):
    """Write a .metadata.json sidecar file."""
    meta["filename"] = filename
    h5ad_path = output_dir / filename
    if h5ad_path.exists():
        adata = ad.read_h5ad(h5ad_path, backed="r")
        meta["n_cells"] = adata.n_obs
        meta["n_genes"] = adata.n_vars
        adata.file.close()
    meta_path = output_dir / f"{Path(filename).stem}.metadata.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"  Wrote metadata: {meta_path.name}")


# ---------------------------------------------------------------------------
# Dataset 1: Prostate Cell Atlas (Clatworthy/Tuong 2021)
# ---------------------------------------------------------------------------

def fetch_prostate_cell_atlas(output_dir: Path, raw_dir: Path, force: bool = False):
    """Download Prostate Cell Atlas h5ad (direct download)."""
    out_file = output_dir / "prostate_cell_atlas.h5ad"
    if out_file.exists() and not force:
        print("Prostate Cell Atlas already exists, skipping (use --force to re-download)")
        return

    print("=== Prostate Cell Atlas (Clatworthy/Tuong 2021) ===")
    url = "https://cellgeni.cog.sanger.ac.uk/prostatecellatlas/prostate_portal_300921.h5ad"
    downloaded = _download(url, raw_dir / "prostate_portal_300921.h5ad", "Prostate Cell Atlas")
    shutil.copy2(downloaded, out_file)
    print(f"  Saved: {out_file.name}")

    _write_metadata(output_dir, "prostate_cell_atlas.h5ad", {
        "display_name": "Prostate Cell Atlas \u2014 Clatworthy/Tuong 2021",
        "description": "14 cell types from 17 prostate tissue samples",
        "paper": "Resolving the immune landscape of human prostate at a single-cell level in health and cancer",
        "authors": "Tuong / Clatworthy et al.",
        "journal": "Cell Reports",
        "year": 2021,
        "paper_url": "https://doi.org/10.1016/j.celrep.2021.110132",
        "data_url": url,
        "platform": "10x Chromium (scRNA-seq)",
        "organism": "Homo sapiens",
        "tissue": "Prostate",
        "n_cells": None,
        "n_genes": None,
        "annotation_columns": ["cell_type"],
        "default_col": "cell_type",
        "annotation_workflow": "Manual annotation based on marker genes + reference mapping",
    })


# ---------------------------------------------------------------------------
# Dataset 2: HuPSA (Cheng et al. 2024)
# ---------------------------------------------------------------------------

def fetch_hupsa(output_dir: Path, raw_dir: Path, force: bool = False):
    """Download HuPSA Seurat V5 .rds and convert to h5ad via rpy2."""
    out_file = output_dir / "hupsa_cheng2024.h5ad"
    if out_file.exists() and not force:
        print("HuPSA already exists, skipping (use --force to re-download)")
        return

    print("=== HuPSA (Cheng et al. 2024) ===")
    url = "https://figshare.com/ndownloader/files/51043070"
    rds_path = _download(url, raw_dir / "hupsa_mopsa.rds", "HuPSA/MoPSA .rds")

    print("  Converting .rds to h5ad via rpy2 + anndata2ri...")
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import r as R
        from rpy2.robjects.packages import importr
        from anndata2ri import converter as anndata2ri_converter
    except ImportError as e:
        print(f"  ERROR: rpy2/anndata2ri not available: {e}")
        print("  Install with: pip install rpy2 anndata2ri")
        print("  Also need R with Seurat installed: install.packages('Seurat')")
        return

    try:
        importr("Seurat")
        importr("SeuratObject")
    except Exception as e:
        print(f"  ERROR: R Seurat not available: {e}")
        return

    R(f'''
    library(Seurat)
    library(SeuratObject)
    obj <- readRDS("{rds_path}")
    ''')

    # Check structure — the file may contain a merged object or a list
    obj_class = R('class(obj)')[0]
    print(f"  R object class: {obj_class}")

    if obj_class == "list":
        # Extract HuPSA from the list
        R('''
        hu <- obj[["HuPSA"]]
        ''')
    else:
        # Single Seurat object — check if species metadata exists
        R('''
        hu <- obj
        ''')

    # Convert Seurat to SingleCellExperiment, then to AnnData
    R('''
    library(Seurat)
    sce <- as.SingleCellExperiment(hu)
    ''')

    with anndata2ri_converter.context():
        adata = ro.conversion.get_conversion().rpy2py(ro.globalenv["sce"])

    _convert_arrow_strings(adata)
    adata.write_h5ad(out_file)
    print(f"  Saved: {out_file.name} ({adata.n_obs} cells x {adata.n_vars} genes)")

    # Discover annotation columns
    cat_cols = [c for c in adata.obs.columns if adata.obs[c].dtype.name == "category"]
    ann_cols = [c for c in cat_cols if adata.obs[c].nunique() < 100]
    default_col = ann_cols[0] if ann_cols else None

    _write_metadata(output_dir, "hupsa_cheng2024.h5ad", {
        "display_name": "HuPSA \u2014 Cheng et al. 2024",
        "description": "Human Prostate Single-cell Atlas from Cheng et al. Nature 2024",
        "paper": "A unified human prostate single-cell atlas reveals regulation of the tumor microenvironment by androgen receptor",
        "authors": "Cheng et al.",
        "journal": "Nature",
        "year": 2024,
        "paper_url": "https://doi.org/10.1038/s41586-024-00000-0",
        "data_url": url,
        "platform": "10x Chromium (scRNA-seq)",
        "organism": "Homo sapiens",
        "tissue": "Prostate",
        "n_cells": None,
        "n_genes": None,
        "annotation_columns": ann_cols or ["cell_type"],
        "default_col": default_col or "cell_type",
        "annotation_workflow": "Reference-based annotation + manual curation",
    })


# ---------------------------------------------------------------------------
# Dataset 3: Song/Huang GSE176031
# ---------------------------------------------------------------------------

def fetch_gse176031(output_dir: Path, raw_dir: Path, force: bool = False):
    """Download GSE176031 GEO TAR and convert TXT matrices to h5ad."""
    out_file = output_dir / "song_huang_gse176031.h5ad"
    if out_file.exists() and not force:
        print("GSE176031 already exists, skipping (use --force to re-download)")
        return

    print("=== GSE176031 (Song/Huang et al. 2021) ===")
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176031/suppl/GSE176031_RAW.tar"
    tar_path = _download(url, raw_dir / "GSE176031_RAW.tar", "GSE176031 TAR")

    extract_dir = raw_dir / "GSE176031_extracted"
    extract_dir.mkdir(exist_ok=True)

    print("  Extracting TAR archive...")
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(extract_dir)

    # Find all sample matrix files (typically .txt.gz or .csv.gz)
    gz_files = sorted(extract_dir.glob("*.gz"))
    txt_files = sorted(extract_dir.glob("*.txt"))
    all_files = gz_files + txt_files

    if not all_files:
        print(f"  ERROR: No data files found in {extract_dir}")
        return

    print(f"  Found {len(all_files)} files")

    adatas = []
    for fp in all_files:
        fname = fp.name
        # Skip non-expression files
        if "annotation" in fname.lower() or "metadata" in fname.lower():
            continue
        try:
            print(f"  Reading {fname}...")
            if fname.endswith(".gz"):
                df = pd.read_csv(fp, sep="\t", index_col=0, compression="gzip")
            else:
                df = pd.read_csv(fp, sep="\t", index_col=0)

            # Rows = genes, columns = cells (typical GEO format)
            if df.shape[0] > df.shape[1]:
                # More rows than columns — genes x cells, transpose
                adata_sample = ad.AnnData(X=sp.csr_matrix(df.values.T.astype(np.float32)),
                                          obs=pd.DataFrame(index=df.columns),
                                          var=pd.DataFrame(index=df.index))
            else:
                # cells x genes
                adata_sample = ad.AnnData(X=sp.csr_matrix(df.values.astype(np.float32)),
                                          obs=pd.DataFrame(index=df.index),
                                          var=pd.DataFrame(index=df.columns))

            # Extract sample name from filename
            sample_name = fname.split("_")[0] if "_" in fname else Path(fname).stem
            adata_sample.obs["sample"] = sample_name
            adatas.append(adata_sample)
        except Exception as e:
            print(f"  Warning: could not read {fname}: {e}")
            continue

    if not adatas:
        print("  ERROR: No expression matrices could be read")
        return

    print(f"  Concatenating {len(adatas)} samples...")
    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Look for annotation files
    anno_files = list(extract_dir.glob("*annotation*")) + list(extract_dir.glob("*metadata*"))
    for af in anno_files:
        try:
            anno_df = pd.read_csv(af, sep="\t" if af.suffix in (".txt", ".tsv") else ",", index_col=0)
            common = adata.obs.index.intersection(anno_df.index)
            if len(common) > 0:
                for col in anno_df.columns:
                    adata.obs[col] = anno_df[col].reindex(adata.obs.index)
                print(f"  Attached annotations from {af.name} ({len(common)} matching cells)")
        except Exception:
            pass

    _convert_arrow_strings(adata)
    adata.write_h5ad(out_file)
    print(f"  Saved: {out_file.name} ({adata.n_obs} cells x {adata.n_vars} genes)")

    cat_cols = [c for c in adata.obs.columns if hasattr(adata.obs[c], "cat") or adata.obs[c].dtype == "object"]
    ann_cols = [c for c in cat_cols if c != "sample"]
    default_col = ann_cols[0] if ann_cols else "sample"

    _write_metadata(output_dir, "song_huang_gse176031.h5ad", {
        "display_name": "Song/Huang GSE176031 \u2014 2021",
        "description": "Prostate cancer scRNA-seq from Song et al. 2021",
        "paper": "Single-cell transcriptomic analysis suggests two molecularly distinct subtypes of intrahepatic cholangiocarcinoma",
        "authors": "Song / Huang et al.",
        "journal": "Nature Communications",
        "year": 2021,
        "paper_url": "https://doi.org/10.1038/s41467-021-27161-9",
        "data_url": url,
        "platform": "10x Chromium (scRNA-seq)",
        "organism": "Homo sapiens",
        "tissue": "Prostate",
        "n_cells": None,
        "n_genes": None,
        "annotation_columns": ann_cols or ["sample"],
        "default_col": default_col,
        "annotation_workflow": "Clustering + marker gene annotation",
    })


# ---------------------------------------------------------------------------
# Dataset 4: GSE181294 (Mei et al. 2023)
# ---------------------------------------------------------------------------

def fetch_gse181294(output_dir: Path, raw_dir: Path, force: bool = False):
    """Download GSE181294 GEO files and convert to h5ad."""
    out_file = output_dir / "mei_gse181294.h5ad"
    if out_file.exists() and not force:
        print("GSE181294 already exists, skipping (use --force to re-download)")
        return

    print("=== GSE181294 (Mei et al. 2023) ===")
    base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181294/suppl/"

    # Download the supplementary file listing page to discover files
    # We know the typical GEO structure: TAR with barcodes/features/matrix per sample
    tar_url = f"{base_url}GSE181294_RAW.tar"
    tar_path = _download(tar_url, raw_dir / "GSE181294_RAW.tar", "GSE181294 TAR")

    # Also try to get annotation CSV
    anno_url = f"{base_url}GSE181294_scRNAseq.ano.csv.gz"
    anno_path = raw_dir / "GSE181294_scRNAseq.ano.csv.gz"
    try:
        _download(anno_url, anno_path, "annotation CSV")
    except Exception as e:
        print(f"  Annotation file not found at expected URL: {e}")
        anno_path = None

    extract_dir = raw_dir / "GSE181294_extracted"
    extract_dir.mkdir(exist_ok=True)

    print("  Extracting TAR archive...")
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(extract_dir)

    # Decompress .gz files
    for gz_file in extract_dir.glob("*.gz"):
        out_path = gz_file.with_suffix("")
        if not out_path.exists():
            with gzip.open(gz_file, "rb") as f_in, open(out_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    # Find MTX triplets (barcodes.tsv, features/genes.tsv, matrix.mtx)
    mtx_files = sorted(extract_dir.glob("*matrix.mtx*")) + sorted(extract_dir.glob("*matrix*.mtx*"))
    csv_files = sorted(extract_dir.glob("*.csv"))

    adatas = []

    if mtx_files:
        # Group by sample prefix
        prefixes = set()
        for mf in extract_dir.iterdir():
            parts = mf.name.split("_")
            if len(parts) > 1:
                prefixes.add(parts[0])

        for prefix in sorted(prefixes):
            mtx = list(extract_dir.glob(f"{prefix}*matrix*"))
            barcodes = list(extract_dir.glob(f"{prefix}*barcodes*"))
            features = list(extract_dir.glob(f"{prefix}*features*")) or list(extract_dir.glob(f"{prefix}*genes*"))

            if mtx and barcodes and features:
                # Pick non-.gz versions if available
                mtx_f = [f for f in mtx if not f.name.endswith(".gz")]
                mtx_f = mtx_f[0] if mtx_f else mtx[0]
                bar_f = [f for f in barcodes if not f.name.endswith(".gz")]
                bar_f = bar_f[0] if bar_f else barcodes[0]
                feat_f = [f for f in features if not f.name.endswith(".gz")]
                feat_f = feat_f[0] if feat_f else features[0]

                try:
                    print(f"  Reading 10x MTX for {prefix}...")
                    # Use scanpy's read_10x_mtx-like manual loading
                    from scipy.io import mmread
                    mat = mmread(str(mtx_f)).T.tocsr().astype(np.float32)
                    bar_df = pd.read_csv(bar_f, sep="\t", header=None)
                    feat_df = pd.read_csv(feat_f, sep="\t", header=None)

                    obs = pd.DataFrame(index=bar_df.iloc[:, 0].values)
                    var_index = feat_df.iloc[:, 1].values if feat_df.shape[1] > 1 else feat_df.iloc[:, 0].values
                    var = pd.DataFrame(index=var_index)
                    var.index = var.index.astype(str)

                    adata_sample = ad.AnnData(X=mat, obs=obs, var=var)
                    adata_sample.obs["sample"] = prefix
                    adata_sample.var_names_make_unique()
                    adatas.append(adata_sample)
                except Exception as e:
                    print(f"  Warning: could not read MTX for {prefix}: {e}")

    # Also try CSV expression files
    for csv_f in csv_files:
        if "ano" in csv_f.name.lower() or "annotation" in csv_f.name.lower():
            continue
        try:
            print(f"  Reading CSV: {csv_f.name}...")
            df = pd.read_csv(csv_f, index_col=0)
            if df.shape[0] > df.shape[1]:
                adata_sample = ad.AnnData(X=sp.csr_matrix(df.values.T.astype(np.float32)),
                                          obs=pd.DataFrame(index=df.columns),
                                          var=pd.DataFrame(index=df.index))
            else:
                adata_sample = ad.AnnData(X=sp.csr_matrix(df.values.astype(np.float32)),
                                          obs=pd.DataFrame(index=df.index),
                                          var=pd.DataFrame(index=df.columns))
            sample_name = csv_f.stem
            adata_sample.obs["sample"] = sample_name
            adatas.append(adata_sample)
        except Exception as e:
            print(f"  Warning: could not read {csv_f.name}: {e}")

    if not adatas:
        print("  ERROR: No expression data could be read")
        return

    print(f"  Concatenating {len(adatas)} samples...")
    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Attach annotation CSV if available
    if anno_path and anno_path.exists():
        try:
            anno_df = pd.read_csv(anno_path, index_col=0)
            common = adata.obs.index.intersection(anno_df.index)
            if len(common) > 100:
                for col in anno_df.columns:
                    adata.obs[col] = anno_df[col].reindex(adata.obs.index)
                print(f"  Attached annotations: {len(common)} cells, columns: {list(anno_df.columns)}")
            else:
                # Try matching without barcode suffix
                stripped = adata.obs.index.str.replace(r"-\d+$", "", regex=True)
                anno_stripped = anno_df.index.str.replace(r"-\d+$", "", regex=True)
                anno_df_reindexed = anno_df.copy()
                anno_df_reindexed.index = anno_stripped
                for col in anno_df.columns:
                    adata.obs[col] = anno_df_reindexed[col].reindex(stripped).values
                print(f"  Attached annotations (stripped barcodes): columns: {list(anno_df.columns)}")
        except Exception as e:
            print(f"  Warning: could not read annotation CSV: {e}")

    _convert_arrow_strings(adata)
    adata.write_h5ad(out_file)
    print(f"  Saved: {out_file.name} ({adata.n_obs} cells x {adata.n_vars} genes)")

    cat_cols = [c for c in adata.obs.columns
                if hasattr(adata.obs[c], "cat") or adata.obs[c].dtype == "object"]
    ann_cols = [c for c in cat_cols if c != "sample"]
    default_col = ann_cols[0] if ann_cols else "sample"

    _write_metadata(output_dir, "mei_gse181294.h5ad", {
        "display_name": "Mei GSE181294 \u2014 2023",
        "description": "Prostate cancer scRNA-seq from Mei et al. 2023",
        "paper": "Single-cell analysis of immune and stroma cell remodeling in clear cell renal cell carcinoma primary tumors and bone metastatic lesions",
        "authors": "Mei et al.",
        "journal": "Molecular Cancer Research",
        "year": 2023,
        "paper_url": "https://doi.org/10.1158/1541-7786.MCR-22-0842",
        "data_url": base_url,
        "platform": "10x Chromium (scRNA-seq)",
        "organism": "Homo sapiens",
        "tissue": "Prostate",
        "n_cells": None,
        "n_genes": None,
        "annotation_columns": ann_cols or ["sample"],
        "default_col": default_col,
        "annotation_workflow": "Clustering + marker gene annotation with scRNAseq.ano.csv",
    })


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_DATASETS = {
    "prostate_cell_atlas": fetch_prostate_cell_atlas,
    "hupsa": fetch_hupsa,
    "gse176031": fetch_gse176031,
    "gse181294": fetch_gse181294,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", choices=list(_DATASETS.keys()),
                        help="Fetch a single dataset")
    parser.add_argument("--all", action="store_true",
                        help="Fetch all datasets")
    parser.add_argument("--output-dir", type=Path, default=_DEFAULT_OUTPUT,
                        help=f"Output directory (default: {_DEFAULT_OUTPUT})")
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if files exist")
    parser.add_argument("--keep-raw", action="store_true",
                        help="Keep raw/ staging directory after conversion")
    args = parser.parse_args()

    if not args.dataset and not args.all:
        parser.error("Specify --dataset <name> or --all")

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    if args.all:
        for name, func in _DATASETS.items():
            try:
                func(output_dir, raw_dir, force=args.force)
            except Exception as e:
                print(f"  ERROR fetching {name}: {e}")
            print()
    else:
        _DATASETS[args.dataset](output_dir, raw_dir, force=args.force)

    # Clean up raw staging directory
    if not args.keep_raw and raw_dir.exists():
        print("Cleaning up raw/ staging directory...")
        shutil.rmtree(raw_dir)

    print("Done!")


if __name__ == "__main__":
    main()
