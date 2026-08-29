#!/usr/bin/env Rscript
# extract_seurat_segmentation.R
#
# Stage 1: Extract custom cell segmentation from a Seurat v5 RDS file.
# Accesses object slots directly — does NOT require a working Seurat installation.
# Outputs intermediate files consumed by build_custom_segmentation.py (Stage 2).
#
# Usage:
#   Rscript extract_seurat_segmentation.R <rds_path> <output_dir> [<image_slot>]
#
# Arguments:
#   rds_path    Path to the Seurat RDS file (e.g. xenium_v3.rds)
#   output_dir  Directory to write intermediate files (created if absent)
#   image_slot  Name of the image slot containing the segmentation
#               (default: auto-detect first slot with a @segmentation child)
#
# Outputs written to <output_dir>/:
#   segmentation_polygons.csv   cell_id, x, y, vertex_i  (one row per polygon vertex)
#   counts.mtx                  Sparse gene-by-cell count matrix (Matrix Market)
#   genes.txt                   One gene name per line (row names of counts)
#   barcodes.txt                One barcode per line (col names of counts)
#   cell_metadata.csv           All @meta.data columns + cell_id (= rownames)

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript extract_seurat_segmentation.R <rds_path> <output_dir> [<image_slot>]\n")
  quit(status = 1)
}

rds_path   <- args[1]
output_dir <- args[2]
image_slot <- if (length(args) >= 3) args[3] else NULL

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
message(sprintf("Output directory: %s", output_dir))

# ── Load RDS ──────────────────────────────────────────────────────────────────
message("Loading RDS (slot access only)...")
obj <- readRDS(rds_path)
message(sprintf("Object class: %s", paste(class(obj), collapse = ", ")))

# ── Locate segmentation ────────────────────────────────────────────────────────
images <- obj@images
image_names <- names(images)
message(sprintf("Image slots: %s", paste(image_names, collapse = ", ")))

if (is.null(image_slot)) {
  # Auto-detect: first slot that contains polygons (handles old and v5 FOV styles)
  for (nm in image_names) {
    img <- images[[nm]]
    # Old style: direct @segmentation slot
    if ("segmentation" %in% slotNames(img)) {
      image_slot <- nm
      message(sprintf("Auto-detected image slot with @segmentation: %s", nm))
      break
    }
    # Seurat v5 FOV: polygons live in @boundaries (a named list of Segmentation objects)
    if ("boundaries" %in% slotNames(img)) {
      bnds <- tryCatch(img@boundaries, error = function(e) NULL)
      if (!is.null(bnds) && length(bnds) > 0) {
        bobj <- bnds[[1]]
        if (!is.null(bobj) && "polygons" %in% slotNames(bobj) && length(bobj@polygons) > 0) {
          image_slot <- nm
          message(sprintf("Auto-detected FOV image slot via @boundaries: %s", nm))
          break
        }
      }
    }
  }
}

if (is.null(image_slot)) {
  stop("Could not auto-detect an image slot with segmentation. Specify it as the third argument.")
}

seg <- images[[image_slot]]
message(sprintf("Image slot class: %s", paste(class(seg), collapse = ", ")))

# Unwrap to the SpatialPolygons-like object containing @polygons
if ("segmentation" %in% slotNames(seg)) {
  # Old Seurat style
  seg <- seg@segmentation
} else if ("boundaries" %in% slotNames(seg)) {
  # Seurat v5 FOV: @boundaries is a named list; prefer "cell", else first entry
  bnds <- seg@boundaries
  preferred <- c("cell", "segmentation")
  hit <- intersect(preferred, names(bnds))
  bname <- if (length(hit) > 0) hit[1] else names(bnds)[1]
  seg <- bnds[[bname]]
  message(sprintf("Using FOV@boundaries[['%s']]", bname))
}

if (!"polygons" %in% slotNames(seg)) {
  stop(sprintf(
    "Cannot find @polygons on object of class '%s'. ",
    paste(class(seg), collapse = ", "),
    "Try printing slotNames(obj@images[['%s']]) and report the output.", image_slot
  ))
}

n_cells <- length(seg@polygons)
message(sprintf("Found %d cell polygons in slot '%s'", n_cells, image_slot))

# ── Extract polygon vertices ──────────────────────────────────────────────────
message("Extracting polygon vertices (this may take a few minutes)...")

# Pre-allocate by estimating ~25 vertices per cell
est_rows <- n_cells * 25L
cell_ids_vec <- character(est_rows)
x_vec        <- numeric(est_rows)
y_vec        <- numeric(est_rows)
vi_vec       <- integer(est_rows)
write_pos    <- 1L

for (i in seq_len(n_cells)) {
  poly   <- seg@polygons[[i]]
  coords <- poly@Polygons[[1L]]@coords  # Nx2, cols x y
  nv     <- nrow(coords)
  end    <- write_pos + nv - 1L
  if (end > length(cell_ids_vec)) {
    # Grow vectors
    extra <- max(nv, est_rows)
    cell_ids_vec <- c(cell_ids_vec, character(extra))
    x_vec        <- c(x_vec,        numeric(extra))
    y_vec        <- c(y_vec,        numeric(extra))
    vi_vec       <- c(vi_vec,       integer(extra))
  }
  cell_ids_vec[write_pos:end] <- poly@ID
  x_vec[write_pos:end]        <- coords[, 1L]
  y_vec[write_pos:end]        <- coords[, 2L]
  vi_vec[write_pos:end]       <- 0L:(nv - 1L)
  write_pos <- end + 1L

  if (i %% 50000L == 0L) message(sprintf("  %d / %d cells processed", i, n_cells))
}

polys_df <- data.frame(
  cell_id  = cell_ids_vec[seq_len(write_pos - 1L)],
  x        = x_vec[seq_len(write_pos - 1L)],
  y        = y_vec[seq_len(write_pos - 1L)],
  vertex_i = vi_vec[seq_len(write_pos - 1L)],
  stringsAsFactors = FALSE
)

out_polys <- file.path(output_dir, "segmentation_polygons.csv")
write.csv(polys_df, out_polys, row.names = FALSE, quote = FALSE)
message(sprintf("Polygons saved: %s (%d rows)", out_polys, nrow(polys_df)))

# ── Extract counts matrix ─────────────────────────────────────────────────────
message("Saving counts matrix (Matrix Market format)...")

# Try Xenium assay first, fall back to first assay
assay_name <- if ("Xenium" %in% names(obj@assays)) "Xenium" else names(obj@assays)[1]
message(sprintf("Using assay: %s", assay_name))
counts <- obj@assays[[assay_name]]$counts  # genes × cells

writeMM(counts, file.path(output_dir, "counts.mtx"))
writeLines(rownames(counts), file.path(output_dir, "genes.txt"))
writeLines(colnames(counts), file.path(output_dir, "barcodes.txt"))
message(sprintf("Counts saved: %d genes x %d cells", nrow(counts), ncol(counts)))

# ── Extract metadata ──────────────────────────────────────────────────────────
message("Saving cell metadata...")
meta <- obj@meta.data
meta$cell_id <- rownames(meta)
write.csv(meta, file.path(output_dir, "cell_metadata.csv"), row.names = FALSE, quote = FALSE)
message(sprintf("Metadata saved: %d cells, %d columns", nrow(meta), ncol(meta)))

message("Stage 1 complete.")
message(sprintf("Next step: conda activate palms && python scripts/build_custom_segmentation.py <xenium_dir> %s", output_dir))
