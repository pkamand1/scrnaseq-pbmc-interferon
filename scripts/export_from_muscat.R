# export_from_muscat.R
# Exports the Kang et al. (2018) PBMC subset that ships with the Bioconductor package muscat
# (object `example_sce`) to plain files that Python can read.
#
# Source data: Kang, H. M., et al. (2018). Multiplexed droplet single-cell RNA-sequencing using natural
# genetic variation. Nature Biotechnology, 36, 89-94. GEO accession GSE96583.
# Subset as distributed in muscat (Crowell et al., 2020, Nature Communications, 11, 6077).
#
# Requires: SingleCellExperiment, Matrix, and the muscat package (or its data/example_sce.rda file).
# Usage (from the repository root):  Rscript scripts/export_from_muscat.R [path/to/example_sce.rda]

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) load(args[1]) else data("example_sce", package = "muscat")

counts <- as(counts(example_sce), "CsparseMatrix")
dir.create("data", showWarnings = FALSE)
writeMM(counts, "data/counts.mtx")
system("gzip -f data/counts.mtx")

genes <- data.frame(symbol = rownames(example_sce), ensembl = rowData(example_sce)$ENSEMBL)
cells <- data.frame(barcode = colnames(example_sce),
                    sample = example_sce$sample_id,
                    donor = sub("^(ctrl|stim)", "donor_", example_sce$sample_id),
                    condition = ifelse(example_sce$group_id == "stim", "IFN-beta", "control"),
                    author_cell_type = example_sce$cluster_id)
write.csv(genes, "data/genes.csv", row.names = FALSE)
write.csv(cells, "data/cells.csv", row.names = FALSE)
cat(nrow(counts), "genes x", ncol(counts), "cells written to data/\n")
