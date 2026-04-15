#!/usr/bin/env Rscript

library(argparse)
library(anndataR)
library(Matrix)

# Parse command line arguments
parser <- ArgumentParser(description="OmniBenchmark module")

# Required by OmniBenchmark
parser$add_argument("--output_dir", dest="output_dir", type="character", required=TRUE,
                   help="Output directory for results")
parser$add_argument("--name", dest="name", type="character", required=TRUE,
                   help="Module name/identifier")

# Add your custom input arguments here
# Example:

parser$add_argument("--normalization_type", dest="normalization_type", type="character", help="Input file")

# parameter for 'input_h5' (comes from 1st stage)
parser$add_argument("--rawdata.h5ad", dest="input_h5", type="character", help="input file")

args <- parser$parse_args()

cat("Full command: ", paste0(commandArgs(), collapse = " "))
cat("output_dir:", args$output_dir, "\n")
cat("name:", args$name, "\n")
cat("normalization_type:", args$normalization_type, "\n")
cat("input_h5:", args$input_h5, "\n")

if (args$normalization_type == "seurat_log1pCP10k") {
  require(Seurat)
  so <- read_h5ad(args$input_h5, as = "Seurat")
  message("Running log1pCP10k")
  so <- NormalizeData(so, scale.factor = 10^4)
  d <- so[["RNA"]]$data
} else if (args$normalization_type == "scuttle_geometricSizeFactors") {
  require(SingleCellExperiment)
  require(scuttle)
  sce <- read_h5ad(args$input_h5, as = "SingleCellExperiment")
  sce <- logNormCounts(sce, geometricSizeFactors(sce))
  d <- Matrix(logcounts(sce), sparse = TRUE)
} else {
  errorCondition("incorrect 'normalization_type' specified")
}

output_file <- file.path(args$output_dir, paste0(args$name, "_normalized.mtx"))
writeMM(d, output_file)
