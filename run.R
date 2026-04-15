#!/usr/bin/env Rscript

library(argparse)

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

args <- parser$parse_args()

cat("Full command: ", paste0(commandArgs(), collapse = " "))
cat("Output directory:", args$output_dir, "\n")
cat("Module name:", args$name, "\n")
cat("Normalization type:", args$name, "\n")


if (args$filter_type == "manual") {
  qc <- metadata(sce)$qc_thresholds
  mt_percent <- rna.qc.metrics$subsets$mt * 100
  keep <- rna.qc.metrics$detected >= qc[qc$metric == "nFeature", "min"] &
    rna.qc.metrics$detected <= qc[qc$metric == "nFeature", "max"] &
    mt_percent < qc[qc$metric == "percent.mt", "max"] &
    rna.qc.metrics$sum <= qc[qc$metric == "nCount", "max"]
} else if (args$filter_type == "scrapper-auto") {
  require(DelayedArray)
  rna.qc.thresholds <- suggestRnaQcThresholds(rna.qc.metrics)
  keep <- filterRnaQcMetrics(rna.qc.thresholds, rna.qc.metrics)
}

output_file <- file.path(args$output_dir, paste0(args$name, "_cellids.txt.gz"))
writeLines(colnames(sce)[keep], gzfile(output_file))
