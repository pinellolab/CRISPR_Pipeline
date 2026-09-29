#!/usr/bin/env Rscript
#
# BarcodeRanks knee/inflection cell calling (DropletUtils), same method and package
# used by 1.perturb_pipeline_output/perturb_pipeline/analysis/scripts/dropletutils_r.R,
# but run on a plain per-barcode total-UMI-count vector instead of the full count
# matrix -- barcodeRanks() only ever uses colSums(m), so a 1-row matrix of totals is
# sufficient and avoids exporting the full gene x cell matrix from Python.
#
# The perturb_pipeline ran this per sub-library (~65-70k barcodes), where no real
# knee exists (see tf_perturb project memory); this script is meant to be called
# per day/rep tag (~800-900k barcodes), where the knee/inflection are well-defined.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
  library(argparse)
})

parser <- ArgumentParser(description = "BarcodeRanks knee/inflection from a total-UMI-counts vector")
parser$add_argument("--counts_file", required = TRUE, help = "One total UMI count per line, no header")
parser$add_argument("--output_file", required = TRUE, help = "Output tsv with knee and inflection thresholds")

args <- parser$parse_args()

totals <- scan(args$counts_file, what = numeric(), quiet = TRUE)
m <- Matrix(totals, nrow = 1, sparse = TRUE)

br_out <- barcodeRanks(m)
knee <- metadata(br_out)$knee
inflection <- metadata(br_out)$inflection

cat("n_barcodes:", length(totals), "\n")
cat("knee:", knee, "\n")
cat("inflection:", inflection, "\n")

write.table(
  data.frame(knee = knee, inflection = inflection),
  file = args$output_file, sep = "\t", row.names = FALSE, quote = FALSE
)
