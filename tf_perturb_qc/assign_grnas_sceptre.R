#!/usr/bin/env Rscript
#
# Guide-to-cell assignment (sceptre::assign_grnas(), mixture method) + optional calibration
# check and on-target discovery analysis, adapted from the lab's tuned implementation for this
# dataset (CRISPR_pipeline/bin/assign_grnas_sceptre.R, see conf/tf_2400_tri_sceptre.config:
# mixture method, probability_threshold/n_em_rep left at sceptre defaults, moi="high"). That
# version builds the sceptre_object from a MuData read via the R "MuData" package; op_pah has
# sceptre itself but not MuData/MultiAssayExperiment, so this version takes plain
# MatrixMarket/text inputs instead (same lightweight pattern as barcode_ranks_knee.R) and
# builds the sceptre_object directly -- avoids adding two more heavy Bioconductor dependencies
# for what's otherwise just a data-loading step.
#
# Named arguments (--flag value), not positional -- there are too many for positional args to
# stay readable. Required: --response_matrix --gene_ids --grna_matrix --grna_ids --grna_target
# --cell_barcodes --moi --output_mtx. Optional: --extra_covariates --discovery_pairs
# --probability_threshold --n_em_rep --save_plots --plot_dir --n_grnas_to_plot
# --n_calibration_pairs --calibration_group_size --run_discovery.
suppressPackageStartupMessages({
  library(Matrix)
})

parse_named_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    parsed[[key]] <- args[i + 1]
    i <- i + 2
  }
  parsed
}

a <- parse_named_args(commandArgs(trailingOnly = TRUE))
get_arg <- function(name, default = NULL) if (!is.null(a[[name]])) a[[name]] else default

moi <- get_arg("moi")
output_mtx_path <- get_arg("output_mtx")
probability_threshold <- get_arg("probability_threshold", "default")
n_em_rep <- get_arg("n_em_rep", "default")
save_plots <- as.logical(get_arg("save_plots", "FALSE"))
plot_dir <- get_arg("plot_dir", NA)
n_grnas_to_plot <- as.integer(get_arg("n_grnas_to_plot", 9))
n_calibration_pairs <- as.integer(get_arg("n_calibration_pairs", 5000))
# 3 guides/target (this is the "TF tri-guide" design) -- grouping negative-control guides the
# same way lets the calibration check's null test mimic the real positive-control test setup.
calibration_group_size <- as.integer(get_arg("calibration_group_size", 3))
run_discovery <- as.logical(get_arg("run_discovery", "FALSE"))

cat("Reading response (gene) matrix...\n")
response_matrix <- readMM(get_arg("response_matrix"))
cell_barcodes <- readLines(get_arg("cell_barcodes"))
rownames(response_matrix) <- readLines(get_arg("gene_ids"))
colnames(response_matrix) <- cell_barcodes

cat("Reading gRNA matrix...\n")
grna_matrix <- readMM(get_arg("grna_matrix"))
rownames(grna_matrix) <- readLines(get_arg("grna_ids"))
colnames(grna_matrix) <- cell_barcodes

grna_target_data_frame <- read.table(get_arg("grna_target"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# extra_covariates: user-requested (2026-09-04) sub-library, replicate, %ribo, %mt, on top of
# sceptre's own automatically-derived covariates (response/grna n_umis, n_nonzero, etc.) --
# these are per-cell technical covariates sceptre doesn't compute on its own.
extra_covariates_path <- get_arg("extra_covariates")
if (!is.null(extra_covariates_path)) {
  extra_covariates <- read.table(extra_covariates_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(extra_covariates) <- cell_barcodes
  extra_covariates$sub <- as.factor(extra_covariates$sub)
  extra_covariates$rep <- as.factor(extra_covariates$rep)
} else {
  extra_covariates <- data.frame()
}

cat("Building sceptre object (moi=", moi, ", n_cells=", length(cell_barcodes), ")...\n")
sceptre_object <- sceptre::import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = moi,
  extra_covariates = extra_covariates
)

# discovery_pairs: on-target only (each real TF target vs its own gene, ~2400 pairs for this
# dataset) -- a knockdown-efficacy check, not genome-wide trans discovery (2026-09-04 decision:
# genome-wide trans would be tens of millions of tests, out of scope for now). Falls back to
# empty (assignment/calibration-check-only run) if not supplied.
discovery_pairs_path <- get_arg("discovery_pairs")
discovery_pairs <- if (!is.null(discovery_pairs_path)) {
  read.table(discovery_pairs_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
} else {
  data.frame(grna_target = character(0), response_id = character(0))
}
cat("Discovery pairs:", nrow(discovery_pairs), "\n")

sceptre_object <- sceptre_object |>
  sceptre::set_analysis_parameters(discovery_pairs = discovery_pairs)

assign_grnas_args <- list(sceptre_object, method = "mixture")
if (probability_threshold != "default") {
  assign_grnas_args$probability_threshold <- as.numeric(probability_threshold)
}
if (n_em_rep != "default") {
  assign_grnas_args$n_em_rep <- as.integer(n_em_rep)
}

cat("Running assign_grnas (mixture)...\n")
sceptre_object <- do.call(sceptre::assign_grnas, assign_grnas_args)

grna_assignment_matrix <- sceptre_object |>
  sceptre::get_grna_assignments() |>
  methods::as("dsparseMatrix")

n_assigned <- sum(Matrix::colSums(grna_assignment_matrix) > 0)
cat("Cells with >=1 assigned guide:", n_assigned, "/", ncol(grna_assignment_matrix), "\n")

writeMM(grna_assignment_matrix, output_mtx_path)
cat("Wrote:", output_mtx_path, "\n")

if (isTRUE(save_plots)) {
  if (is.na(plot_dir)) stop("--save_plots requires --plot_dir to be set")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Calibration check: tests random genes against synthetic negative-control gRNA groups (built
  # from "non-targeting" guides, grouped calibration_group_size at a time to mimic the real
  # multi-guide-per-target design) -- p-values should be ~uniform if the testing procedure
  # itself isn't inflating false positives. Also auto-runs run_qc() the first time.
  cat("Running calibration check (n_calibration_pairs=", n_calibration_pairs,
      ", calibration_group_size=", calibration_group_size, ")...\n")
  sceptre_object <- sceptre::run_calibration_check(
    sceptre_object,
    n_calibration_pairs = n_calibration_pairs,
    calibration_group_size = calibration_group_size,
    parallel = FALSE
  )
  calibration_results <- sceptre::get_result(sceptre_object, "run_calibration_check")
  write.table(
    calibration_results, file.path(plot_dir, "run_calibration_check_results.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  if (isTRUE(run_discovery) && nrow(discovery_pairs) > 0) {
    cat("Running discovery analysis (", nrow(discovery_pairs), "on-target pairs)...\n")
    sceptre_object <- sceptre::run_discovery_analysis(sceptre_object, parallel = FALSE)
    discovery_results <- sceptre::get_result(sceptre_object, "run_discovery_analysis")
    write.table(
      discovery_results, file.path(plot_dir, "run_discovery_analysis_results.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }

  # Dumps the full standard sceptre output bundle (per-gRNA assignment diagnostic, gRNA count
  # distributions, QC plot, calibration check plot+results, discovery analysis plot+results,
  # analysis_summary.txt) -- covers everything above plus more, so no need to ggsave() each
  # plot manually here.
  cat("Writing full sceptre output bundle to", plot_dir, "...\n")
  sceptre::write_outputs_to_directory(sceptre_object, plot_dir)
  cat("Saved plots/tables to:", plot_dir, "\n")
}
