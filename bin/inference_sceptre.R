#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tibble)
})

# Declare known globals used in dplyr pipelines to avoid R CMD check notes
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "grna_id",
      "guide_id",
      "response_id",
      "grna_target",
      "p_value",
      "log_2_fold_change",
      "gene_id",
      "intended_target_name",
      "intended_target_chr",
      "intended_target_start",
      "intended_target_end",
      "intended_target_key",
      "log2_fc"
    )
  )
}

normalize_covariate_column <- function(column) {
  if (is.factor(column)) {
    return(droplevels(column))
  }
  if (is.character(column)) {
    return(factor(column))
  }
  if (is.logical(column)) {
    return(factor(column, levels = c(FALSE, TRUE)))
  }
  if (is.numeric(column) || is.integer(column)) {
    return(as.numeric(column))
  }
  factor(as.character(column))
}

column_has_variation <- function(column) {
  non_missing <- column[!is.na(column)]
  if (length(non_missing) == 0) {
    return(FALSE)
  }
  length(unique(non_missing)) > 1
}

make_model_matrix_data <- function(covariates_df) {
  as.data.frame(
    lapply(
      covariates_df,
      function(column) {
        if (is.factor(column)) {
          # Keep missing values explicit for model-matrix construction.
          return(stats::relevel(addNA(column), ref = levels(addNA(column))[1]))
        }
        if (all(is.na(column))) {
          return(rep(0, length(column)))
        }
        replacement <- stats::median(column, na.rm = TRUE)
        column[is.na(column)] <- replacement
        column
      }
    )
  )
}

retain_non_collinear_covariates <- function(covariates_df) {
  if (ncol(covariates_df) <= 1) {
    return(covariates_df)
  }

  kept <- character(0)
  dropped <- character(0)
  for (covariate_name in names(covariates_df)) {
    candidate_names <- c(kept, covariate_name)
    candidate_df <- covariates_df[, candidate_names, drop = FALSE]
    candidate_model_df <- make_model_matrix_data(candidate_df)
    model_matrix <- stats::model.matrix(object = ~., data = candidate_model_df)
    is_full_rank <- as.integer(Matrix::rankMatrix(model_matrix)) == ncol(model_matrix)
    if (is_full_rank) {
      kept <- candidate_names
    } else {
      dropped <- c(dropped, covariate_name)
    }
  }

  if (length(dropped) > 0) {
    warning(sprintf("Dropping collinear covariates: %s", paste(dropped, collapse = ", ")))
  }
  covariates_df[, kept, drop = FALSE]
}

prepare_extra_covariates <- function(mudata, remove_collinear_covariates = TRUE) {
  if (is.null(SummarizedExperiment::colData(mudata))) {
    return(data.frame())
  }
  covariates <- SummarizedExperiment::colData(mudata) |> as.data.frame(stringsAsFactors = FALSE)
  if (ncol(covariates) == 0) {
    return(data.frame())
  }

  covariates[] <- lapply(covariates, normalize_covariate_column)
  varying_columns <- vapply(covariates, column_has_variation, logical(1))
  covariates <- covariates[, varying_columns, drop = FALSE]
  if (ncol(covariates) == 0) {
    return(data.frame())
  }

  if (!remove_collinear_covariates) {
    return(covariates)
  }
  retain_non_collinear_covariates(covariates)
}
convert_mudata_to_sceptre_object_v1 <- function(mudata, remove_collinear_covariates = FALSE) {
  # extract information from MuData
  moi <- MultiAssayExperiment::metadata(mudata[["guide"]])$moi
  if (is.null(SummarizedExperiment::assayNames(mudata[["gene"]]))) {
    SummarizedExperiment::assayNames(mudata[["gene"]]) <- "counts"
  } else {
    SummarizedExperiment::assayNames(mudata[["gene"]])[[1]] <- "counts"
  }
  if (is.null(SummarizedExperiment::assayNames(mudata[["guide"]]))) {
    SummarizedExperiment::assayNames(mudata[["guide"]]) <- "counts"
  } else {
    SummarizedExperiment::assayNames(mudata[["guide"]])[[1]] <- "counts"
  }

  scRNA_data <- mudata@ExperimentList$gene
  guides_data <- mudata@ExperimentList$guide
  response_matrix <- scRNA_data@assays@data@listData[["counts"]]

  extra_covariates <- prepare_extra_covariates(
    mudata = mudata,
    remove_collinear_covariates = remove_collinear_covariates
  )

  # if guide assignments not present, then extract guide counts
  if (length(guides_data@assays@data@listData) == 1) {
    grna_matrix <- guides_data@assays@data@listData[["counts"]]
    # otherwise, extract guide assignments
  } else {
    grna_matrix <- guides_data@assays@data@listData[["guide_assignment"]]
  }

  grna_ids <- rownames(SingleCellExperiment::rowData(mudata[["guide"]]))
  rownames(grna_matrix) <- grna_ids

  gene_ids <- rownames(SingleCellExperiment::rowData(mudata[["gene"]]))
  rownames(response_matrix) <- gene_ids
  guide_row_data <- SingleCellExperiment::rowData(mudata[["guide"]]) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "grna_id")

  required_cols <- c(
    "intended_target_key",
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end"
  )
  missing_cols <- setdiff(required_cols, colnames(guide_row_data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Guide rowData is missing required columns: %s. Upstream metadata preparation did not run.",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (any(guide_row_data$intended_target_name == "non-targeting", na.rm = TRUE)) {
    stop("Found guides assigned to exact 'non-targeting'. Expected bucketed non-targeting groups (e.g., non-targeting|1).")
  }

  if (any(is.na(guide_row_data$intended_target_key))) {
    stop("Guide rowData contains NA intended_target_key values.")
  }

  grna_target_data_frame <- guide_row_data |>
    dplyr::transmute(grna_id = grna_id, grna_target = intended_target_key)

  # assemble information into sceptre object
  sceptre_object <- sceptre::import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = moi,
    extra_covariates = extra_covariates
  )

  # return sceptre object
  return(sceptre_object)
}


inference_sceptre_m <- function(mudata, n_processors = NA, ...) {
  # convert MuData object to sceptre object
  sceptre_object <- convert_mudata_to_sceptre_object_v1(mudata, remove_collinear_covariates = TRUE)
  guide_lookup <- SingleCellExperiment::rowData(mudata[["guide"]]) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "grna_id")
  if (!"guide_id" %in% colnames(guide_lookup)) {
    guide_lookup$guide_id <- guide_lookup$grna_id
  }
  target_lookup <- guide_lookup |>
    dplyr::select(
      intended_target_key,
      intended_target_name,
      intended_target_chr,
      intended_target_start,
      intended_target_end
    ) |>
    dplyr::distinct()

  args_list <- list(...)
  requested_control_group <- if ("control_group" %in% names(args_list)) as.character(args_list$control_group)[1] else NA_character_
  if (!is.na(requested_control_group) && requested_control_group != "complement") {
    warning(sprintf("Overriding control_group='%s' with 'complement' for SCEPTRE DE analysis.", requested_control_group))
  }
  args_list$control_group <- "complement"

  # Check if pairs_to_test exists in metadata
  if (!is.null(MultiAssayExperiment::metadata(mudata)$pairs_to_test)) {
    pairs_to_test <- MultiAssayExperiment::metadata(mudata)$pairs_to_test |>
      as.data.frame()

    if (!"gene_id" %in% colnames(pairs_to_test) && "gene_name" %in% colnames(pairs_to_test)) {
      pairs_to_test <- pairs_to_test |>
        dplyr::rename(gene_id = gene_name)
    }
    if (!"gene_id" %in% colnames(pairs_to_test)) {
      stop("pairs_to_test metadata must contain gene_id or gene_name.")
    }

    if (!"intended_target_key" %in% colnames(pairs_to_test)) {
      if (!"guide_id" %in% colnames(pairs_to_test)) {
        stop("pairs_to_test metadata is missing intended_target_key and guide_id columns.")
      }
      pairs_to_test <- pairs_to_test |>
        dplyr::mutate(guide_id = as.character(guide_id))
      pairs_to_test <- pairs_to_test |>
        dplyr::left_join(
          guide_lookup |>
            dplyr::transmute(guide_id = as.character(guide_id), intended_target_key = intended_target_key) |>
            dplyr::distinct(),
          by = "guide_id"
        )
    }

    discovery_pairs <- pairs_to_test |>
      dplyr::transmute(
        grna_target = as.character(intended_target_key),
        response_id = as.character(gene_id)
      ) |>
      dplyr::mutate(
        grna_target = replace(grna_target, tolower(trimws(grna_target)) == "nan", NA_character_)
      ) |>
      dplyr::filter(!is.na(grna_target), !is.na(response_id)) |>
      dplyr::distinct()

    args_list[["discovery_pairs"]] <- discovery_pairs
  } else {
    # No pairs_to_test found - use SCEPTRE's construct_trans_pairs for trans analysis
    args_list[["discovery_pairs"]] <- sceptre::construct_trans_pairs(sceptre_object)
  }

  # construct formula excluding gRNA covariates to avoid multicollinearity
  # (gRNA assignments are binary, making grna_n_nonzero and grna_n_umis identical)
  formula_object <- sceptre:::auto_construct_formula_object(
    cell_covariates = sceptre_object@covariate_data_frame,
    include_grna_covariates = FALSE
  )
  args_list[["formula_object"]] <- formula_object

  # We'll run two analyses on the same sceptre_object to exploit caching:
  # 1) union (grouped / per-element)
  # 2) singleton (per-guide)

  # Prepare a copy of the sceptre_object reference for chaining
  args_list$sceptre_object <- sceptre_object

  ## 1) Union / grouped (per-element) analysis
  args_union <- args_list
  args_union$grna_integration_strategy <- "union"
  # set analysis parameters for union
  sceptre_object <- do.call(sceptre::set_analysis_parameters, args_union)
  # assign grnas and run QC (relaxed thresholds to keep all cells; mirror prior behaviour)
  sceptre_object <- sceptre_object |>
    sceptre::assign_grnas(
      method = "thresholding",
      threshold = 1,
      parallel = (!is.na(n_processors) && as.integer(n_processors) > 1),
      n_processors = if (!is.na(n_processors) && as.integer(n_processors) > 1) as.integer(n_processors) else NULL
    ) |>
    sceptre::run_qc(
      n_nonzero_trt_thresh = 0L,
      n_nonzero_cntrl_thresh = 0L,
      p_mito_threshold = 1
    ) |>
    sceptre::run_discovery_analysis(
      parallel = (!is.na(n_processors) && as.integer(n_processors) > 1),
      n_processors = if (!is.na(n_processors) && as.integer(n_processors) > 1) as.integer(n_processors) else NULL
    )

  # get union (per-element) results
  union_results <- sceptre_object |>
    sceptre::get_result(analysis = "run_discovery_analysis") |>
    dplyr::select(response_id, grna_target, p_value, log_2_fold_change) |>
    dplyr::rename(
      gene_id = response_id,
      intended_target_key = grna_target,
      log2_fc = log_2_fold_change
    ) |>
    dplyr::left_join(target_lookup, by = "intended_target_key") |>
    dplyr::select(
      gene_id,
      intended_target_name,
      intended_target_chr,
      intended_target_start,
      intended_target_end,
      log2_fc,
      p_value
    )
  if (any(is.na(union_results$intended_target_name))) {
    stop("Unable to decode intended_target_key values back to intended target metadata in SCEPTRE union output.")
  }
  union_test_results <- union_results

  # store union results in mudata metadata as requested
  MultiAssayExperiment::metadata(mudata)$per_element_results <- union_test_results

  # also write the per-element (union) results to file
  try(
    write.table(
      union_test_results,
      file = gzfile("per_element_output.tsv.gz"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    ),
    silent = TRUE
  )

  ## 2) Singleton (per-guide) analysis — reuse sceptre_object to exploit caching
  args_singleton <- args_list
  args_singleton$grna_integration_strategy <- "singleton"
  # set analysis parameters for singleton (reuse same sceptre_object reference)
  sceptre_object <- do.call(sceptre::set_analysis_parameters, args_singleton)

  # run assignment, qc and discovery for singleton
  sceptre_object <- sceptre_object |>
    sceptre::assign_grnas(
      method = "thresholding",
      threshold = 1,
      parallel = (!is.na(n_processors) && as.integer(n_processors) > 1),
      n_processors = if (!is.na(n_processors) && as.integer(n_processors) > 1) as.integer(n_processors) else NULL
    ) |>
    sceptre::run_qc(
      n_nonzero_trt_thresh = 0L,
      n_nonzero_cntrl_thresh = 0L,
      p_mito_threshold = 1
    ) |>
    sceptre::run_discovery_analysis(
      parallel = (!is.na(n_processors) && as.integer(n_processors) > 1),
      n_processors = if (!is.na(n_processors) && as.integer(n_processors) > 1) as.integer(n_processors) else NULL
    )

  # extract singleton (per-guide) results, preserve grna_id and rename to guide_id
  singleton_results <- sceptre_object |>
    sceptre::get_result(analysis = "run_discovery_analysis") |>
    dplyr::select(response_id, grna_id, p_value, log_2_fold_change) |>
    dplyr::rename(
      gene_id = response_id,
      guide_id = grna_id,
      log2_fc = log_2_fold_change
    )
  singleton_test_results <- singleton_results
  MultiAssayExperiment::metadata(mudata)$per_guide_results <- singleton_test_results

  try(
    write.table(
      singleton_test_results,
      file = gzfile("per_guide_output.tsv.gz"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    ),
    silent = TRUE
  )

  # return mudata (with union results in metadata) and both result tables
  return(list(
    mudata = mudata,
    union_test_results = union_test_results,
    singleton_test_results = singleton_test_results
  ))
}

### Run Command (only execute when script is run directly, not when sourced)
if (!exists(".sourced_from_test")) {
  args <- commandArgs(trailingOnly = TRUE)
  args <- readLines(commandArgs(trailingOnly = TRUE)[1])

  # obtain the command line arguments
  mudata_fp <- args[1]
  side <- args[2]
  grna_integration_strategy <- args[3]
  resampling_approximation <- args[4]
  control_group <- args[5]
  resampling_mechanism <- args[6]
  n_processors <- if (length(args) >= 7) suppressWarnings(as.integer(args[7])) else NA

  # read MuData
  mudata_in <- MuData::readH5MU(mudata_fp)

  # run sceptre inference
  results <- inference_sceptre_m(
    mudata = mudata_in,
    side = side,
    grna_integration_strategy = grna_integration_strategy,
    resampling_approximation = resampling_approximation,
    control_group = control_group,
    resampling_mechanism = resampling_mechanism,
    n_processors = n_processors
  )
  # write outputs: per-element (union) and per-guide (singleton)
  if (!is.null(results$union_test_results)) {
    try(write.table(results$union_test_results, file = "per_element_output.tsv", sep = "\t", row.names = FALSE, quote = FALSE), silent = TRUE)
  }
  if (!is.null(results$singleton_test_results)) {
    try(write.table(results$singleton_test_results, file = "per_guide_output.tsv", sep = "\t", row.names = FALSE, quote = FALSE), silent = TRUE)
  }

  # write the modified MuData (contains union results in metadata as 'test_results')
  try(MuData::writeH5MU(object = results$mudata, file = "inference_mudata.h5mu"), silent = TRUE)
}
