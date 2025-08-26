#!/usr/bin/env Rscript

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
      "log2_fc"
    )
  )
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

  if (!is.null(SummarizedExperiment::colData(mudata))) {
    covariates <- SummarizedExperiment::colData(mudata) |> as.data.frame()
    covariates[] <- lapply(covariates, as.factor)
    # Check and remove factor variables with fewer than two levels
    number_of_levels <- sapply(covariates, function(x) length(unique(x)))
    multi_level_factors <- number_of_levels > 1
    covariates_clean <- covariates[, multi_level_factors, drop = FALSE]

    if (ncol(covariates_clean) == 0) {
      remove_collinear_covariates <- FALSE
    }

    if (remove_collinear_covariates) {
      model_matrix <- stats::model.matrix(object = ~., data = covariates_clean)
      multicollinear <- Matrix::rankMatrix(model_matrix) < ncol(model_matrix)
      if (multicollinear) {
        extra_covariates <- data.frame()
      } else {
        extra_covariates <- covariates_clean
      }
    } else {
      extra_covariates <- covariates_clean
    }
  } else {
    extra_covariates <- data.frame()
  }

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
  grna_target_data_frame <- SingleCellExperiment::rowData(mudata[["guide"]]) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "grna_id") |>
    dplyr::rename(grna_target = intended_target_name) |>
    # dplyr::mutate(grna_target = ifelse(!targeting & (moi=='low'), "non-targeting", grna_target)) |>
    dplyr::select(grna_id, grna_target)

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


inference_sceptre_m <- function(mudata, ...) {
  # convert MuData object to sceptre object
  sceptre_object <- convert_mudata_to_sceptre_object_v1(mudata, remove_collinear_covariates = TRUE)

  # extract set of discovery pairs to test
  pairs_to_test <- MultiAssayExperiment::metadata(mudata)$pairs_to_test |>
    as.data.frame()
  moi <- MultiAssayExperiment::metadata(mudata[["guide"]])$moi

  discovery_pairs <- pairs_to_test |>
    dplyr::rename(
      grna_target = intended_target_name,
      response_id = gene_id
    ) |>
    dplyr::filter(grna_target != "non-targeting")

  # assemble base arguments to set_analysis_parameters()
  args_list <- list(...)

  if ("discovery_pairs" %in% names(args_list)) {
    warning("The `discovery_pairs` argument is ignored. The `discovery_pairs` are set from the `pairs_to_test` metadata.")
  }
  # always use discovery pairs from metadata
  args_list[["discovery_pairs"]] <- discovery_pairs
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
    sceptre::assign_grnas(method = "thresholding", threshold = 1) |>
    sceptre::run_qc(
      n_nonzero_trt_thresh = 0L,
      n_nonzero_cntrl_thresh = 0L,
      p_mito_threshold = 1
    )

  # run discovery analysis (grouped)
  sceptre_object <- sceptre_object |>
    sceptre::run_discovery_analysis()

  # get union (per-element) results
  union_results <- sceptre_object |>
    sceptre::get_result(analysis = "run_discovery_analysis") |>
    dplyr::select(response_id, grna_target, p_value, log_2_fold_change) |>
    dplyr::rename(
      gene_id = response_id,
      intended_target_name = grna_target,
      log2_fc = log_2_fold_change
    )
  # use union_results directly (no extra distinct/left_join)
  union_test_results <- union_results

  # store union results in mudata metadata as requested
  MultiAssayExperiment::metadata(mudata)$test_results <- union_test_results

  # also write the per-element (union) results to file
  try(
    write.table(
      union_test_results,
      file = "per_element_output.tsv",
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
    sceptre::assign_grnas(method = "thresholding", threshold = 1) |>
    sceptre::run_qc(
      n_nonzero_trt_thresh = 0L,
      n_nonzero_cntrl_thresh = 0L,
      p_mito_threshold = 1
    ) |>
    sceptre::run_discovery_analysis()

  # extract singleton (per-guide) results, preserve grna_id and rename to guide_id
  singleton_results <- sceptre_object |>
    sceptre::get_result(analysis = "run_discovery_analysis") |>
    dplyr::select(response_id, grna_id, grna_target, p_value, log_2_fold_change) |>
    dplyr::rename(
      gene_id = response_id,
      guide_id = grna_id,
      intended_target_name = grna_target,
      log2_fc = log_2_fold_change
    )

  # use singleton_results directly
  singleton_test_results <- singleton_results

  try(
    write.table(
      singleton_test_results,
      file = "per_guide_output.tsv",
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

  # read MuData
  mudata_in <- MuData::readH5MU(mudata_fp)

  # run sceptre inference
  results <- inference_sceptre_m(
    mudata = mudata_in,
    side = side,
    grna_integration_strategy = grna_integration_strategy,
    resampling_approximation = resampling_approximation,
    control_group = control_group,
    resampling_mechanism = resampling_mechanism
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
