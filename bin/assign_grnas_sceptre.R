#!/usr/bin/env Rscript
library(Matrix)

assign_grnas_sceptre_v1 <- function(mudata, probability_threshold = "default", n_em_rep = "default", n_processors = NA) {
  # convert MuData object to sceptre object, removing multicollinear covariates
  sceptre_object <- convert_mudata_to_sceptre_object_v1(
    mudata,
    remove_collinear_covariates = TRUE
  )

  # set analysis parameters
  sceptre_object <- sceptre_object |>
    sceptre::set_analysis_parameters(
      discovery_pairs = data.frame(
        grna_target = character(0),
        response_id = character(0)
      )
    )

  # assign gRNAs
  assign_grnas_args <- list(sceptre_object, method = "mixture")
  
  # Add optional parameters if they are not "default"
  if (probability_threshold != "default") {
    assign_grnas_args$probability_threshold <- as.numeric(probability_threshold)
  }
  if (n_em_rep != "default") {
    assign_grnas_args$n_em_rep <- as.integer(n_em_rep)
  }
  # Parallelization controls inferred from n_processors
  if (!is.na(n_processors) && as.integer(n_processors) > 1) {
    assign_grnas_args$parallel <- TRUE
    assign_grnas_args$n_processors <- as.integer(n_processors)
  }
  
  sceptre_object <- do.call(sceptre::assign_grnas, assign_grnas_args)

  # extract sparse logical matrix of gRNA assignments
  grna_assignment_matrix <- sceptre_object |>
    sceptre::get_grna_assignments() |>
    methods::as("dsparseMatrix")
  colnames(grna_assignment_matrix) <- colnames(MultiAssayExperiment::assay(mudata[['guide']]))

  # add gRNA assignment matrix to MuData
  SummarizedExperiment::assays(mudata[['guide']])[['guide_assignment']] <- grna_assignment_matrix

  return(list(mudata = mudata, guide_assignment = grna_assignment_matrix))
}

convert_mudata_to_sceptre_object_v1 <- function(mudata, remove_collinear_covariates = FALSE){
  # extract information from MuData
  moi <- MultiAssayExperiment::metadata(mudata[['guide']])$moi
  if(is.null(SummarizedExperiment::assayNames(mudata[['gene']]))){
    SummarizedExperiment::assayNames(mudata[['gene']]) <- 'counts'
  } else{
    SummarizedExperiment::assayNames(mudata[['gene']])[[1]] <- 'counts'
  }
  if(is.null(SummarizedExperiment::assayNames(mudata[['guide']]))){
    SummarizedExperiment::assayNames(mudata[['guide']]) <- 'counts'
  } else{
    SummarizedExperiment::assayNames(mudata[['guide']])[[1]] <- 'counts'
  }

  scRNA_data <- mudata@ExperimentList$gene
  guides_data <- mudata@ExperimentList$guide
  response_matrix <- scRNA_data@assays@data@listData[["counts"]]

  if(!is.null(SummarizedExperiment::colData(mudata))){
    covariates <- SummarizedExperiment::colData(mudata) |> as.data.frame()
    covariates[] <- lapply(covariates, as.factor)
    # Check and remove factor variables with fewer than two levels
    number_of_levels <- sapply(covariates, function(x) length(unique(x)))
    multi_level_factors <- number_of_levels > 1
    covariates_clean <- covariates[, multi_level_factors, drop = FALSE]

    if(ncol(covariates_clean) == 0){
      remove_collinear_covariates <- FALSE
    }
    
    if(remove_collinear_covariates){
      model_matrix <- stats::model.matrix(object = ~ ., data = covariates_clean)
      multicollinear <- Matrix::rankMatrix(model_matrix) < ncol(model_matrix)
      if(multicollinear){
        extra_covariates <- data.frame()
      } else{
        extra_covariates <- covariates_clean
      }
    } else{
      extra_covariates <- covariates_clean
    }
  } else{
    extra_covariates <- data.frame()
  }

  # if guide assignments not present, then extract guide counts
  if(length(guides_data@assays@data@listData) == 1){
    grna_matrix <- guides_data@assays@data@listData[["counts"]]
    # otherwise, extract guide assignments
  } else{
    grna_matrix <- guides_data@assays@data@listData[["guide_assignment"]]
  }

  grna_ids <- rownames(SingleCellExperiment::rowData(mudata[['guide']]))
  rownames(grna_matrix) <- grna_ids

  gene_ids <- rownames(SingleCellExperiment::rowData(mudata[['gene']]))
  rownames(response_matrix) <- gene_ids
  grna_target_data_frame <- SingleCellExperiment::rowData(mudata[['guide']]) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "grna_id") |>
    dplyr::rename(grna_target = intended_target_name) |>
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

# Run Command (only execute when script is run directly, not when sourced)
if (!exists(".sourced_from_test")) {
  
  args <- commandArgs(trailingOnly = TRUE)
  
  mudata_input <- args[1]
  probability_threshold <- if(length(args) >= 2) args[2] else "default"
  n_em_rep <- if(length(args) >= 3) args[3] else "default"
  n_processors <- if(length(args) >= 4) suppressWarnings(as.integer(args[4])) else NA
  
  mudata_in <- MuData::readH5MU(mudata_input)
  results <- assign_grnas_sceptre_v1(mudata = mudata_in, 
                                     probability_threshold = probability_threshold,
                                     n_em_rep = n_em_rep,
                                     n_processors = n_processors)
  guide_assignment <- results$guide_assignment
  
  writeMM(guide_assignment, "guide_assignment.mtx")
}
