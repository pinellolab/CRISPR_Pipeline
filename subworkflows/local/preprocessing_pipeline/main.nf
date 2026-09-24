include { PreprocessAnnData } from '../../../modules/local/PreprocessAnnData'
include { concat_preprocessed_rna } from '../../../modules/local/concat_preprocessed_rna'
include { collect_measurement_set_qc } from '../../../modules/local/collect_measurement_set_qc'
include { skipGTFDownload } from '../../../modules/local/skipGTFDownload'
include { downloadGTF } from '../../../modules/local/downloadGTF'

workflow preprocessing_pipeline {

    take:
    concat_anndata_rna
    trans_out_dir
    parsed_covariate_file

    main:
    // Convert the single emitted covariate file to a value channel. Otherwise
    // it would be consumed as a second queue input and only the first RNA
    // measurement set would launch a QC task.
    parsed_covariate_value = parsed_covariate_file.first()
    // HTO singlet selection is the mutually exclusive doublet/multiplet filter
    // for hashing-enabled experiments. Scrublet is only meaningful when HTO is
    // unavailable, even if ENABLE_SCRUBLET remains true in a shared config.
    scrublet_enabled = params.ENABLE_SCRUBLET && !params.ENABLE_DATA_HASHING
    scrublet_expected_doublet_rate = params.SCRUBLET_expected_doublet_rate != null ? (params.SCRUBLET_expected_doublet_rate as Double) : (params.SCRUBLET_assay_type == 'cc-perturb-seq' ? 0.025 : 0.08)
    Preprocessed_AnnData = PreprocessAnnData(
        trans_out_dir,
        parsed_covariate_value,
        params.QC_min_counts_per_cell,
        params.REFERENCE_transcriptome,
        params.QC_barcode_filter,
        params.QC_MAD_total_counts,
        params.QC_MAD_n_genes,
        scrublet_enabled,
        scrublet_expected_doublet_rate,
        params.SCRUBLET_n_prin_comps,
        params.SCRUBLET_adaptive_pca_fallback,
        params.SCRUBLET_failure_policy
    )

    filtered_measurement_sets = Preprocessed_AnnData.filtered_measurement_set
        .collect()
        .map { files -> files.sort { a, b -> a.getName() <=> b.getName() } }
    measurement_set_qc_dirs = Preprocessed_AnnData.measurement_set_qc
        .collect()
        .map { dirs -> dirs.sort { a, b -> a.getName() <=> b.getName() } }

    Concatenated_QC = concat_preprocessed_rna(
        filtered_measurement_sets,
        params.QC_pct_mito,
        params.QC_min_cells_per_gene
    )
    Collected_QC = collect_measurement_set_qc(
        measurement_set_qc_dirs,
        Concatenated_QC.post_concat_qc
    )

    if (file(params.REFERENCE_gtf_local_path).exists()) {
        GTF_Reference = skipGTFDownload(file(params.REFERENCE_gtf_local_path))
    }
    else {
        GTF_Reference = downloadGTF(params.REFERENCE_gtf_download_path)
    }

    emit:
    filtered_anndata_rna = Concatenated_QC.filtered_anndata_rna
    figures_dir = Collected_QC.figures_dir
    adata_rna = concat_anndata_rna
    gencode_gtf = GTF_Reference.gencode_gtf
}
