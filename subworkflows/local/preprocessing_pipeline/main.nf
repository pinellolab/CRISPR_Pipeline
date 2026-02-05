include { PreprocessAnnData } from '../../../modules/local/PreprocessAnnData'
include { skipGTFDownload } from '../../../modules/local/skipGTFDownload'
include { downloadGTF } from '../../../modules/local/downloadGTF'

workflow preprocessing_pipeline {

    take:
    concat_anndata_rna
    trans_out_dir

    main:
    Preprocessed_AnnData = PreprocessAnnData(
        concat_anndata_rna,
        trans_out_dir.flatten().first(),
        params.QC_min_genes_per_cell,
        params.QC_min_cells_per_gene,
        params.QC_pct_mito,
        params.REFERENCE_transcriptome,
        params.QC_barcode_filter
    )

    if (file(params.REFERENCE_gtf_local_path).exists()) {
        GTF_Reference = skipGTFDownload(file(params.REFERENCE_gtf_local_path))
    }
    else {
        GTF_Reference = downloadGTF(params.REFERENCE_gtf_download_path)
    }

    emit:
    filtered_anndata_rna = Preprocessed_AnnData.filtered_anndata_rna
    figures_dir = Preprocessed_AnnData.figures_dir
    adata_rna = Preprocessed_AnnData.adata_rna
    gencode_gtf = GTF_Reference.gencode_gtf
}
