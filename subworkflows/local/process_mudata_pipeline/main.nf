nextflow.enable.dsl=2

include { PreprocessAnnData } from '../../../modules/local/PreprocessAnnData'
include { CreateMuData } from '../../../modules/local/CreateMuData'
include { doublets_scrub } from '../../../modules/local/doublets_scrub'
include { guide_assignment_pipeline } from '../guide_assignment_pipeline'
include { skipGTFDownload } from '../../../modules/local/skipGTFDownload'
include { downloadGTF } from '../../../modules/local/downloadGTF'
include { inference_pipeline } from '../inference_pipeline'

workflow process_mudata_pipeline {

    take:
    concat_anndata_rna
    trans_out_dir
    concat_anndata_guide
    guide_out_dir
    covariate_string
    ch_guide_design

    main:

    Preprocessed_AnnData = PreprocessAnnData(
        concat_anndata_rna,
        trans_out_dir.flatten().first(),
        params.QC_min_genes_per_cell,
        params.QC_min_cells_per_gene,
        params.QC_pct_mito,
        params.REFERENCE_transcriptome
        )

    if (file(params.REFERENCE_gtf_local_path).exists()) {
        GTF_Reference = skipGTFDownload(file(params.REFERENCE_gtf_local_path))
    }
    else {
        GTF_Reference = downloadGTF(params.REFERENCE_gtf_download_path)
    }

    MuData = CreateMuData(
        Preprocessed_AnnData.filtered_anndata_rna,
        concat_anndata_guide,
        ch_guide_design,
        GTF_Reference.gencode_gtf,
        params.Multiplicity_of_infection,
        params.GUIDE_ASSIGNMENT_capture_method
        )

    MuData_Doublets = doublets_scrub(MuData.mudata)

    Mudata_concat = guide_assignment_pipeline(MuData_Doublets.mudata_doublet)

    GuideInference = inference_pipeline(
        Mudata_concat.concat_mudata,
        GTF_Reference.gencode_gtf
    )


    emit:
    inference_mudata = GuideInference.inference_mudata
    gencode_gtf = GTF_Reference.gencode_gtf
    figures_dir = Preprocessed_AnnData.figures_dir
    adata_rna = Preprocessed_AnnData.adata_rna
    filtered_anndata_rna = Preprocessed_AnnData.filtered_anndata_rna
    adata_guide = MuData.adata_guide

}
