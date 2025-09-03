
nextflow.enable.dsl=2

include { PreprocessAnnData } from '../../../modules/local/PreprocessAnnData'
include { CreateMuData_HASHING } from '../../../modules/local/CreateMuData_HASHING'
include { demultiplex } from '../../../modules/local/demultiplex'
include { filter_hashing } from '../../../modules/local/filter_hashing'
include { hashing_concat } from '../../../modules/local/hashing_concat'
include { guide_assignment_pipeline } from '../guide_assignment_pipeline'
include { skipGTFDownload } from '../../../modules/local/skipGTFDownload'
include { downloadGTF } from '../../../modules/local/downloadGTF'
include { inference_pipeline } from '../inference_pipeline'

workflow process_mudata_pipeline_HASHING {

    take:
    concat_anndata_rna
    trans_out_dir
    concat_anndata_guide
    guide_out_dir
    concat_anndata_hashing
    hashing_out_dir
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

    Hashing_Filtered = filter_hashing(
        Preprocessed_AnnData.filtered_anndata_rna,
        concat_anndata_hashing
        )

    Demultiplex = demultiplex(Hashing_Filtered.hashing_filtered_anndata.flatten())

    hashing_demux_anndata_collected =Demultiplex.hashing_demux_anndata.collect()
    hashing_demux_anndata_collected.view()

    hashing_demux_unfiltered_anndata_collected =Demultiplex.hashing_demux_unfiltered_anndata.collect()
    hashing_demux_unfiltered_anndata_collected.view()

    Hashing_Concat = hashing_concat(hashing_demux_anndata_collected, hashing_demux_unfiltered_anndata_collected)

    if (file(params.REFERENCE_gtf_local_path).exists()) {
        GTF_Reference = skipGTFDownload(file(params.REFERENCE_gtf_local_path))
    }
    else {
        GTF_Reference = downloadGTF(params.REFERENCE_gtf_download_path)
    }

    MuData = CreateMuData_HASHING(
        Preprocessed_AnnData.filtered_anndata_rna,
        concat_anndata_guide,
        Hashing_Concat.concatenated_hashing_demux,
        ch_guide_design,
        GTF_Reference.gencode_gtf,
        params.Multiplicity_of_infection,
        params.GUIDE_ASSIGNMENT_capture_method
        )

    Mudata_concat = guide_assignment_pipeline(MuData.mudata)

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
    adata_hashing = Hashing_Filtered.adata_hashing
    adata_demux = Hashing_Concat.concatenated_hashing_demux
    adata_unfiltered_demux = Hashing_Concat.concatenated_hashing_unfiltered_demux
}
