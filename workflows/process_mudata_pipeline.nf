nextflow.enable.dsl=2

include { PreprocessAnnData } from '../processes/PreprocessAnnData.nf'
include { CreateMuData } from '../processes/CreateMuData.nf' 
include { doublets_scrub } from '../processes/doublets_scrub.nf'
include { prepare_assignment } from '../processes/prepare_assignment.nf'
include { mudata_concat } from '../processes/mudata_concat.nf'
include { guide_assignment_cleanser } from '../processes/guide_assignment_cleanser.nf'
include { guide_assignment_sceptre } from '../processes/guide_assignment_sceptre.nf'
include { skipGTFDownload } from '../processes/skipGTFDownload.nf'
include { downloadGTF } from '../processes/downloadGTF.nf'
include { prepare_guide_inference } from '../processes/prepare_guide_inference.nf'
include { prepare_all_guide_inference } from '../processes/prepare_all_guide_inference.nf'
include { prepare_user_guide_inference } from '../processes/prepare_user_guide_inference.nf'
include { inference_sceptre } from '../processes/inference_sceptre.nf'
include { inference_perturbo } from '../processes/inference_perturbo.nf'
include { inference_mudata } from '../processes/inference_mudata.nf'
include { mergedResults } from '../processes/mergedResults.nf'

workflow process_mudata_pipeline {

    take:
    concat_anndata_rna
    trans_out_dir
    concat_anndata_guide
    guide_out_dir
    covariate_string

    main:

    Preprocessed_AnnData = PreprocessAnnData(
        concat_anndata_rna,
        trans_out_dir.flatten().first(),
        params.min_genes,
        params.min_cells,
        params.pct_mito,
        params.transcriptome
        )
    
    if (file(params.gtf_local_path).exists()) {
        GTF_Reference = skipGTFDownload(file(params.gtf_local_path))
    }
    else {
        GTF_Reference = downloadGTF(params.gtf_download_path)
    }

    MuData = CreateMuData(
        Preprocessed_AnnData.filtered_anndata_rna,
        concat_anndata_guide, 
        file(params.guide_metadata),
        GTF_Reference.gencode_gtf,
        params.moi,
        params.capture_method
        )

    MuData_Doublets = doublets_scrub(MuData.mudata) 

    Prepare_assignment = prepare_assignment{MuData_Doublets.mudata_doublet}
    prepare_assignment_collected = Prepare_assignment.prepare_assignment_mudata.collect()

    if (params.assignment_method == "cleanser") {
        Guide_Assignment = guide_assignment_cleanser(prepare_assignment_collected, params.THRESHOLD)
        guide_assignment_collected =  Guide_Assignment.guide_assignment_mudata_output.collect()
        Mudata_concat = mudata_concat(guide_assignment_collected)
        }

    else if (params.assignment_method == "sceptre") {
        Guide_Assignment = guide_assignment_sceptre(prepare_assignment_collected)
        guide_assignment_collected =  Guide_Assignment.guide_assignment_mudata_output.collect()
        Mudata_concat = mudata_concat(guide_assignment_collected)
        }

    if (params.inference_option == 'predefined_pairs') {
        PrepareInference = prepare_user_guide_inference(
            Mudata_concat.concat_mudata,
            file(params.user_inference)
        )}
    else if (params.inference_option == 'by_distance') {
        PrepareInference = prepare_guide_inference(
            Mudata_concat.concat_mudata,
            GTF_Reference.gencode_gtf,
            params.distance_from_center
        )}
    else if (params.inference_option == 'all_by_all') {
        PrepareInference = prepare_all_guide_inference(
            Mudata_concat.concat_mudata,
            GTF_Reference.gencode_gtf
        )}

    if (params.inference_method == "sceptre"){
        TestResults = inference_sceptre(PrepareInference.mudata_inference_input, covariate_string)
        GuideInference = inference_mudata(TestResults.test_results, PrepareInference.mudata_inference_input, params.inference_method)
    }
    else if (params.inference_method == "perturbo"){
        GuideInference = inference_perturbo(PrepareInference.mudata_inference_input, params.inference_method)
    }
    else if (params.inference_method == "sceptre,perturbo") {
        SceptreResults = inference_sceptre(PrepareInference.mudata_inference_input, covariate_string)
        PerturboResults = inference_perturbo(PrepareInference.mudata_inference_input,  "perturbo")
        GuideInference = mergedResults(SceptreResults.test_results, PerturboResults.inference_mudata)
    }


    emit:
    inference_mudata = GuideInference.inference_mudata
    gencode_gtf = GTF_Reference.gencode_gtf
    figures_dir = Preprocessed_AnnData.figures_dir
    adata_rna = Preprocessed_AnnData.adata_rna
    filtered_anndata_rna = Preprocessed_AnnData.filtered_anndata_rna
    adata_guide = MuData.adata_guide

}
