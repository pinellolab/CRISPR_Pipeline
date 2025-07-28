nextflow.enable.dsl=2

include { PreprocessAnnData } from '../../../modules/local/PreprocessAnnData'
include { CreateMuData } from '../../../modules/local/CreateMuData'
include { doublets_scrub } from '../../../modules/local/doublets_scrub'
include { prepare_assignment } from '../../../modules/local/prepare_assignment'
include { mudata_concat } from '../../../modules/local/mudata_concat'
include { guide_assignment_cleanser } from '../../../modules/local/guide_assignment_cleanser'
include { guide_assignment_sceptre } from '../../../modules/local/guide_assignment_sceptre'
include { skipGTFDownload } from '../../../modules/local/skipGTFDownload'
include { downloadGTF } from '../../../modules/local/downloadGTF'
include { prepare_guide_inference } from '../../../modules/local/prepare_guide_inference'
include { prepare_user_guide_inference } from '../../../modules/local/prepare_user_guide_inference'
include { inference_sceptre } from '../../../modules/local/inference_sceptre'
include { inference_perturbo } from '../../../modules/local/inference_perturbo'
include { inference_perturbo_trans} from '../../../modules/local/inference_perturbo_trans'
include { inference_mudata } from '../../../modules/local/inference_mudata'
include { mergedResults } from '../../../modules/local/mergedResults'
include { publishFiles } from '../../../modules/local/publishFiles'
include { mergeMudata } from '../../../modules/local/mergeMudata'

workflow process_mudata_pipeline {

    take:
    concat_anndata_rna
    trans_out_dir
    concat_anndata_guide
    guide_out_dir
    //covariate_string
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

    Prepare_assignment = prepare_assignment{MuData_Doublets.mudata_doublet}

    if (params.GUIDE_ASSIGNMENT_method == "cleanser") {
        Guide_Assignment = guide_assignment_cleanser(Prepare_assignment.prepare_assignment_mudata.flatten(), params.GUIDE_ASSIGNMENT_cleanser_probability_threshold)
        guide_assignment_collected =  Guide_Assignment.guide_assignment_mudata_output.collect()
        Mudata_concat = mudata_concat(guide_assignment_collected)
        }

    else if (params.GUIDE_ASSIGNMENT_method == "sceptre") {
        Guide_Assignment = guide_assignment_sceptre(Prepare_assignment.prepare_assignment_mudata.flatten(), params.GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold, params.GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep)
        guide_assignment_collected =  Guide_Assignment.guide_assignment_mudata_output.collect()
        Mudata_concat = mudata_concat(guide_assignment_collected)
        }

    if (params.INFERENCE_target_guide_pairing_strategy == 'predefined_pairs') {
        PrepareInference = prepare_user_guide_inference(
            Mudata_concat.concat_mudata,
            file(params.INFERENCE_predefined_pairs_to_test)
        )}
    else if (params.INFERENCE_target_guide_pairing_strategy == 'by_distance') {
        PrepareInference = prepare_guide_inference(
            Mudata_concat.concat_mudata,
            GTF_Reference.gencode_gtf,
            params.INFERENCE_max_target_distance_bp
        )}
    else if (params.INFERENCE_target_guide_pairing_strategy == 'default') {
        PrepareInference_cis = prepare_guide_inference(
            Mudata_concat.concat_mudata,
            GTF_Reference.gencode_gtf,
            params.INFERENCE_max_target_distance_bp
        )
    }

    if (params.INFERENCE_method == "sceptre"){
        TestResults = inference_sceptre(PrepareInference.mudata_inference_input)
        GuideInference = inference_mudata(TestResults.test_results, PrepareInference.mudata_inference_input, params.INFERENCE_method)
    }
    else if (params.INFERENCE_method == "perturbo"){
        GuideInference = inference_perturbo(PrepareInference.mudata_inference_input, params.INFERENCE_method, params.Multiplicity_of_infection)
    }
    else if (params.INFERENCE_method == "sceptre,perturbo") {
        SceptreResults = inference_sceptre(PrepareInference.mudata_inference_input)
        PerturboResults = inference_perturbo(PrepareInference.mudata_inference_input,  "perturbo", params.Multiplicity_of_infection)
        GuideInference = mergedResults(SceptreResults.test_results, PerturboResults.inference_mudata)
    }
    else if (params.INFERENCE_method == "default"){
        if (params.INFERENCE_target_guide_pairing_strategy != 'default') {
            error "INFERENCE_method='default' requires INFERENCE_target_guide_pairing_strategy='default'"
        }
        // Process cis results
        SceptreResults_cis = inference_sceptre(PrepareInference_cis.mudata_inference_input)
        PerturboResults_cis = inference_perturbo(PrepareInference_cis.mudata_inference_input, "perturbo", params.Multiplicity_of_infection)
        GuideInference_cis = mergedResults(SceptreResults_cis.test_results, PerturboResults_cis.inference_mudata)
        // Process trans results
        GuideInference_trans = inference_perturbo_trans(Mudata_concat.concat_mudata, "perturbo", params.Multiplicity_of_infection)

        // Rename tsv outputs to avoid conflicts
        cis_per_element = GuideInference_cis.per_element_output.map { file -> file.copyTo(file.parent.resolve("cis-${file.name}")) }
        cis_per_guide = GuideInference_cis.per_guide_output.map { file -> file.copyTo(file.parent.resolve("cis-${file.name}")) }

        trans_per_element = GuideInference_trans.per_element_output.map { file -> file.copyTo(file.parent.resolve("trans-${file.name}")) }
        trans_per_guide = GuideInference_trans.per_guide_output.map { file -> file.copyTo(file.parent.resolve("trans-${file.name}")) }

        PublishFiles = publishFiles(cis_per_element, cis_per_guide, trans_per_element, trans_per_guide)

        // Rename h5mu outputs to avoid conflicts
        cis_file = GuideInference_cis.inference_mudata.map { file ->
            file.copyTo(file.parent.resolve('cis_inference_mudata.h5mu'))
        }
        trans_file = GuideInference_trans.inference_mudata.map { file ->
            file.copyTo(file.parent.resolve('trans_inference_mudata.h5mu'))
        }

        GuideInference = mergeMudata(cis_file, trans_file)

    }


    emit:
    inference_mudata = GuideInference.inference_mudata
    gencode_gtf = GTF_Reference.gencode_gtf
    figures_dir = Preprocessed_AnnData.figures_dir
    adata_rna = Preprocessed_AnnData.adata_rna
    filtered_anndata_rna = Preprocessed_AnnData.filtered_anndata_rna
    adata_guide = MuData.adata_guide

}
