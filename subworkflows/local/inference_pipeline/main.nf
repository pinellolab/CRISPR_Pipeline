nextflow.enable.dsl=2

include { prepare_guide_inference } from '../../../modules/local/prepare_guide_inference'
include { prepare_all_guide_inference } from '../../../modules/local/prepare_all_guide_inference'
include { prepare_user_guide_inference } from '../../../modules/local/prepare_user_guide_inference'
include { inference_sceptre } from '../../../modules/local/inference_sceptre'
include { inference_perturbo } from '../../../modules/local/inference_perturbo'
include { inference_perturbo_trans } from '../../../modules/local/inference_perturbo_trans'
include { mergedResults } from '../../../modules/local/mergedResults'
include { publishFiles } from '../../../modules/local/publishFiles'
include { mergeMudata } from '../../../modules/local/mergeMudata'

workflow inference_pipeline {

    take:
    mudata_concat
    gtf_reference

    main:

    if (params.INFERENCE_target_guide_pairing_strategy == 'predefined_pairs') {
        PrepareInference = prepare_user_guide_inference(
            mudata_concat,
            file(params.INFERENCE_predefined_pairs_to_test)
        )}
    else if (params.INFERENCE_target_guide_pairing_strategy == 'by_distance') {
        PrepareInference = prepare_guide_inference(
            mudata_concat,
            gtf_reference,
            params.INFERENCE_max_target_distance_bp,
            false
        )}
    else if (params.INFERENCE_target_guide_pairing_strategy == 'default') {
        PrepareInference = prepare_guide_inference(
            mudata_concat,
            gtf_reference,
            params.INFERENCE_max_target_distance_bp,
            true
        )
    }

    // Determine the mudata input once (avoid duplicate variable definitions)
    def mudata_input

    if (params.INFERENCE_target_guide_pairing_strategy != 'all_by_all') {
        mudata_input = PrepareInference.mudata_inference_input
    } else {
        mudata_input = mudata_concat
    }

    if (params.INFERENCE_method == "sceptre"){
        TestResults = inference_sceptre(mudata_input)
        GuideInference = TestResults.inference_mudata
    }
    else if (params.INFERENCE_method == "perturbo"){
        TestResults = inference_perturbo(mudata_input, params.INFERENCE_method, params.Multiplicity_of_infection)
        GuideInference = TestResults.inference_mudata
    }
    else if (params.INFERENCE_method == "sceptre,perturbo") {
        SceptreResults = inference_sceptre(mudata_input)
        PerturboResults = inference_perturbo(mudata_input,  "perturbo", params.Multiplicity_of_infection)
        GuideInference = mergedResults(
            SceptreResults.per_guide_output,
            SceptreResults.per_element_output,
            PerturboResults.per_guide_output,
            PerturboResults.per_element_output,
            mudata_input
        )
    }
    else if (params.INFERENCE_method == "default"){
        if (params.INFERENCE_target_guide_pairing_strategy != 'default') {
            error "INFERENCE_method='default' requires INFERENCE_target_guide_pairing_strategy='default'"
        }
        // Process cis results
        SceptreResults_cis = inference_sceptre(PrepareInference.mudata_inference_input)
        PerturboResults_cis = inference_perturbo(PrepareInference.mudata_inference_input, "perturbo", params.Multiplicity_of_infection)
        MergedInference_cis = mergedResults(
            SceptreResults_cis.per_guide_output,
            SceptreResults_cis.per_element_output,
            PerturboResults_cis.per_guide_output,
            PerturboResults_cis.per_element_output,
            PrepareInference.mudata_inference_input
        )
        // Process trans results - use concat_mudata directly
        MergedInference_trans = inference_perturbo_trans(mudata_concat, "perturbo", params.Multiplicity_of_infection, PerturboResults_cis.inference_mudata)

        // Rename tsv outputs to avoid conflicts
        cis_per_element = MergedInference_cis.per_element_output.map { file -> file.copyTo(file.parent.resolve("cis-${file.name}")) }
        cis_per_guide = MergedInference_cis.per_guide_output.map { file -> file.copyTo(file.parent.resolve("cis-${file.name}")) }

        trans_per_element = MergedInference_trans.per_element_output.map { file -> file.copyTo(file.parent.resolve("trans-${file.name}")) }
        trans_per_guide = MergedInference_trans.per_guide_output.map { file -> file.copyTo(file.parent.resolve("trans-${file.name}")) }

        PublishFiles = publishFiles(cis_per_element, cis_per_guide, trans_per_element, trans_per_guide)

        // Rename h5mu outputs to avoid conflicts
        cis_file = MergedInference_cis.inference_mudata.map { file ->
            file.copyTo(file.parent.resolve('cis_inference_mudata.h5mu'))
        }
        trans_file = MergedInference_trans.inference_mudata.map { file ->
            file.copyTo(file.parent.resolve('trans_inference_mudata.h5mu'))
        }

        MergedInference = mergeMudata(
            MergedInference_cis.per_guide_output,
            MergedInference_cis.per_element_output,
            MergedInference_trans.per_guide_output,
            MergedInference_trans.per_element_output,
            PrepareInference.mudata_concat,
        )

    }

    emit:
    inference_mudata = MergedInference.inference_mudata

}
