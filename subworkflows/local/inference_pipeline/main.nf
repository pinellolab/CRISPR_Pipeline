nextflow.enable.dsl=2

include { prepare_guide_inference } from '../../../modules/local/prepare_guide_inference'
include { prepare_all_guide_inference } from '../../../modules/local/prepare_all_guide_inference'
include { prepare_user_guide_inference } from '../../../modules/local/prepare_user_guide_inference'
include { inference_sceptre } from '../../../modules/local/inference_sceptre'
include { inference_perturbo } from '../../../modules/local/inference_perturbo'
include { inference_perturbo as inference_perturbo_trans } from '../../../modules/local/inference_perturbo'
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
    else if (params.INFERENCE_target_guide_pairing_strategy == 'all_by_all') {
        // Skip prepare_all_guide_inference and use mudata directly
        PrepareInference = Channel.empty()
    }
    else if (params.INFERENCE_target_guide_pairing_strategy == 'default') {
        PrepareInference_cis = prepare_guide_inference(
            mudata_concat,
            gtf_reference,
            params.INFERENCE_max_target_distance_bp,
            true
        )
        // Skip prepare_all_guide_inference for trans analysis
        PrepareInference_trans = Channel.empty()
    }

    if (params.INFERENCE_method == "sceptre"){
        def mudata_input = params.INFERENCE_target_guide_pairing_strategy == 'all_by_all' ? mudata_concat : PrepareInference.mudata_inference_input
        TestResults = inference_sceptre(mudata_input)
        GuideInference = TestResults.inference_mudata
    }
    else if (params.INFERENCE_method == "perturbo"){
        def mudata_input = params.INFERENCE_target_guide_pairing_strategy == 'all_by_all' ? mudata_concat : PrepareInference.mudata_inference_input
        TestResults = inference_perturbo(mudata_input, params.INFERENCE_method, params.Multiplicity_of_infection)
        GuideInference = TestResults.inference_mudata
    }
    else if (params.INFERENCE_method == "sceptre,perturbo") {
        def mudata_input = params.INFERENCE_target_guide_pairing_strategy == 'all_by_all' ? mudata_concat : PrepareInference.mudata_inference_input
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
        SceptreResults_cis = inference_sceptre(PrepareInference_cis.mudata_inference_input)
        PerturboResults_cis = inference_perturbo(PrepareInference_cis.mudata_inference_input, "perturbo", params.Multiplicity_of_infection)
        GuideInference_cis = mergedResults(
            SceptreResults_cis.per_guide_output,
            SceptreResults_cis.per_element_output,
            PerturboResults_cis.per_guide_output,
            PerturboResults_cis.per_element_output,
            PrepareInference_cis.mudata_inference_input
        )
        // Process trans results - use concat_mudata directly
        GuideInference_trans = inference_perturbo_trans(mudata_concat, "perturbo", params.Multiplicity_of_infection)

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

        GuideInference = mergeMudata(
            GuideInference_cis.per_guide_output,
            GuideInference_cis.per_element_output,
            GuideInference_trans.per_guide_output,
            GuideInference_trans.per_element_output,
            PrepareInference_cis.mudata_inference_input
        )

    }

    emit:
    inference_mudata = GuideInference.inference_mudata

}
