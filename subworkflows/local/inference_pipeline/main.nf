nextflow.enable.dsl=2

include { prepare_guide_inference } from '../../../modules/local/prepare_guide_inference'
include { prepare_all_guide_inference } from '../../../modules/local/prepare_all_guide_inference'
include { prepare_user_guide_inference } from '../../../modules/local/prepare_user_guide_inference'
include { inference_sceptre } from '../../../modules/local/inference_sceptre'
include { sceptre_chunk_prepare } from '../../../modules/local/sceptre_chunk_prepare'
include { sceptre_chunk_merge } from '../../../modules/local/sceptre_chunk_merge'
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
    sort_paths = { paths -> paths.sort { a, b -> a.toString() <=> b.toString() } }

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
    } else{
        error("Invalid INFERENCE_target_guide_pairing_strategy: ${params.INFERENCE_target_guide_pairing_strategy}")
    }

    // Determine the mudata input once (avoid duplicate variable definitions)
    def mudata_input

    if (params.INFERENCE_target_guide_pairing_strategy != 'all_by_all') {
        mudata_input = PrepareInference.mudata_inference_input
    } else {
        mudata_input = mudata_concat
    }

    if (params.INFERENCE_method == "sceptre"){
        SceptreChunkInput = sceptre_chunk_prepare(mudata_input)
        SceptreChunkResults = inference_sceptre(SceptreChunkInput.mudata_chunks.flatten())
        TestResults = sceptre_chunk_merge(
            SceptreChunkResults.per_guide_output.collect().map(sort_paths),
            SceptreChunkResults.per_element_output.collect().map(sort_paths),
            mudata_input,
            SceptreChunkInput.chunk_manifest
        )
        FinalInference = TestResults.inference_mudata
    }
    else if (params.INFERENCE_method == "perturbo"){
        TestResults = inference_perturbo(mudata_input, params.INFERENCE_method, params.Multiplicity_of_infection)
        FinalInference = TestResults.inference_mudata
    }
    else if (params.INFERENCE_method == "sceptre,perturbo") {
        SceptreChunkInput = sceptre_chunk_prepare(mudata_input)
        SceptreChunkResults = inference_sceptre(SceptreChunkInput.mudata_chunks.flatten())
        SceptreResults = sceptre_chunk_merge(
            SceptreChunkResults.per_guide_output.collect().map(sort_paths),
            SceptreChunkResults.per_element_output.collect().map(sort_paths),
            mudata_input,
            SceptreChunkInput.chunk_manifest
        )
        PerturboResults = inference_perturbo(mudata_input,  "perturbo", params.Multiplicity_of_infection)
        MergedInference = mergedResults(
            SceptreResults.per_guide_output,
            SceptreResults.per_element_output,
            PerturboResults.per_guide_output,
            PerturboResults.per_element_output,
            mudata_input
        )
        FinalInference = MergedInference.inference_mudata
    }
    else if (params.INFERENCE_method == "default"){
        if (params.INFERENCE_target_guide_pairing_strategy != 'default') {
            error "INFERENCE_method='default' requires INFERENCE_target_guide_pairing_strategy='default'"
        }
        // Process cis results
        SceptreChunkInput_cis = sceptre_chunk_prepare(PrepareInference.mudata_inference_input)
        SceptreChunkResults_cis = inference_sceptre(SceptreChunkInput_cis.mudata_chunks.flatten())
        SceptreResults_cis = sceptre_chunk_merge(
            SceptreChunkResults_cis.per_guide_output.collect().map(sort_paths),
            SceptreChunkResults_cis.per_element_output.collect().map(sort_paths),
            PrepareInference.mudata_inference_input,
            SceptreChunkInput_cis.chunk_manifest
        )
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

        MergedInference = mergeMudata(
            MergedInference_cis.per_guide_output,
            MergedInference_cis.per_element_output,
            MergedInference_trans.per_guide_output,
            MergedInference_trans.per_element_output,
            mudata_concat,
        )
        FinalInference = MergedInference.inference_mudata
    } else {
        error("Invalid INFERENCE_method: ${params.INFERENCE_method}. Valid options: sceptre, perturbo, sceptre,perturbo, default")
    }

    emit:
    inference_mudata = FinalInference

}
