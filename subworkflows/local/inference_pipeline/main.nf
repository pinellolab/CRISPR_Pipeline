nextflow.enable.dsl=2

include { prepare_guide_inference } from '../../../modules/local/prepare_guide_inference'
include { prepare_all_guide_inference } from '../../../modules/local/prepare_all_guide_inference'
include { prepare_user_guide_inference } from '../../../modules/local/prepare_user_guide_inference'
include { inference_sceptre } from '../../../modules/local/inference_sceptre'
include { sceptre_chunk_prepare } from '../../../modules/local/sceptre_chunk_prepare'
include { sceptre_chunk_merge } from '../../../modules/local/sceptre_chunk_merge'
include { inference_perturbo } from '../../../modules/local/inference_perturbo'
include { mergedResults } from '../../../modules/local/mergedResults'
include { publishFiles } from '../../../modules/local/publishFiles'
include { mergeMudata } from '../../../modules/local/mergeMudata'
include { buildCatalogElement } from '../../../modules/local/buildCatalogElement'
include { buildCatalogGuide } from '../../../modules/local/buildCatalogGuide'

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
        TestResults = inference_perturbo(mudata_input, mudata_input, params.INFERENCE_method)
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
        PerturboResults = inference_perturbo(mudata_input, mudata_input, "perturbo")
        MergedInference = mergedResults(
            SceptreResults.per_guide_output,
            SceptreResults.per_element_output,
            PerturboResults.local_per_guide_output,
            PerturboResults.local_per_element_output,
            mudata_input,
            true  // this workflow ends here, so it needs the MuData
        )
        FinalInference = MergedInference.inference_mudata
    }
    else if (params.INFERENCE_method == "default"){
        if (params.INFERENCE_target_guide_pairing_strategy != 'default') {
            error "INFERENCE_method='default' requires INFERENCE_target_guide_pairing_strategy='default'"
        }
        // Process local-analysis results
        SceptreChunkInput_local = sceptre_chunk_prepare(PrepareInference.mudata_inference_input)
        SceptreChunkResults_local = inference_sceptre(SceptreChunkInput_local.mudata_chunks.flatten())
        SceptreResults_local = sceptre_chunk_merge(
            SceptreChunkResults_local.per_guide_output.collect().map(sort_paths),
            SceptreChunkResults_local.per_element_output.collect().map(sort_paths),
            PrepareInference.mudata_inference_input,
            SceptreChunkInput_local.chunk_manifest
        )
        // One PerTurbo run on every gene of concat_mudata yields both the local
        // (requested-pair) and the global (transcriptome-wide) tables; the pairs
        // come from the prepared inference input.
        PerturboResults = inference_perturbo(mudata_concat, PrepareInference.mudata_inference_input, "perturbo")
        MergedInference_local = mergedResults(
            SceptreResults_local.per_guide_output,
            SceptreResults_local.per_element_output,
            PerturboResults.local_per_guide_output,
            PerturboResults.local_per_element_output,
            PrepareInference.mudata_inference_input,
            false  // mergeMudata assembles the published MuData from these tables
        )

        MergedInference = mergeMudata(
            MergedInference_local.per_guide_output,
            MergedInference_local.per_element_output,
            PerturboResults.global_per_guide_output,
            PerturboResults.global_per_element_output,
            mudata_concat,
        )
        // Catalog construction is intentionally independent from mergeMudata:
        // each large table is cacheable and resumable on its own.
        buildCatalogElement(
            MergedInference.local_analysis_per_element_output,
            MergedInference.global_analysis_per_element_output,
            MergedInference.inference_mudata,
        )
        buildCatalogGuide(
            MergedInference.local_analysis_per_guide_output,
            MergedInference.global_analysis_per_guide_output,
            MergedInference.inference_mudata,
        )
        FinalInference = MergedInference.inference_mudata
    } else {
        error("Invalid INFERENCE_method: ${params.INFERENCE_method}. Valid options: sceptre, perturbo, sceptre,perturbo, default")
    }

    emit:
    inference_mudata = FinalInference

}
