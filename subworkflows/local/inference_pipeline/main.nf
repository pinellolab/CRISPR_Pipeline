nextflow.enable.dsl=2

include { prepare_guide_inference } from '../../../modules/local/prepare_guide_inference'
include { prepare_guide_inference as prepare_global_guide_inference } from '../../../modules/local/prepare_guide_inference'
include { prepare_all_guide_inference } from '../../../modules/local/prepare_all_guide_inference'
include { prepare_user_guide_inference } from '../../../modules/local/prepare_user_guide_inference'
include { inference_sceptre } from '../../../modules/local/inference_sceptre'
include { sceptre_chunk_prepare } from '../../../modules/local/sceptre_chunk_prepare'
include { sceptre_chunk_merge } from '../../../modules/local/sceptre_chunk_merge'
include { inference_perturbo } from '../../../modules/local/inference_perturbo'
include { mergedResults } from '../../../modules/local/mergedResults'
include { publishFiles } from '../../../modules/local/publishFiles'
include { mergeMudata } from '../../../modules/local/mergeMudata'
include { export_igvf_inference_tsv } from '../../../modules/local/export_igvf_inference_tsv'
include { buildCatalogElement } from '../../../modules/local/buildCatalogElement'
include { buildCatalogGuide } from '../../../modules/local/buildCatalogGuide'
include { resolveControlGroupFromParams } from '../../../modules/local/control_group'
include { resolve_inference_mode } from '../../../modules/local/resolve_inference_mode'
include { ntc_calibration_pipeline } from '../ntc_calibration_pipeline'

workflow inference_pipeline {

    take:
    mudata_concat
    gtf_reference

    main:
    sort_paths = { paths -> paths.sort { a, b -> a.toString() <=> b.toString() } }

    // The cells a perturbation is compared against, decided once for both methods
    // and logged once. This is the line to read when you wonder why a screen used
    // one pool. It aborts here if the setting is unusable, rather than three hours
    // into a run.
    def control_group = resolveControlGroupFromParams()
    log.info(control_group.log_line)

    // The parameters alone leave 'auto' half-decided: whether the screen is
    // low-MOI and has a control pool is a property of the assignments. Each
    // method branch below measures the full guide matrix once
    // (modules/local/resolve_inference_mode) and both processes read the same
    // map: the control group and MOI mode for SCEPTRE, the CRT pool for
    // PerTurbo, never 'from-moi', so neither driver settles the design alone.
    def modeFromJson = { f ->
        def j = new groovy.json.JsonSlurper().parseText(f.text)
        [
            setting: j.setting,
            moi: j.analysis_moi,
            analysis_moi: j.analysis_moi,
            declared_moi: j.declared_moi,
            sceptre_control_group: j.sceptre_control_group,
            perturbo_crt_pool: j.perturbo_crt_pool,
            perturbo_pool_effective: j.perturbo_pool_effective,
            sceptre_provenance: j.sceptre_provenance,
            perturbo_provenance: j.perturbo_provenance,
            overrides: j.overrides,
            reason: j.reason,
            log_line: j.log_line,
        ]
    }
    def resolution_setting = control_group.setting
    def resolution_moi = params.Multiplicity_of_infection ?: ''
    def resolution_perturbo_pool = params.INFERENCE_PERTURBO_CRT_POOL ?: ''
    def resolution_sceptre_group = params.INFERENCE_SCEPTRE_control_group ?: ''

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

    def calibration_placeholder = Channel.value(file("${workflow.projectDir}/assets/calibration_empty"))
    // A command-line `--ENABLE_NTC_CALIBRATION false` reaches here as the text
    // 'false' on Nextflow 26 (a boolean on 25), and Groovy reads any non-empty
    // text as true. Decide on the text, so the flag means what it says.
    def calibration_flag = params.ENABLE_NTC_CALIBRATION
    def enable_calibration = calibration_flag == true ||
        (calibration_flag != null && calibration_flag.toString().trim().toLowerCase() in ['true', '1', 'yes'])
    def run_calibration = enable_calibration && params.INFERENCE_method in ['default', 'sceptre,perturbo']
    if (enable_calibration && !run_calibration) {
        log.info("Non-targeting calibration: skipped, it needs both methods (INFERENCE_method='${params.INFERENCE_method}').")
    }

    if (params.INFERENCE_method == "sceptre"){
        mode = resolve_inference_mode(mudata_input, resolution_setting, resolution_moi, resolution_perturbo_pool, resolution_sceptre_group)
            .json.map(modeFromJson).first()
        mode.subscribe { resolved -> log.info(resolved.log_line) }
        calibration_output = calibration_placeholder
        SceptreChunkInput = sceptre_chunk_prepare(mudata_input)
        SceptreChunkResults = inference_sceptre(SceptreChunkInput.mudata_chunks.flatten(), mode)
        TestResults = sceptre_chunk_merge(
            SceptreChunkResults.per_guide_output.collect().map(sort_paths),
            SceptreChunkResults.per_element_output.collect().map(sort_paths),
            mudata_input,
            SceptreChunkInput.chunk_manifest
        )
        FinalInference = TestResults.inference_mudata
    }
    else if (params.INFERENCE_method == "perturbo"){
        mode = resolve_inference_mode(mudata_input, resolution_setting, resolution_moi, resolution_perturbo_pool, resolution_sceptre_group)
            .json.map(modeFromJson).first()
        mode.subscribe { resolved -> log.info(resolved.log_line) }
        calibration_output = calibration_placeholder
        TestResults = inference_perturbo(mudata_input, mudata_input, params.INFERENCE_method, mode)
        FinalInference = TestResults.inference_mudata
    }
    else if (params.INFERENCE_method == "sceptre,perturbo") {
        ModeJson = resolve_inference_mode(mudata_input, resolution_setting, resolution_moi, resolution_perturbo_pool, resolution_sceptre_group)
        mode = ModeJson.json.map(modeFromJson).first()
        mode.subscribe { resolved -> log.info(resolved.log_line) }
        SceptreChunkInput = sceptre_chunk_prepare(mudata_input)
        SceptreChunkResults = inference_sceptre(SceptreChunkInput.mudata_chunks.flatten(), mode)
        SceptreResults = sceptre_chunk_merge(
            SceptreChunkResults.per_guide_output.collect().map(sort_paths),
            SceptreChunkResults.per_element_output.collect().map(sort_paths),
            mudata_input,
            SceptreChunkInput.chunk_manifest
        )
        PerturboResults = inference_perturbo(mudata_input, mudata_input, "perturbo", mode)
        MergedInference = mergedResults(
            SceptreResults.per_guide_output,
            SceptreResults.per_element_output,
            PerturboResults.local_per_guide_output,
            PerturboResults.local_per_element_output,
            mudata_input,
            true  // this workflow ends here, so it needs the MuData
        )
        FinalInference = MergedInference.inference_mudata
        if (run_calibration) {
            Calibration = ntc_calibration_pipeline(
                mudata_input, mudata_input, mode, ModeJson.json,
                SceptreResults.per_guide_output, SceptreResults.per_element_output,
                PerturboResults.global_per_guide_output, PerturboResults.global_per_element_output
            )
            calibration_output = Calibration.calibration_output
        } else {
            calibration_output = calibration_placeholder
        }
    }
    else if (params.INFERENCE_method == "default"){
        if (params.INFERENCE_target_guide_pairing_strategy != 'default') {
            error "INFERENCE_method='default' requires INFERENCE_target_guide_pairing_strategy='default'"
        }
        // PerTurbo fits the full gene universe, but it must receive the same
        // configured cell design as the local SCEPTRE input. Passing the raw
        // concatenated MuData here silently bypassed INFERENCE_COV_* and could
        // retain a constant batch that SCEPTRE correctly skipped. Prepare a
        // second, un-subset copy with the identical pair universe instead.
        PrepareInferenceGlobal = prepare_global_guide_inference(
            mudata_concat,
            gtf_reference,
            params.INFERENCE_max_target_distance_bp,
            false
        )
        // The mode is measured on the full-guide copy PerTurbo fits; SCEPTRE's
        // cis subset drops untested guides and could misread the design.
        ModeJson = resolve_inference_mode(PrepareInferenceGlobal.mudata_inference_input, resolution_setting, resolution_moi, resolution_perturbo_pool, resolution_sceptre_group)
        mode = ModeJson.json.map(modeFromJson).first()
        mode.subscribe { resolved -> log.info(resolved.log_line) }
        // Process local-analysis results
        SceptreChunkInput_local = sceptre_chunk_prepare(PrepareInference.mudata_inference_input)
        SceptreChunkResults_local = inference_sceptre(SceptreChunkInput_local.mudata_chunks.flatten(), mode)
        SceptreResults_local = sceptre_chunk_merge(
            SceptreChunkResults_local.per_guide_output.collect().map(sort_paths),
            SceptreChunkResults_local.per_element_output.collect().map(sort_paths),
            PrepareInference.mudata_inference_input,
            SceptreChunkInput_local.chunk_manifest
        )
        // One PerTurbo run on every gene of concat_mudata yields both the local
        // (requested-pair) and the global (transcriptome-wide) tables; the pairs
        // come from the prepared inference input.
        PerturboResults = inference_perturbo(PrepareInferenceGlobal.mudata_inference_input, PrepareInference.mudata_inference_input, "perturbo", mode)
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
            PrepareInferenceGlobal.mudata_inference_input,
        )
        if (params.OUTPUT_IGVF_PACKAGE && params.OUTPUT_IGVF_TSV && params.ANALYSIS_mode == 'full' && !params.INFERENCE_input_mudata) {
            export_igvf_inference_tsv(
                MergedInference.local_analysis_per_element_output,
                MergedInference.local_analysis_per_guide_output,
                MergedInference.global_analysis_per_element_output,
                MergedInference.global_analysis_per_guide_output,
            )
        }
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
        // The non-targeting calibration pass: SCEPTRE on the NT pseudo-elements,
        // joined to PerTurbo's global rows. Separate from the main run and its
        // families; see subworkflows/local/ntc_calibration_pipeline.
        if (run_calibration) {
            Calibration = ntc_calibration_pipeline(
                PrepareInferenceGlobal.mudata_inference_input, PrepareInference.mudata_inference_input, mode, ModeJson.json,
                SceptreResults_local.per_guide_output, SceptreResults_local.per_element_output,
                PerturboResults.global_per_guide_output, PerturboResults.global_per_element_output
            )
            calibration_output = Calibration.calibration_output
        } else {
            calibration_output = calibration_placeholder
        }
    } else {
        error("Invalid INFERENCE_method: ${params.INFERENCE_method}. Valid options: sceptre, perturbo, sceptre,perturbo, default")
    }

    emit:
    inference_mudata = FinalInference
    inference_mode = mode
    calibration_output = calibration_output

}
