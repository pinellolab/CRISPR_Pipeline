#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/crispr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/crispr
    Website: https://nf-co.re/crispr
    Slack  : https://nfcore.slack.com/channels/crispr
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CRISPR_PIPELINE }  from './workflows/crispr_pipeline'
include { inference_pipeline } from './subworkflows/local/inference_pipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_crispr_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_crispr_pipeline'
include { PACKAGE_EXISTING_IGVF_RESULTS } from './subworkflows/local/utils_nfcore_crispr_pipeline'
include { skipGTFDownload } from './modules/local/skipGTFDownload'
include { downloadGTF } from './modules/local/downloadGTF'
include { remove_clonal_cells } from './modules/local/remove_clonal_cells'
include { stage_checkpoint_dashboard; render_continuation_dashboard } from './modules/local/checkpoint_dashboard'
include { prepare_external_mudata } from './modules/local/prepare_external_mudata'
include { evaluation_pipeline } from './subworkflows/local/evaluation_pipeline'
include { BALANCED_BENCHMARK } from './workflows/balanced_benchmark'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_CRISPR {

    take:
    samplesheet // channel: samplesheet read in from --input
    samplesheet_file

    main:
    //
    // WORKFLOW: Run pipeline
    //
    CRISPR_PIPELINE(
        samplesheet,
        samplesheet_file
        )
}

//
// WORKFLOW: Run default inference starting from an existing MuData file
//
workflow INFERENCE_FROM_MUDATA {

    main:
    if (!params.INFERENCE_input_mudata) {
        error "INFERENCE_FROM_MUDATA requires --INFERENCE_input_mudata <path_to_h5mu>."
    }
    if (params.INFERENCE_method != 'default') {
        error "INFERENCE_FROM_MUDATA requires --INFERENCE_method 'default'."
    }
    if (params.INFERENCE_target_guide_pairing_strategy != 'default') {
        error "INFERENCE_FROM_MUDATA requires --INFERENCE_target_guide_pairing_strategy 'default'."
    }

    mudata_input = file(params.INFERENCE_input_mudata)
    if (!mudata_input.exists()) {
        error "INFERENCE_input_mudata file was not found: ${params.INFERENCE_input_mudata}"
    }

    if (file(params.REFERENCE_gtf_local_path).exists()) {
        GTF_Reference = skipGTFDownload(file(params.REFERENCE_gtf_local_path))
    }
    else {
        if (!params.REFERENCE_gtf_download_path) {
            error "REFERENCE_gtf_download_path is not set and REFERENCE_gtf_local_path does not exist."
        }
        GTF_Reference = downloadGTF(params.REFERENCE_gtf_download_path)
    }

    if (params.ENABLE_CLONE_REMOVAL) {
        CloneRemoval = remove_clonal_cells(mudata_input)
        mudata_for_inference = CloneRemoval.filtered_mudata
    } else {
        mudata_for_inference = mudata_input
    }

    Inference = inference_pipeline(
        mudata_for_inference,
        GTF_Reference.gencode_gtf
    )

    emit:
    inference_mudata = Inference.inference_mudata
}

// A second, explicit invocation after the successful pre-inference checkpoint.
// Neither mapping nor clone removal is repeated on this path.
workflow CONTINUE_FROM_CHECKPOINT {
    main:
    if (!params.ANALYSIS_parent_manifest) {
        error "ANALYSIS_mode='continue' requires --ANALYSIS_parent_manifest."
    }
    if (params.INFERENCE_method != 'default' || params.INFERENCE_target_guide_pairing_strategy != 'default') {
        error "The first checkpoint pilot requires default inference and pairing settings."
    }
    parent_manifest = file(params.ANALYSIS_parent_manifest)
    if (!parent_manifest.exists()) {
        error "Parent manifest not found: ${parent_manifest}"
    }
    parent_mudata = parent_manifest.parent.resolve('preinference_mudata.h5mu')
    if (!parent_mudata.exists()) {
        error "Parent MuData not found beside the manifest: ${parent_mudata}"
    }
    checkpoint_dir = parent_manifest.toRealPath().parent
    parent_outdir = checkpoint_dir.parent
    parent_sources = []
    guide_report = []
    if (checkpoint_dir.parent?.name == 'preinference_checkpoint') {
        parent_outdir = checkpoint_dir.parent.parent
        continuation_outdir = file(params.outdir)
        if (continuation_outdir.exists() && continuation_outdir.toRealPath() == parent_outdir) {
            error "Checkpoint continuation requires a different output directory from the parent run."
        }
        def validation = parent_outdir.resolve("pipeline_info/${checkpoint_dir.name}/guide_metadata.validation.json")
        guide_report = validation.exists() ? [validation] : []
        if (!checkpoint_dir.resolve('checkpoint_qc').exists()) {
            parent_sources = [
                'measurement_set_qc', 'guide_mapping_qc', 'guide_assignment_qc',
                'postconcat_embedding_qc', 'clone_removal', 'hashing_qc', 'sequencing_saturation',
                'pipeline_outputs/seqspeccheck', 'seqspeccheck'
            ].collect { parent_outdir.resolve(it) }.findAll { it.exists() }.unique { it.name }
        }
    }
    def snapshot = checkpoint_dir.resolve('checkpoint_qc')
    if (snapshot.exists()) {
        parent_sources = [snapshot]
    }

    edited_mudata = params.ANALYSIS_input_mudata ? file(params.ANALYSIS_input_mudata) : parent_mudata
    keep_cells = file(params.ANALYSIS_keep_cells_tsv ?: "${projectDir}/dummy_hash.txt")
    covariates = file(params.ANALYSIS_covariates_tsv ?: "${projectDir}/dummy_hash.txt")
    if (!edited_mudata.exists() || !keep_cells.exists() || !covariates.exists()) {
        error "One or more checkpoint continuation inputs do not exist."
    }

    if (file(params.REFERENCE_gtf_local_path).exists()) {
        GTF_Reference = skipGTFDownload(file(params.REFERENCE_gtf_local_path))
    } else {
        if (!params.REFERENCE_gtf_download_path) {
            error "REFERENCE_gtf_download_path is unset and local GTF does not exist."
        }
        GTF_Reference = downloadGTF(params.REFERENCE_gtf_download_path)
    }
    External = prepare_external_mudata(
        parent_manifest, parent_mudata, edited_mudata,
        GTF_Reference.gencode_gtf, keep_cells, covariates
    )
    CheckpointDashboard = stage_checkpoint_dashboard(External.report, parent_sources, guide_report)
    Inference = inference_pipeline(External.mudata, GTF_Reference.gencode_gtf)
    Evaluation = evaluation_pipeline(GTF_Reference.gencode_gtf, Inference.inference_mudata)
    render_continuation_dashboard(
        CheckpointDashboard.dashboard, Inference.inference_mudata, Evaluation.evaluation_output_dir,
        Evaluation.control_output_dir, Inference.calibration_output,
        Inference.cell_filter_report.ifEmpty { file("${projectDir}/dummy_hash.txt") }
    )

    emit:
    analysis_input = External.mudata
    external_report = External.report
    inference_mudata = Inference.inference_mudata
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    if (params.OUTPUT_IGVF_PACKAGE && params.ANALYSIS_mode == 'full' && params.INFERENCE_method == 'default' && !params.INFERENCE_input_mudata && !params.DASHBOARD_PLUGIN?.enabled) {
        error 'OUTPUT_IGVF_PACKAGE requires the bundled nf-crispr-dashboard plugin (DASHBOARD_PLUGIN.enabled=true).'
    }
    if (params.OUTPUT_IGVF_PACKAGE && params.ANALYSIS_mode == 'full' && params.INFERENCE_method != 'default') {
        log.warn('IGVF upload package skipped: four canonical local/global tables require INFERENCE_method=default.')
    }
    if (!(params.ANALYSIS_mode in ['checkpoint', 'continue', 'full'])) {
        error "ANALYSIS_mode must be checkpoint, continue, or full."
    }
    // A supplied MuData bypasses FASTQ mapping and guide assignment, allowing
    // reproducible inference recovery from an already validated dataset.
    if (params.BALANCED_source_mudata) {
        if (params.ANALYSIS_mode != 'full' || params.INFERENCE_input_mudata || params.INFERENCE_method != 'default' || params.INFERENCE_target_guide_pairing_strategy != 'default') {
            error 'Balanced benchmarks require ANALYSIS_mode=full, default inference/pairing, and no INFERENCE_input_mudata.'
        }
        if (!params.BALANCED_parent_inference || !params.BALANCED_dataset || params.BALANCED_cells <= 0 || params.BALANCED_rna_umis <= 0) {
            error 'Balanced benchmarks require parent inference, dataset label, positive cell and RNA UMI targets.'
        }
        PIPELINE_INITIALISATION(
            params.version, params.validate_params, params.monochrome_logs, args,
            params.outdir, params.input, params.DEMO_MODE, params.ENABLE_DATA_HASHING
        )
        BALANCED_BENCHMARK()
    }
    else if (params.ANALYSIS_mode == 'continue') {
        CONTINUE_FROM_CHECKPOINT()
    }
    else if (params.INFERENCE_input_mudata) {
        INFERENCE_FROM_MUDATA()
    }
    else {
        //
        // Run initialisation tasks
        //
        PIPELINE_INITIALISATION (
            params.version,
            params.validate_params,
            params.monochrome_logs,
            args,
            params.outdir,
            params.input, //updated_samplesheet
            params.DEMO_MODE,
            params.ENABLE_DATA_HASHING
        )

        if (params.DEBUG_VAR) {
            PIPELINE_INITIALISATION.out.samplesheet.view { meta, fastqs ->
                "Sample: ${meta.id}, Single-end: ${meta.single_end}, Files: ${fastqs}"
            }
        }

        //
        // WORKFLOW: Run main workflow
        //
        NFCORE_CRISPR (
            PIPELINE_INITIALISATION.out.samplesheet,
            PIPELINE_INITIALISATION.out.samplesheet_file
        )
    }

    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        params.DEMO_MODE,

    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
