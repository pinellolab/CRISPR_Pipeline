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
include { prepare_external_mudata } from './modules/local/prepare_external_mudata'
include { evaluation_pipeline } from './subworkflows/local/evaluation_pipeline'
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
    Inference = inference_pipeline(External.mudata, GTF_Reference.gencode_gtf)
    evaluation_pipeline(GTF_Reference.gencode_gtf, Inference.inference_mudata)

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
    if (params.ANALYSIS_mode == 'continue') {
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
