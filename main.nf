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
include { skipGTFDownload } from './modules/local/skipGTFDownload'
include { downloadGTF } from './modules/local/downloadGTF'
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

    main:
    //
    // WORKFLOW: Run pipeline
    //
    CRISPR_PIPELINE(
        samplesheet
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

    Inference = inference_pipeline(
        mudata_input,
        GTF_Reference.gencode_gtf
    )

    emit:
    inference_mudata = Inference.inference_mudata
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input //updated_samplesheet
    )

    PIPELINE_INITIALISATION.out.samplesheet.view { meta, fastqs ->
        "Sample: ${meta.id}, Single-end: ${meta.single_end}, Files: ${fastqs}"
    }

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_CRISPR (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,

    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
