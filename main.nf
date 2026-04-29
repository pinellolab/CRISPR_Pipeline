#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/perturbseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/perturbseq
    Website: https://nf-co.re/perturbseq
    Slack  : https://nfcore.slack.com/channels/perturbseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PERTURBSEQ  } from './workflows/perturbseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_perturbseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_perturbseq_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_perturbseq_pipeline'
include { PREPAREREFERENCERESOURCES } from './subworkflows/local/preparereferenceresources'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.reference_fasta = params.reference_fasta ?: getGenomeAttribute('fasta')
params.reference_gtf   = params.reference_gtf   ?: getGenomeAttribute('gtf')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PERTURBSEQ {

    take:
    samplesheet // channel: samplesheet read in from --input
    transcriptome_index
    transcriptome_t2g
    reference_gtf
    reference_versions

    main:

    //
    // WORKFLOW: Run pipeline
    //
    PERTURBSEQ (
        samplesheet,
        transcriptome_index,
        transcriptome_t2g,
        reference_gtf,
        reference_versions,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = PERTURBSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.help,
        params.help_full,
        params.show_hidden
    )

    def ch_input = file(params.input, checkIfExists: true)
    def ch_samplesheet = channel
        .fromPath(ch_input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: "sample_${row.file_modality}_${row.measurement_sets}",
                measurement_sets: row.measurement_sets,
                modality: (row.file_modality ?: '').toLowerCase(),
                seqspec: row.seqspec,
                barcode_onlist: row.barcode_onlist,
                guide_design: row.guide_design,
                barcode_hashtag_map: row.barcode_hashtag_map,
            ]
            def reads = [file(row.R1_path, checkIfExists: true)]
            if (row.R2_path) {
                reads += file(row.R2_path, checkIfExists: true)
            }
            [meta, reads]
        }

    def ch_rna_workflow = channel.value('standard')
    def ch_reference_fasta = params.reference_fasta
        ? channel.value(file(params.reference_fasta, checkIfExists: true))
        : channel.empty()
    def ch_reference_gtf = params.reference_gtf
        ? channel.value(file(params.reference_gtf, checkIfExists: true))
        : channel.empty()

    def ReferenceResources = PREPAREREFERENCERESOURCES(
        ch_reference_fasta,
        ch_reference_gtf,
        ch_rna_workflow
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PERTURBSEQ (
        ch_samplesheet,
        ReferenceResources.out.transcriptome_index,
        ReferenceResources.out.transcriptome_t2g,
        ReferenceResources.out.reference_gtf,
        ReferenceResources.out.versions
    )
    // //
    // // SUBWORKFLOW: Run completion tasks
    // //
    // PIPELINE_COMPLETION (
    //     params.email,
    //     params.email_on_fail,
    //     params.plaintext_email,
    //     params.outdir,
    //     params.monochrome_logs,
    //     NFCORE_PERTURBSEQ.out.multiqc_report
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
