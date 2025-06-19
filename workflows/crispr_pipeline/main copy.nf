/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { seqSpecCheck_pipeline } from '../../subworkflows/local/seqSpecCheck_pipeline'
include { seqSpecCheck_pipeline_HASHING } from '../../subworkflows/local/seqSpecCheck_pipeline_HASHING'
include { prepare_mapping_pipeline } from '../../subworkflows/local/prepare_mapping_pipeline'
include { mapping_rna_pipeline } from '../../subworkflows/local/mapping_rna_pipeline'
include { mapping_guide_pipeline } from '../../subworkflows/local/mapping_guide_pipeline'
include { mapping_hashing_pipeline } from '../../subworkflows/local/mapping_hashing_pipeline'
include { process_mudata_pipeline_HASHING } from '../../subworkflows/local/process_mudata_pipeline_HASHING'
include { process_mudata_pipeline } from '../../subworkflows/local/process_mudata_pipeline'
include { evaluation_pipeline } from '../../subworkflows/local/evaluation_pipeline'
include { dashboard_pipeline_HASHING } from '../../subworkflows/local/dashboard_pipeline_HASHING'
include { dashboard_pipeline } from '../../subworkflows/local/dashboard_pipeline'

include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_crispr_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW FOR CRISPR PERTURBED-SEQ PIPELINE
//

workflow CRISPR_PIPELINE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()

    // Parse the samplesheet and create channels for each modality
    ch_samples = ch_samplesheet.map { meta, fastqs ->
        // Modify meta to include all necessary information
        meta.modality = meta.modality.toLowerCase()
        [meta, fastqs]
    }

    ch_rna = ch_samples.filter { meta, _fastqs -> meta.modality == 'scrna' }
    //ch_rna.view()
    ch_guide = ch_samples.filter { meta, _fastqs -> meta.modality == 'grna' }
    //ch_guide.view()
    ch_hash = ch_samples.filter { meta, _fastqs -> meta.modality == 'hash' }
    //ch_hash.view()

    ch_rna_seqspec = ch_rna
        .map { meta, _fastqs -> meta.seqspec }
        .unique()
        .first()

    ch_guide_seqspec = ch_guide
        .map { meta, _fastqs -> meta.seqspec }
        .unique()
        .first()

    ch_hash_seqspec = ch_hash
        .map { meta, _fastqs -> meta.seqspec }
        .unique()
        .first()

    // View the seqspec paths
    ch_rna_seqspec.view { "RNA seqspec: $it" }
    ch_guide_seqspec.view { "Guide seqspec: $it" }
    ch_hash_seqspec.view { "Hash seqspec: $it" }

    // barcode_onlist
    ch_barcode_onlist = ch_rna
        .map { meta, _fastqs -> meta.barcode_onlist }
        .unique()
        .first()

    // Run seqSpecCheck pipeline
    if (params.ENABLE_DATA_HASHING == "true") {
        seqSpecCheck_pipeline_HASHING(ch_guide.first(), ch_hash.first())
    } else {
        seqSpecCheck_pipeline(ch_guide.first())
    }

    prepare_mapping_pipeline(ch_samples)

    // Run mapping pipelines for each modality
    mapping_rna_pipeline(
            ch_rna,
            ch_rna_seqspec,
            ch_barcode_onlist,
            prepare_mapping_pipeline.out.parsed_covariate_file
        )

    mapping_guide_pipeline(
        ch_guide,
        ch_guide_seqspec,
        ch_barcode_onlist,
        prepare_mapping_pipeline.out.parsed_covariate_file
        )

    if (params.ENABLE_DATA_HASHING == "true") {
        mapping_hashing_pipeline(
            ch_hash,
            ch_hash_seqspec,
            ch_barcode_onlist,
            prepare_mapping_pipeline.out.parsed_covariate_file
            )

        process_mudata_pipeline_HASHING(
            mapping_rna_pipeline.out.concat_anndata_rna,
            mapping_rna_pipeline.out.trans_out_dir,
            mapping_guide_pipeline.out.concat_anndata_guide,
            mapping_guide_pipeline.out.guide_out_dir,
            mapping_hashing_pipeline.out.concat_anndata_hashing,
            mapping_hashing_pipeline.out.hashing_out_dir,
            prepare_mapping_pipeline.out.covariate_string
            )

        evaluation_pipeline (
            process_mudata_pipeline_HASHING.out.gencode_gtf,
            process_mudata_pipeline_HASHING.out.inference_mudata
            )

        dashboard_pipeline_HASHING (
            seqSpecCheck_pipeline_HASHING.out.guide_seqSpecCheck_plots,
            seqSpecCheck_pipeline_HASHING.out.guide_position_table,
            seqSpecCheck_pipeline_HASHING.out.hashing_seqSpecCheck_plots,
            seqSpecCheck_pipeline_HASHING.out.hashing_position_table,
            process_mudata_pipeline_HASHING.out.adata_rna,
            process_mudata_pipeline_HASHING.out.filtered_anndata_rna,
            mapping_rna_pipeline.out.ks_transcripts_out_dir_collected,
            process_mudata_pipeline_HASHING.out.adata_guide,
            mapping_guide_pipeline.out.ks_guide_out_dir_collected,
            process_mudata_pipeline_HASHING.out.adata_hashing,
            mapping_hashing_pipeline.out.ks_hashing_out_dir_collected,
            process_mudata_pipeline_HASHING.out.adata_demux,
            process_mudata_pipeline_HASHING.out.adata_unfiltered_demux,
            process_mudata_pipeline_HASHING.out.inference_mudata,
            process_mudata_pipeline_HASHING.out.figures_dir,
            evaluation_pipeline.out.evaluation_output_dir
            )
    }
    else {
        process_mudata_pipeline(
            mapping_rna_pipeline.out.concat_anndata_rna,
            mapping_rna_pipeline.out.trans_out_dir,
            mapping_guide_pipeline.out.concat_anndata_guide,
            mapping_guide_pipeline.out.guide_out_dir,
            prepare_mapping_pipeline.out.covariate_string
            )
        evaluation_pipeline (
            process_mudata_pipeline.out.gencode_gtf,
            process_mudata_pipeline.out.inference_mudata
            )

        dashboard_pipeline (
            seqSpecCheck_pipeline.out.guide_seqSpecCheck_plots,
            seqSpecCheck_pipeline.out.guide_position_table,
            process_mudata_pipeline.out.adata_rna,
            process_mudata_pipeline.out.filtered_anndata_rna,
            mapping_rna_pipeline.out.ks_transcripts_out_dir_collected,
            process_mudata_pipeline.out.adata_guide,
            mapping_guide_pipeline.out.ks_guide_out_dir_collected,
            process_mudata_pipeline.out.inference_mudata,
            process_mudata_pipeline.out.figures_dir,
            evaluation_pipeline.out.evaluation_output_dir
            )
    }



    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
