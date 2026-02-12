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
// Import modular subworkflows
include { preprocessing_pipeline } from '../../subworkflows/local/preprocessing_pipeline'
include { guide_assignment_pipeline } from '../../subworkflows/local/guide_assignment_pipeline'
include { inference_pipeline } from '../../subworkflows/local/inference_pipeline'
include { additional_qc_plots } from '../../modules/local/additional_qc_plots'

// Import hashing-specific modules
include { CreateMuData } from '../../modules/local/CreateMuData'
include { doublets_scrub } from '../../modules/local/doublets_scrub'
include { demultiplex } from '../../modules/local/demultiplex'
include { filter_hashing } from '../../modules/local/filter_hashing'
include { hashing_concat } from '../../modules/local/hashing_concat'
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
    ch_guide.view()
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

    //guide_design
    ch_guide_design = ch_guide
        .map { meta, _fastqs -> meta.guide_design }
        .unique()
        .first()

    //barcode_hashtag_map
    ch_barcode_hashtag_map = ch_hash
        .map { meta, _fastqs -> meta.barcode_hashtag_map }
        .unique()
        .first()

    // Run seqSpecCheck pipeline
    if (params.ENABLE_DATA_HASHING) {
        seqSpecCheck_pipeline_HASHING(ch_guide.first(), ch_hash.first(), ch_guide_design, ch_barcode_hashtag_map)
    } else {
        seqSpecCheck_pipeline(ch_guide.first(), ch_guide_design)
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
        ch_guide_design,
        prepare_mapping_pipeline.out.parsed_covariate_file,
        params.reverse_complement_guides,
        params.spacer_tag
        )

    // Common preprocessing for both workflows
    Preprocessing = preprocessing_pipeline(
        mapping_rna_pipeline.out.concat_anndata_rna,
        mapping_rna_pipeline.out.trans_out_dir
    )

    if (params.ENABLE_DATA_HASHING) {
        mapping_hashing_pipeline(
            ch_hash,
            ch_hash_seqspec,
            ch_barcode_onlist,
            ch_barcode_hashtag_map,
            prepare_mapping_pipeline.out.parsed_covariate_file
            )

        // Hashing-specific processing
        Hashing_Filtered = filter_hashing(
            Preprocessing.filtered_anndata_rna,
            mapping_hashing_pipeline.out.concat_anndata_hashing
        )

        Demultiplex = demultiplex(Hashing_Filtered.hashing_filtered_anndata.flatten())

        hashing_demux_anndata_collected = Demultiplex.hashing_demux_anndata.collect()
        hashing_demux_unfiltered_anndata_collected = Demultiplex.hashing_demux_unfiltered_anndata.collect()

        Hashing_Concat = hashing_concat(hashing_demux_anndata_collected, hashing_demux_unfiltered_anndata_collected)

        // Create MuData with hashing
        MergeMuData = CreateMuData(
            Preprocessing.filtered_anndata_rna,
            mapping_guide_pipeline.out.concat_anndata_guide,
            ch_guide_design,
            Preprocessing.gencode_gtf,
            params.Multiplicity_of_infection,
            params.GUIDE_ASSIGNMENT_capture_method,
            Hashing_Concat.concatenated_hashing_demux
        )

        // Shared processing pipeline
        GuideAssignment = guide_assignment_pipeline(MergeMuData.mudata)
        Inference = inference_pipeline(GuideAssignment.concat_mudata, Preprocessing.gencode_gtf)

        evaluation_pipeline (
            Preprocessing.gencode_gtf,
            Inference.inference_mudata
            )

        AdditionalQC = additional_qc_plots(
            Inference.inference_mudata
        )

        dashboard_pipeline_HASHING (
            seqSpecCheck_pipeline_HASHING.out.guide_seqSpecCheck_plots,
            seqSpecCheck_pipeline_HASHING.out.guide_position_table,
            seqSpecCheck_pipeline_HASHING.out.hashing_seqSpecCheck_plots,
            seqSpecCheck_pipeline_HASHING.out.hashing_position_table,
            Preprocessing.adata_rna,
            Preprocessing.filtered_anndata_rna,
            mapping_rna_pipeline.out.ks_transcripts_out_dir_collected,
            MergeMuData.adata_guide,
            mapping_guide_pipeline.out.ks_guide_out_dir_collected,
            Hashing_Filtered.adata_hashing,
            mapping_hashing_pipeline.out.ks_hashing_out_dir_collected,
            Hashing_Concat.concatenated_hashing_demux,
            Hashing_Concat.concatenated_hashing_unfiltered_demux,
            Inference.inference_mudata,
            AdditionalQC.additional_qc,
            Preprocessing.figures_dir,
            evaluation_pipeline.out.evaluation_output_dir,
            evaluation_pipeline.out.control_output_dir

            )
    }
    else {
        // Create MuData without hashing
        MergeMuData = CreateMuData(
            Preprocessing.filtered_anndata_rna,
            mapping_guide_pipeline.out.concat_anndata_guide,
            ch_guide_design,
            Preprocessing.gencode_gtf,
            params.Multiplicity_of_infection,
            params.GUIDE_ASSIGNMENT_capture_method,
            file('dummy_hash.txt') // Dummy file for hashing parameter when not using hashing
        )

        // Conditionally run scrublet based on ENABLE_SCRUBLET parameter (defaults to false)
        if (params.ENABLE_SCRUBLET ?: false) {
            MuData_Doublets = doublets_scrub(MergeMuData.mudata)
            mudata_for_processing = MuData_Doublets.mudata_doublet
        } else {
            mudata_for_processing = MergeMuData.mudata
        }

        // Shared processing pipeline
        GuideAssignment = guide_assignment_pipeline(mudata_for_processing)
        Inference = inference_pipeline(GuideAssignment.concat_mudata, Preprocessing.gencode_gtf)

        evaluation_pipeline (
            Preprocessing.gencode_gtf,
            Inference.inference_mudata
            )

        AdditionalQC = additional_qc_plots(
            Inference.inference_mudata
        )

        dashboard_pipeline (
            seqSpecCheck_pipeline.out.guide_seqSpecCheck_plots,
            seqSpecCheck_pipeline.out.guide_position_table,
            Preprocessing.adata_rna,
            Preprocessing.filtered_anndata_rna,
            mapping_rna_pipeline.out.ks_transcripts_out_dir_collected,
            MergeMuData.adata_guide,
            mapping_guide_pipeline.out.ks_guide_out_dir_collected,
            Inference.inference_mudata,
            AdditionalQC.additional_qc,
            Preprocessing.figures_dir,
            evaluation_pipeline.out.evaluation_output_dir,
            evaluation_pipeline.out.control_output_dir
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
