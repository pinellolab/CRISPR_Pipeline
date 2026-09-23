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
include { remove_clonal_cells } from '../../modules/local/remove_clonal_cells'
include { sequencing_saturation } from '../../modules/local/sequencing_saturation'
include { filter_hto_post_clone } from '../../modules/local/filter_hto_post_clone'
include { postconcat_embedding_qc as embedding_before_clone } from '../../modules/local/postconcat_embedding_qc'
include { postconcat_embedding_qc as embedding_after_clone } from '../../modules/local/postconcat_embedding_qc'

// Import hashing-specific modules
include { CreateMuData } from '../../modules/local/CreateMuData'
include { demultiplex } from '../../modules/local/demultiplex'
include { filter_hashing } from '../../modules/local/filter_hashing'
include { hashing_concat } from '../../modules/local/hashing_concat'
include { evaluation_pipeline } from '../../subworkflows/local/evaluation_pipeline'
include { tf_benchmark } from '../../modules/local/tf_benchmark'

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
    ch_guide = ch_samples.filter { meta, _fastqs -> meta.modality == 'grna' }
    ch_hash = ch_samples.filter { meta, _fastqs -> meta.modality == 'hash' }
    if (params.DEBUG_VAR) {
        ch_guide.view()
    }

    ch_rna_seqspec = ch_rna
        .map { meta, _fastqs -> file(meta.seqspec) }
        .unique()
        .first()

    ch_guide_seqspec = ch_guide
        .map { meta, _fastqs -> file(meta.seqspec) }
        .unique()
        .first()

    ch_hash_seqspec = ch_hash
        .map { meta, _fastqs -> file(meta.seqspec) }
        .unique()
        .first()

    if (params.DEBUG_VAR) {
        ch_rna_seqspec.view { "RNA seqspec: $it" }
        ch_guide_seqspec.view { "Guide seqspec: $it" }
        ch_hash_seqspec.view { "Hash seqspec: $it" }
    }

    // barcode_onlist
    ch_barcode_onlist = ch_rna
        .map { meta, _fastqs -> file(meta.barcode_onlist) }
        .unique()
        .first()

    //guide_design
    ch_guide_design = ch_guide
        .map { meta, _fastqs -> file(meta.guide_design) }
        .unique()
        .first()

    //barcode_hashtag_map
    ch_barcode_hashtag_map = ch_hash
        .map { meta, _fastqs -> file(meta.barcode_hashtag_map) }
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
        mapping_rna_pipeline.out.trans_out_dir,
        prepare_mapping_pipeline.out.parsed_covariate_file
    )

    if (params.ENABLE_SEQUENCING_SATURATION) {
        SaturationQC = sequencing_saturation(
            mapping_rna_pipeline.out.ks_transcripts_out_dir_collected,
            Preprocessing.filtered_anndata_rna,
            mapping_rna_pipeline.out.transcriptome_t2g,
            prepare_mapping_pipeline.out.parsed_covariate_file
        )
        saturation_qc_dir = SaturationQC.saturation_qc
    } else {
        saturation_qc_dir = file("${workflow.projectDir}/assets/saturation_qc_empty")
    }

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

        hashing_demux_anndata_collected = Demultiplex.hashing_demux_anndata
            .collect()
            .map { files -> files.sort { a, b -> a.toString() <=> b.toString() } }
        hashing_demux_unfiltered_anndata_collected = Demultiplex.hashing_demux_unfiltered_anndata
            .collect()
            .map { files -> files.sort { a, b -> a.toString() <=> b.toString() } }

        Hashing_Concat = hashing_concat(hashing_demux_anndata_collected, hashing_demux_unfiltered_anndata_collected)

        // Create MuData with hashing
        MergeMuData = CreateMuData(
            Preprocessing.filtered_anndata_rna,
            mapping_guide_pipeline.out.concat_anndata_guide,
            ch_guide_design,
            Preprocessing.gencode_gtf,
            params.Multiplicity_of_infection,
            params.GUIDE_ASSIGNMENT_capture_method,
            params.REFERENCE_restrict_genes_to_gtf,
            // Preserve classes until GEX/guide QC; HTO filtering then establishes
            // the common raw-count population for clone calling and embedding QC.
            Hashing_Concat.concatenated_hashing_unfiltered_demux
        )

        // Shared processing pipeline
        GuideAssignment = guide_assignment_pipeline(MergeMuData.mudata)
        // HTO support is assessed on GEX/guide-qualified cells before the parallel branches.
        HTOFilter = filter_hto_post_clone(GuideAssignment.concat_mudata)
        qualified_mudata = HTOFilter.filtered_mudata
        if (params.ENABLE_POSTCONCAT_EMBEDDING_QC) {
            BeforeEmbedding = embedding_before_clone(qualified_mudata, 'before_clone', file("${projectDir}/assets/cell_cycle/regev_lab_cell_cycle_genes.txt"))
        }
        if (params.ENABLE_CLONE_REMOVAL) {
            CloneRemoval = remove_clonal_cells(qualified_mudata)
            mudata_before_hto = CloneRemoval.filtered_mudata
            clone_qc_dir = CloneRemoval.clone_qc
        } else {
            mudata_before_hto = qualified_mudata
            clone_qc_dir = file("${workflow.projectDir}/assets/clone_qc_empty")
        }
        mudata_for_inference = mudata_before_hto
        if (params.ENABLE_POSTCONCAT_EMBEDDING_QC) {
            if (params.ENABLE_CLONE_REMOVAL) {
                AfterEmbedding = embedding_after_clone(mudata_before_hto, 'after_clone', file("${projectDir}/assets/cell_cycle/regev_lab_cell_cycle_genes.txt"))
                mudata_for_inference = AfterEmbedding.filtered_mudata
                embedding_dirs = BeforeEmbedding.qc_dir.mix(AfterEmbedding.qc_dir).collect()
            } else {
                mudata_for_inference = BeforeEmbedding.filtered_mudata
                embedding_dirs = BeforeEmbedding.qc_dir.collect()
            }
        } else {
            embedding_dirs = file("${workflow.projectDir}/assets/embedding_qc_empty")
        }
        Inference = inference_pipeline(mudata_for_inference, Preprocessing.gencode_gtf)

        evaluation_pipeline (
            Preprocessing.gencode_gtf,
            Inference.inference_mudata
            )

        AdditionalQC = additional_qc_plots(
            Inference.inference_mudata,
            clone_qc_dir,
            saturation_qc_dir,
            GuideAssignment.guide_assignment_qc,
            HTOFilter.hto_qc,
            embedding_dirs
        )

        if (params.ENABLE_BENCHMARK) {
            Benchmark = tf_benchmark(
                Inference.inference_mudata,
                Preprocessing.gencode_gtf,
                file(params.ENCODE_BED_DIR),
                params.DEMO_MODE
            )
            benchmark_output_dir = Benchmark.benchmark_output
        } else {
            benchmark_output_dir = Channel.value(file("${workflow.projectDir}/assets/benchmark_empty"))
        }

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
            HTOFilter.filtered_hashing,
            HTOFilter.unfiltered_hashing,
            Inference.inference_mudata,
            AdditionalQC.additional_qc,
            Preprocessing.figures_dir,
            evaluation_pipeline.out.evaluation_output_dir,
            evaluation_pipeline.out.control_output_dir,
            benchmark_output_dir

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
            params.REFERENCE_restrict_genes_to_gtf,
            file("${workflow.projectDir}/dummy_hash.txt") // Dummy file for hashing parameter when not using hashing
        )

        // Scrublet now runs independently for each RNA measurement set before
        // concatenation, so the assembled MuData is already doublet-filtered.
        mudata_for_processing = MergeMuData.mudata

        // Shared processing pipeline
        GuideAssignment = guide_assignment_pipeline(mudata_for_processing)
        if (params.ENABLE_POSTCONCAT_EMBEDDING_QC) {
            BeforeEmbedding = embedding_before_clone(GuideAssignment.concat_mudata, 'before_clone', file("${projectDir}/assets/cell_cycle/regev_lab_cell_cycle_genes.txt"))
        }
        if (params.ENABLE_CLONE_REMOVAL) {
            CloneRemoval = remove_clonal_cells(GuideAssignment.concat_mudata)
            mudata_for_inference = CloneRemoval.filtered_mudata
            clone_qc_dir = CloneRemoval.clone_qc
        } else {
            mudata_for_inference = GuideAssignment.concat_mudata
            clone_qc_dir = file("${workflow.projectDir}/assets/clone_qc_empty")
        }
        if (params.ENABLE_POSTCONCAT_EMBEDDING_QC) {
            if (params.ENABLE_CLONE_REMOVAL) {
                AfterEmbedding = embedding_after_clone(mudata_for_inference, 'after_clone', file("${projectDir}/assets/cell_cycle/regev_lab_cell_cycle_genes.txt"))
                mudata_for_inference = AfterEmbedding.filtered_mudata
                embedding_dirs = BeforeEmbedding.qc_dir.mix(AfterEmbedding.qc_dir).collect()
            } else {
                mudata_for_inference = BeforeEmbedding.filtered_mudata
                embedding_dirs = BeforeEmbedding.qc_dir.collect()
            }
        } else {
            embedding_dirs = file("${workflow.projectDir}/assets/embedding_qc_empty")
        }
        Inference = inference_pipeline(mudata_for_inference, Preprocessing.gencode_gtf)

        evaluation_pipeline (
            Preprocessing.gencode_gtf,
            Inference.inference_mudata
            )

        AdditionalQC = additional_qc_plots(
            Inference.inference_mudata,
            clone_qc_dir,
            saturation_qc_dir,
            GuideAssignment.guide_assignment_qc,
            file("${workflow.projectDir}/assets/hto_qc_empty"),
            embedding_dirs
        )

        if (params.ENABLE_BENCHMARK) {
            Benchmark = tf_benchmark(
                Inference.inference_mudata,
                Preprocessing.gencode_gtf,
                file(params.ENCODE_BED_DIR),
                params.DEMO_MODE
            )
            benchmark_output_dir = Benchmark.benchmark_output
        } else {
            benchmark_output_dir = Channel.value(file("${workflow.projectDir}/assets/benchmark_empty"))
        }

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
            evaluation_pipeline.out.control_output_dir,
            benchmark_output_dir
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
