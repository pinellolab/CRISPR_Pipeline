/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQSPECCHECK                   } from '../modules/local/seqspeccheck'
include { PREPARECOVARIATE               } from '../modules/local/preparecovariate'
include { SEQSPECPARSER as SEQSPECPARSER_RNA       } from '../modules/local/seqspecparser'
include { SEQSPECPARSER as SEQSPECPARSER_GUIDE     } from '../modules/local/seqspecparser'
include { SEQSPECPARSER as SEQSPECPARSER_HASH      } from '../modules/local/seqspecparser'
include { EXTRACTKBCHEMISTRY as EXTRACT_RNA_CHEM   } from '../modules/local/extractkbchemistry'
include { EXTRACTKBCHEMISTRY as EXTRACT_GUIDE_CHEM } from '../modules/local/extractkbchemistry'
include { EXTRACTKBCHEMISTRY as EXTRACT_HASH_CHEM  } from '../modules/local/extractkbchemistry'
include { PREPAREGUIDECHEMISTRY          } from '../modules/local/prepareguidechemistry'
include { CREATEGUIDEREF                 } from '../modules/local/createguideref'
include { CREATEHASHINGREF               } from '../modules/local/createhashingref'
include { ANNDATACONCAT as ANNDATACONCAT_RNA   } from '../modules/local/anndataconcat'
include { ANNDATACONCAT as ANNDATACONCAT_GUIDE } from '../modules/local/anndataconcat'
include { ANNDATACONCAT as ANNDATACONCAT_HASH  } from '../modules/local/anndataconcat'
include { QUANT_TRANSCRIPTOME as QUANT_RNA   } from '../subworkflows/local/quant_transcriptome'
include { QUANT_TRANSCRIPTOME as QUANT_GUIDE } from '../subworkflows/local/quant_transcriptome'
include { QUANT_TRANSCRIPTOME as QUANT_HASH  } from '../subworkflows/local/quant_transcriptome'
include { preprocessing_pipeline         } from '../subworkflows/local/preprocessing_pipeline'
include { guide_assignment_pipeline      } from '../subworkflows/local/guide_assignment_pipeline'
include { inference_pipeline             } from '../subworkflows/local/inference_pipeline'
include { additional_qc_plots            } from '../modules/local/additional_qc_plots'
include { CreateMuData                   } from '../modules/local/CreateMuData'
include { doublets_scrub                 } from '../modules/local/doublets_scrub'
include { demultiplex                    } from '../modules/local/demultiplex'
include { filter_hashing                 } from '../modules/local/filter_hashing'
include { hashing_concat                 } from '../modules/local/hashing_concat'
include { evaluation_pipeline            } from '../subworkflows/local/evaluation_pipeline'
include { tf_benchmark                   } from '../modules/local/tf_benchmark'
include { dashboard_pipeline_HASHING     } from '../subworkflows/local/dashboard_pipeline_HASHING'
include { dashboard_pipeline             } from '../subworkflows/local/dashboard_pipeline'
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText         } from '../subworkflows/local/utils_nfcore_perturbseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PERTURBSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def ch_empty_file = channel.value([])
    def ch_rna_workflow = channel.value('standard')
    def ch_hash_workflow = channel.value('kite')
    def ch_transcriptome_index = channel.value(file(params.transcriptome_index, checkIfExists: true))
    def ch_transcriptome_t2g = channel.value(file(params.transcriptome_t2g, checkIfExists: true))

    // Parse the samplesheet and create channels for each modality
    def ch_samples = ch_samplesheet.map { meta, fastqs ->
        meta.modality = meta.modality.toLowerCase()
        [meta, fastqs]
    }

    def ch_rna   = ch_samples.filter { meta, _fastqs -> meta.modality == 'scrna' }
    def ch_guide = ch_samples.filter { meta, _fastqs -> meta.modality == 'grna' }
    def ch_hash  = ch_samples.filter { meta, _fastqs -> meta.modality == 'hash' }
    def ch_rna_prefix = ch_rna.map { meta, _reads -> "${meta.measurement_sets}_ks_transcripts_out" }
    def ch_guide_prefix = ch_guide.map { meta, _reads -> "${meta.measurement_sets}_ks_guide_out" }
    def ch_hash_prefix = ch_hash.map { meta, _reads -> "${meta.measurement_sets}_ks_hashing_out" }
    def ch_rna_for_count = ch_rna.map { meta, reads -> [meta + [id: meta.measurement_sets], reads] }
    def ch_guide_for_count = ch_guide.map { meta, reads -> [meta + [id: meta.measurement_sets], reads] }
    def ch_hash_for_count = ch_hash.map { meta, reads -> [meta + [id: meta.measurement_sets], reads] }

    def ch_rna_seqspec = ch_rna
        .map { meta, _fastqs -> file(meta.seqspec) }
        .unique()
        .first()

    def ch_guide_seqspec = ch_guide
        .map { meta, _fastqs -> file(meta.seqspec) }
        .unique()
        .first()

    def ch_hash_seqspec = ch_hash
        .map { meta, _fastqs -> file(meta.seqspec) }
        .unique()
        .first()

    // barcode_onlist
    def ch_barcode_onlist = ch_rna
        .map { meta, _fastqs -> file(meta.barcode_onlist) }
        .unique()
        .first()

    // guide_design
    def ch_guide_design = ch_guide
        .map { meta, _fastqs -> file(meta.guide_design) }
        .unique()
        .first()

    // barcode_hashtag_map
    def ch_barcode_hashtag_map = ch_hash
        .map { meta, _fastqs -> file(meta.barcode_hashtag_map) }
        .unique()
        .first()

    guide_seqSpecCheck = SEQSPECCHECK(
        ch_guide,
        ch_guide_design,
        'guide'
    )
    // TODO: add hashing seqSpecCheck version to ch_version

    // Run seqSpecCheck pipeline
    if (params.ENABLE_DATA_HASHING) {
        hash_seqSpecCheck = SEQSPECCHECK(
            ch_hash,
            ch_barcode_hashtag_map,
            'hashing'
        )
        // TODO: add hashing seqSpecCheck version to ch_version
    }

    def ch_covariate_list = ch_samples
        .map { meta, _reads ->
            [meta.measurement_sets]
        }
        .unique()
        .map { measurement_sets ->
            [batch: measurement_sets]
        }
        .collect()
        .map { items ->
            def sorted_covariates = items.sort { a, b -> a.batch.toString() <=> b.batch.toString() }
            [batch: sorted_covariates*.batch.flatten()]
        }

    ch_covariates = PREPARECOVARIATE(ch_covariate_list)
    ch_versions = ch_versions.mix(ch_covariates.out.versions.first())

    def SeqSpecRna = SEQSPECPARSER_RNA(
        ch_rna_seqspec,
        ch_barcode_onlist,
        'rna'
    )
    ch_versions = ch_versions.mix(SeqSpecRna.out.versions.first())

    def RawRnaChem = EXTRACT_RNA_CHEM(SeqSpecRna.out.parsed_seqspec)
    ch_versions = ch_versions.mix(RawRnaChem.out.versions.first())
    def ch_rna_technology = RawRnaChem.out.technology_file.map { file -> file.text.trim() }

    def QuantRna = QUANT_RNA(
        ch_rna_for_count,
        ch_transcriptome_index,
        ch_transcriptome_t2g,
        SeqSpecRna.out.barcode_file,
        ch_empty_file,
        ch_empty_file,
        ch_rna_technology,
        ch_rna_workflow,
        ch_rna_prefix,
    )
    ch_versions = ch_versions.mix(QuantRna.out.versions.first())

    def ks_transcripts_out_dir_collected = QuantRna.out.count
        .map { _meta, count_dir -> count_dir }
        .collect()
        .map { dirs -> dirs.sort { a, b -> a.getName() <=> b.getName() } }

    def ConcatRna = ANNDATACONCAT_RNA(
        ch_covariates.parsed_covariate_file,
        ks_transcripts_out_dir_collected
    )
    ch_versions = ch_versions.mix(ConcatRna.out.versions.first())

    def GuideRef = CREATEGUIDEREF(
        ch_guide_design,
        params.reverse_complement_guides,
        params.spacer_tag
    )
    ch_versions = ch_versions.mix(GuideRef.out.versions.first())

    def SeqSpecGuide = SEQSPECPARSER_GUIDE(
        ch_guide_seqspec,
        ch_barcode_onlist,
        'crispr'
    )
    ch_versions = ch_versions.mix(SeqSpecGuide.out.versions.first())

    def RawGuideChem = EXTRACT_GUIDE_CHEM(SeqSpecGuide.out.parsed_seqspec)
    ch_versions = ch_versions.mix(RawGuideChem.out.versions.first())

    def GuideChem = PREPAREGUIDECHEMISTRY(
        RawGuideChem.out.technology_file,
        params.is_10x3v3,
        params.spacer_tag
    )
    ch_versions = ch_versions.mix(GuideChem.out.versions.first())
    def ch_guide_technology = GuideChem.out.technology_file.map { file -> file.text.trim() }
    def ch_guide_workflow = GuideChem.out.workflow_file.map { file -> file.text.trim() }

    def QuantGuide = QUANT_GUIDE(
        ch_guide_for_count,
        GuideRef.out.guide_index,
        GuideRef.out.t2g_guide,
        SeqSpecGuide.out.barcode_file,
        ch_empty_file,
        ch_empty_file,
        ch_guide_technology,
        ch_guide_workflow,
        ch_guide_prefix,
    )
    ch_versions = ch_versions.mix(QuantGuide.out.versions.first())

    def ks_guide_out_dir_collected = QuantGuide.out.count
        .map { _meta, count_dir -> count_dir }
        .collect()
        .map { dirs -> dirs.sort { a, b -> a.getName() <=> b.getName() } }

    def ConcatGuide = ANNDATACONCAT_GUIDE(
        ch_covariates.parsed_covariate_file,
        ks_guide_out_dir_collected
    )
    ch_versions = ch_versions.mix(ConcatGuide.out.versions.first())

    // Common preprocessing for both workflows
    def Preprocessing = preprocessing_pipeline(
        ConcatRna.out.concat_anndata,
        QuantRna.out.count.map { _meta, count_dir -> count_dir }
    )

    def benchmark_output_dir

    if (params.ENABLE_DATA_HASHING) {
        def HashRef = CREATEHASHINGREF(ch_barcode_hashtag_map)
        ch_versions = ch_versions.mix(HashRef.out.versions.first())

        def SeqSpecHash = SEQSPECPARSER_HASH(
            ch_hash_seqspec,
            ch_barcode_onlist,
            'tag'
        )
        ch_versions = ch_versions.mix(SeqSpecHash.out.versions.first())

        def RawHashChem = EXTRACT_HASH_CHEM(SeqSpecHash.out.parsed_seqspec)
        ch_versions = ch_versions.mix(RawHashChem.out.versions.first())
        def ch_hash_technology = RawHashChem.out.technology_file.map { file -> file.text.trim() }

        def QuantHash = QUANT_HASH(
            ch_hash_for_count,
            HashRef.out.hashing_index,
            HashRef.out.t2g_hashing,
            SeqSpecHash.out.barcode_file,
            ch_empty_file,
            ch_empty_file,
            ch_hash_technology,
            ch_hash_workflow,
            ch_hash_prefix,
        )
        ch_versions = ch_versions.mix(QuantHash.out.versions.first())

        def ks_hashing_out_dir_collected = QuantHash.out.count
            .map { _meta, count_dir -> count_dir }
            .collect()
            .map { dirs -> dirs.sort { a, b -> a.getName() <=> b.getName() } }

        def ConcatHash = ANNDATACONCAT_HASH(
            ch_covariates.parsed_covariate_file,
            ks_hashing_out_dir_collected
        )
        ch_versions = ch_versions.mix(ConcatHash.out.versions.first())

        // Hashing-specific processing
        def Hashing_Filtered = filter_hashing(
            Preprocessing.filtered_anndata_rna,
            ConcatHash.out.concat_anndata
        )

        def Demultiplex = demultiplex(Hashing_Filtered.hashing_filtered_anndata.flatten())

        def hashing_demux_anndata_collected = Demultiplex.hashing_demux_anndata
            .collect()
            .map { files -> files.sort { a, b -> a.toString() <=> b.toString() } }
        def hashing_demux_unfiltered_anndata_collected = Demultiplex.hashing_demux_unfiltered_anndata
            .collect()
            .map { files -> files.sort { a, b -> a.toString() <=> b.toString() } }

        def Hashing_Concat = hashing_concat(hashing_demux_anndata_collected, hashing_demux_unfiltered_anndata_collected)

        // Create MuData with hashing
        def MergeMuData = CreateMuData(
            Preprocessing.filtered_anndata_rna,
            ConcatGuide.out.concat_anndata,
            ch_guide_design,
            Preprocessing.gencode_gtf,
            params.Multiplicity_of_infection,
            params.GUIDE_ASSIGNMENT_capture_method,
            Hashing_Concat.concatenated_hashing_demux
        )

        def GuideAssignment = guide_assignment_pipeline(MergeMuData.mudata)
        def Inference = inference_pipeline(GuideAssignment.concat_mudata, Preprocessing.gencode_gtf)

        evaluation_pipeline(
            Preprocessing.gencode_gtf,
            Inference.inference_mudata
        )

        def AdditionalQC = additional_qc_plots(Inference.inference_mudata)

        if (params.ENABLE_BENCHMARK) {
            def Benchmark = tf_benchmark(
                Inference.inference_mudata,
                Preprocessing.gencode_gtf,
                file(params.ENCODE_BED_DIR)
            )
            benchmark_output_dir = Benchmark.benchmark_output
        } else {
            benchmark_output_dir = channel.value(file("${workflow.projectDir}/assets/benchmark_empty"))
        }

        dashboard_pipeline_HASHING(
            guide_seqSpecCheck.out.seqSpecCheck_plots,
            guide_seqSpecCheck.out.position_table,
            hash_seqSpecCheck.out.seqSpecCheck_plots,
            hash_seqSpecCheck.out.position_table,
            Preprocessing.adata_rna,
            Preprocessing.filtered_anndata_rna,
            ks_transcripts_out_dir_collected,
            MergeMuData.adata_guide,
            ks_guide_out_dir_collected,
            Hashing_Filtered.adata_hashing,
            ks_hashing_out_dir_collected,
            Hashing_Concat.concatenated_hashing_demux,
            Hashing_Concat.concatenated_hashing_unfiltered_demux,
            Inference.inference_mudata,
            AdditionalQC.additional_qc,
            Preprocessing.figures_dir,
            evaluation_pipeline.out.evaluation_output_dir,
            evaluation_pipeline.out.control_output_dir,
            benchmark_output_dir
        )
    } else {
        // Create MuData without hashing
        def MergeMuData = CreateMuData(
            Preprocessing.filtered_anndata_rna,
            ConcatGuide.out.concat_anndata,
            ch_guide_design,
            Preprocessing.gencode_gtf,
            params.Multiplicity_of_infection,
            params.GUIDE_ASSIGNMENT_capture_method,
            file("${workflow.projectDir}/dummy_hash.txt")
        )

        def mudata_for_processing
        if (params.ENABLE_SCRUBLET ?: false) {
            def MuData_Doublets = doublets_scrub(MergeMuData.mudata)
            mudata_for_processing = MuData_Doublets.mudata_doublet
        } else {
            mudata_for_processing = MergeMuData.mudata
        }

        def GuideAssignment = guide_assignment_pipeline(mudata_for_processing)
        def Inference = inference_pipeline(GuideAssignment.concat_mudata, Preprocessing.gencode_gtf)

        evaluation_pipeline(
            Preprocessing.gencode_gtf,
            Inference.inference_mudata
        )

        def AdditionalQC = additional_qc_plots(Inference.inference_mudata)

        if (params.ENABLE_BENCHMARK) {
            def Benchmark = tf_benchmark(
                Inference.inference_mudata,
                Preprocessing.gencode_gtf,
                file(params.ENCODE_BED_DIR)
            )
            benchmark_output_dir = Benchmark.benchmark_output
        } else {
            benchmark_output_dir = channel.value(file("${workflow.projectDir}/assets/benchmark_empty"))
        }

        dashboard_pipeline(
            guide_seqSpecCheck.out.seqSpecCheck_plots,
            guide_seqSpecCheck.out.position_table,
            Preprocessing.adata_rna,
            Preprocessing.filtered_anndata_rna,
            ks_transcripts_out_dir_collected,
            MergeMuData.adata_guide,
            ks_guide_out_dir_collected,
            Inference.inference_mudata,
            AdditionalQC.additional_qc,
            Preprocessing.figures_dir,
            evaluation_pipeline.out.evaluation_output_dir,
            evaluation_pipeline.out.control_output_dir,
            benchmark_output_dir
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_perturbseq_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'perturbseq'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                                                    // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
