include { balance_mudata } from '../../modules/local/balance_mudata'
include { balanced_normalization } from '../../modules/local/balanced_normalization'
include { balanced_dashboard } from '../../modules/local/balanced_dashboard'
include { skipGTFDownload } from '../../modules/local/skipGTFDownload'
include { inference_pipeline } from '../../subworkflows/local/inference_pipeline'
include { evaluation_pipeline } from '../../subworkflows/local/evaluation_pipeline'
include { additional_qc_plots } from '../../modules/local/additional_qc_plots'
include { tf_benchmark } from '../../modules/local/tf_benchmark'
include { collect_benchmark_outputs } from '../../modules/local/collect_benchmark_outputs'
include { check_igvf_dashboard } from '../../modules/local/check_igvf_dashboard'

workflow BALANCED_BENCHMARK {
    main:
    source = file(params.BALANCED_source_mudata, checkIfExists: true)
    parent = file(params.BALANCED_parent_inference, checkIfExists: true)
    gtf = file(params.REFERENCE_gtf_local_path, checkIfExists: true)
    GTF = skipGTFDownload(gtf)
    Balanced = balance_mudata(source, parent)
    Normalized = balanced_normalization(Balanced.full_mudata, file("${projectDir}/assets/cell_cycle/regev_lab_cell_cycle_genes.txt"))
    Inference = inference_pipeline(Balanced.inference_input, GTF.gencode_gtf)
    Evaluation = evaluation_pipeline(GTF.gencode_gtf, Inference.inference_mudata)
    Additional = additional_qc_plots(
        Inference.inference_mudata,
        file("${projectDir}/assets/clone_qc_empty"),
        file("${projectDir}/assets/saturation_qc_empty"),
        file("${projectDir}/assets/clone_qc_empty"),
        file("${projectDir}/assets/hto_qc_empty"),
        Normalized.qc_dir.collect(),
        file("${projectDir}/assets/clone_qc_empty")
    )
    if (params.ENABLE_BENCHMARK) {
        Benchmark = tf_benchmark(Inference.inference_mudata, GTF.gencode_gtf, file(params.ENCODE_BED_DIR), params.DEMO_MODE)
        benchmark_dir = Benchmark.benchmark_output
    } else {
        benchmark_dir = Channel.value(file("${projectDir}/assets/benchmark_empty"))
    }
    Collected = collect_benchmark_outputs(benchmark_dir, Inference.calibration_output)
    Dashboard = balanced_dashboard(Inference.inference_mudata, Balanced.report, Additional.additional_qc,
                                  Evaluation.evaluation_output_dir, Evaluation.control_output_dir, Collected.collected)
    if (params.OUTPUT_IGVF_PACKAGE) {
        check_igvf_dashboard(Dashboard.archive)
    }
    emit:
    inference_mudata = Inference.inference_mudata
    matching_report = Balanced.report
}
