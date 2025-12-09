nextflow.enable.dsl=2

include { evaluation_plot } from '../../../modules/local/evaluation_plot'
include { evaluation_plot_default } from '../../../modules/local/evaluation_plot_default'
include { evaluation_undefined_plot } from '../../../modules/local/evaluation_undefined_plot'
include { evaluation_undefined_plot_default } from '../../../modules/local/evaluation_undefined_plot_default'
include { evaluation_controls } from '../../../modules/local/evaluation_controls'

workflow evaluation_pipeline {

    take:
    gencode_gtf
    inference_mudata

    main:
    if (params.INFERENCE_method == 'default') {
        if (params.NETWORK_custom_central_nodes == 'undefined') {
            Evaluation = evaluation_undefined_plot_default(inference_mudata, gencode_gtf, params.NETWORK_central_nodes_num)
        } else {
            Evaluation = evaluation_plot_default(inference_mudata, params.NETWORK_custom_central_nodes, gencode_gtf)
        }
    } else {
        if (params.NETWORK_custom_central_nodes == 'undefined') {
            Evaluation = evaluation_undefined_plot(inference_mudata, gencode_gtf, params.NETWORK_central_nodes_num)
        } else {
            Evaluation = evaluation_plot(inference_mudata, params.NETWORK_custom_central_nodes, gencode_gtf)
        }
    }
    evaluation_controls(inference_mudata)


    emit:
    evaluation_output_dir = Evaluation.evaluation_output
}
