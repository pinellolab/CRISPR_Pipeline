nextflow.enable.dsl=2

include { evaluation_plot } from '../../../modules/local/evaluation_plot'
//include { evaluation_plot_default } from '../../../modules/local/evaluation_plot_default'
//include { evaluation_undefined_plot } from '../../../modules/local/evaluation_undefined_plot'
//include { evaluation_undefined_plot_default } from '../../../modules/local/evaluation_undefined_plot_default'
include { evaluation_controls } from '../../../modules/local/evaluation_controls'

workflow evaluation_pipeline {

    take:
    gencode_gtf
    inference_mudata

    main:

    Evaluation = evaluation_plot(inference_mudata, gencode_gtf)

    Controls_evaluation = evaluation_controls(inference_mudata)


    emit:
    evaluation_output_dir = Evaluation.evaluation_output
    control_output_dir = Controls_evaluation.evaluation_controls
}
