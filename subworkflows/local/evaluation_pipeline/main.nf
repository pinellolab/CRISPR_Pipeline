process EVALUATION_STUB {
    input:
    path gencode_gtf
    path inference_mudata

    output:
    path 'evaluation_output', emit: evaluation_output_dir
    path 'plots',             emit: control_output_dir

    script:
    """
    mkdir -p evaluation_output plots
    touch evaluation_output/placeholder.txt
    touch plots/placeholder.txt
    """
}

workflow evaluation_pipeline {
    take:
    gencode_gtf
    inference_mudata

    main:
    EVALUATION_STUB(gencode_gtf, inference_mudata)

    emit:
    evaluation_output_dir = EVALUATION_STUB.out.evaluation_output_dir
    control_output_dir    = EVALUATION_STUB.out.control_output_dir
}