process additional_qc_plots {
    input:
    path inference_mudata

    output:
    path 'additional_qc', emit: additional_qc

    script:
    """
    mkdir -p additional_qc
    touch additional_qc/placeholder.txt
    """
}