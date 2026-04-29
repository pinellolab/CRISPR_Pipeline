process tf_benchmark {
    input:
    path inference_mudata
    path gencode_gtf
    path encode_bed_dir

    output:
    path 'benchmark_output', emit: benchmark_output

    script:
    """
    mkdir -p benchmark_output
    touch benchmark_output/placeholder.txt
    """
}