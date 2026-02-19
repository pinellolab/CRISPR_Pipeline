process tf_benchmark {
    cache 'lenient'

    input:
        path mudata
        path gencode_gtf
        path encode_bed_dir

    output:
        path "benchmark_output", emit: benchmark_output

    script:
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        tf_benchmark.py --mudata ${mudata} --gtf ${gencode_gtf} --encode_bed_dir ${encode_bed_dir} --output_dir benchmark_output
        """
}
