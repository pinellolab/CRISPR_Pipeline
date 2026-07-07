process tf_benchmark {
    cache 'lenient'

    input:
        path mudata
        path gencode_gtf
        path encode_bed_dir
        val demo_mode

    output:
        path "benchmark_output", emit: benchmark_output

    script:
        def demo_arg = demo_mode ? '--demo-mode' : ''
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        tf_benchmark.py --mudata ${mudata} --gtf ${gencode_gtf} --encode_bed_dir ${encode_bed_dir} --output_dir benchmark_output ${demo_arg}
        """
}
