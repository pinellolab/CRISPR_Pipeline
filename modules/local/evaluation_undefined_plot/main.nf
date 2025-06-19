
process evaluation_undefined_plot {

    cache 'lenient'

    input:

    path mdata
    path gencode_gtf
    val num_nodes

    output:
    path "evaluation_output" , emit: evaluation_output

    script:
            """
            network_plot_undefined.py ${mdata} --num_nodes ${num_nodes} --min_weight 0.1
            volcano_plot.py ${mdata} --log2_fc 1 --p_value 0.05
            igv.py ${mdata} --gtf ${gencode_gtf}
            """
}
