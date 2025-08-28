
process evaluation_plot_default {

    cache 'lenient'
    input:

    path mdata
    val user_central_nodes
    path gencode_gtf

    output:
    path "evaluation_output" , emit: evaluation_output

    script:
            """
            #network_plot.py ${mdata} --central_nodes ${user_central_nodes} --default --min_weight 0.1
            #volcano_plot.py ${mdata} --log2_fc 1 --p_value 0.05 --default
            igv.py ${mdata} --gtf ${gencode_gtf} --default
            """
}
