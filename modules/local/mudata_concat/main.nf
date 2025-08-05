process mudata_concat {
    cache 'lenient'
    debug true

    input:
        path (mudata_input, name: "?/*")
        val fraction_cells
    output:
        path "concat_mudata.h5mu", emit: concat_mudata

    script:
    """
        mudata_concat.py -i ${mudata_input} -o concat_mudata.h5mu -g ${fraction_cells}
    """
}
