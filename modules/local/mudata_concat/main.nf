process mudata_concat {
    cache 'lenient'
    debug true

    input:
        path (mudata_input, name: "?/*")
        val fraction_cells
        val dual_guide
    output:
        path "concat_mudata.h5mu", emit: concat_mudata

    script:
    """
        mudata_concat.py -i ${mudata_input} -o concat_mudata.h5mu -g ${fraction_cells}

        if [ "${dual_guide}" = "true" ]; then
            echo "Dual guide mode enabled, processing accordingly."
            collapse_guides.py concat_mudata.h5mu  concat_mudata_collapsed.h5mu
            mv concat_mudata_collapsed.h5mu concat_mudata.h5mu

        else
            echo "Single guide mode, proceeding with standard processing."
        fi
    """
}








