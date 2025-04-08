process mudata_concat {
    cache 'lenient'
    debug true

    input:
        path mudata_input
        
    output:
        path "concat_mudata.h5mu", emit: concat_mudata

    script:
    """
        mudata_concat.py -i "*_mudata_output.h5mu" -o concat_mudata.h5mu
    """
}