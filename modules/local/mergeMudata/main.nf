process mergeMudata {
    cache 'lenient'
    debug true
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, mode: 'copy', overwrite: true

    input:
        path cis_per_guide
        path cis_per_element
        path trans_per_guide
        path trans_per_element
        path base_mudata

    output:
        path "inference_mudata.h5mu", emit: inference_mudata
        path "cis_per_guide_output.tsv.gz", emit: cis_per_guide_output
        path "cis_per_element_output.tsv.gz", emit: cis_per_element_output
        path "trans_per_guide_output.tsv.gz", emit: trans_per_guide_output
        path "trans_per_element_output.tsv.gz", emit: trans_per_element_output

    script:
    """
        merge_cis_trans_results.py \\
            --cis_per_guide ${cis_per_guide} \\
            --cis_per_element ${cis_per_element} \\
            --trans_per_guide ${trans_per_guide} \\
            --trans_per_element ${trans_per_element} \\
            --base_mudata ${base_mudata} \\
            --output inference_mudata.h5mu
    """
}
