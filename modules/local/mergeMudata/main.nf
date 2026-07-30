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
        path "cis_per_guide_output.*", emit: cis_per_guide_output
        path "cis_per_element_output.*", emit: cis_per_element_output
        path "trans_per_guide_output.*", emit: trans_per_guide_output
        path "trans_per_element_output.*", emit: trans_per_element_output
        path "catalog_per_element_output.*", emit: catalog_per_element_output

    script:
    def results_ext = params.INFERENCE_PERTURBO_TRANS_RESULTS_FORMAT == 'parquet' ? 'parquet' : 'tsv.gz'
    """
        merge_cis_trans_results.py \\
            --cis_per_guide ${cis_per_guide} \\
            --cis_per_element ${cis_per_element} \\
            --trans_per_guide ${trans_per_guide} \\
            --trans_per_element ${trans_per_element} \\
            --base_mudata ${base_mudata} \\
            --output inference_mudata.h5mu \\
            --results_format ${params.INFERENCE_PERTURBO_TRANS_RESULTS_FORMAT}

        build_catalog_per_element_output.py \\
            --cis_per_element cis_per_element_output.${results_ext} \\
            --trans_per_element trans_per_element_output.${results_ext} \\
            --mudata inference_mudata.h5mu \\
            --output catalog_per_element_output.${results_ext}
    """
}
