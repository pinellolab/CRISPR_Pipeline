process mergeMudata {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, mode: 'copy', overwrite: true

    input:
        path local_analysis_per_guide
        path local_analysis_per_element
        path global_analysis_per_guide
        path global_analysis_per_element
        path base_mudata

    output:
        path "inference_mudata.h5mu", emit: inference_mudata
        path "local_analysis_per_guide_output.*", emit: local_analysis_per_guide_output
        path "local_analysis_per_element_output.*", emit: local_analysis_per_element_output
        path "global_analysis_per_guide_output.*", emit: global_analysis_per_guide_output
        path "global_analysis_per_element_output.*", emit: global_analysis_per_element_output
        path "catalog_per_element_output.*", emit: catalog_per_element_output
        path "catalog_per_guide_output.*", emit: catalog_per_guide_output

    script:
    def results_ext = params.INFERENCE_PERTURBO_TRANS_RESULTS_FORMAT == 'parquet' ? 'parquet' : 'tsv.gz'
    """
        merge_local_global_results.py \\
            --local_analysis_per_guide ${local_analysis_per_guide} \\
            --local_analysis_per_element ${local_analysis_per_element} \\
            --global_analysis_per_guide ${global_analysis_per_guide} \\
            --global_analysis_per_element ${global_analysis_per_element} \\
            --base_mudata ${base_mudata} \\
            --output inference_mudata.h5mu \\
            --results_format ${params.INFERENCE_PERTURBO_TRANS_RESULTS_FORMAT}

        build_catalog_per_element_output.py \\
            --local_analysis_per_element local_analysis_per_element_output.${results_ext} \\
            --global_analysis_per_element global_analysis_per_element_output.${results_ext} \\
            --mudata inference_mudata.h5mu \\
            --output catalog_per_element_output.${results_ext}

        build_catalog_per_guide_output.py \\
            --local_analysis_per_guide local_analysis_per_guide_output.${results_ext} \\
            --global_analysis_per_guide global_analysis_per_guide_output.${results_ext} \\
            --mudata inference_mudata.h5mu \\
            --output catalog_per_guide_output.${results_ext}
    """
}
