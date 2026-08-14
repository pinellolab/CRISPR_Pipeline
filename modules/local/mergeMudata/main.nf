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

    script:
    """
        export POLARS_MAX_THREADS=${task.cpus}

        merge_local_global_results.py \\
            --local_analysis_per_guide ${local_analysis_per_guide} \\
            --local_analysis_per_element ${local_analysis_per_element} \\
            --global_analysis_per_guide ${global_analysis_per_guide} \\
            --global_analysis_per_element ${global_analysis_per_element} \\
            --base_mudata ${base_mudata} \\
            --output inference_mudata.h5mu \\
            --results_format ${params.INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT}
    """
}
