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
        path local_analysis_per_guide
        path local_analysis_per_element
        path global_analysis_per_guide
        path global_analysis_per_element
        path base_mudata

    output:
        path "inference_mudata.h5mu", emit: inference_mudata
        path "local_analysis_per_guide_output.tsv.gz", emit: local_analysis_per_guide_output
        path "local_analysis_per_element_output.tsv.gz", emit: local_analysis_per_element_output
        path "global_analysis_per_guide_output.tsv.gz", emit: global_analysis_per_guide_output
        path "global_analysis_per_element_output.tsv.gz", emit: global_analysis_per_element_output
        path "catalog_per_element_output.tsv.gz", emit: catalog_per_element_output
        path "catalog_per_guide_output.tsv.gz", emit: catalog_per_guide_output

    script:
    """
        merge_local_global_results.py \\
            --local_analysis_per_guide ${local_analysis_per_guide} \\
            --local_analysis_per_element ${local_analysis_per_element} \\
            --global_analysis_per_guide ${global_analysis_per_guide} \\
            --global_analysis_per_element ${global_analysis_per_element} \\
            --base_mudata ${base_mudata} \\
            --output inference_mudata.h5mu

        build_catalog_per_element_output.py \\
            --local_analysis_per_element local_analysis_per_element_output.tsv.gz \\
            --global_analysis_per_element global_analysis_per_element_output.tsv.gz \\
            --mudata inference_mudata.h5mu \\
            --output catalog_per_element_output.tsv.gz

        build_catalog_per_guide_output.py \\
            --local_analysis_per_guide local_analysis_per_guide_output.tsv.gz \\
            --global_analysis_per_guide global_analysis_per_guide_output.tsv.gz \\
            --mudata inference_mudata.h5mu \\
            --output catalog_per_guide_output.tsv.gz
    """
}
