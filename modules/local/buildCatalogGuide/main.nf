process buildCatalogGuide {
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
        path global_analysis_per_guide
        path inference_mudata

    output:
        path "catalog_per_guide_output.*", emit: catalog_per_guide_output

    script:
    def results_ext = params.INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT == 'parquet' ? 'parquet' : 'tsv.gz'
    """
        export POLARS_MAX_THREADS=${task.cpus}

        build_catalog_per_guide_output.py \\
            --local_analysis_per_guide ${local_analysis_per_guide} \\
            --global_analysis_per_guide ${global_analysis_per_guide} \\
            --mudata ${inference_mudata} \\
            --output catalog_per_guide_output.${results_ext}
    """
}
