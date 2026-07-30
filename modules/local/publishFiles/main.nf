process publishFiles {
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, mode: 'copy', overwrite: true

    input:
        path cis_per_element_results
        path cis_per_guide_results
        path trans_per_element_results
        path trans_per_guide_results

    output:
        path "cis_per_element_output.*"
        path "cis_per_guide_output.*"
        path "trans_per_element_output.*"
        path "trans_per_guide_output.*"

    script:
    def results_ext = params.INFERENCE_PERTURBO_TRANS_RESULTS_FORMAT == 'parquet' ? 'parquet' : 'tsv.gz'
    """
        # Check all files exist
        for file in "${cis_per_element_results}" "${cis_per_guide_results}" "${trans_per_element_results}" "${trans_per_guide_results}"; do
            if [[ ! -f "\$file" ]]; then
                echo "ERROR: File not found: \$file"
                exit 1
            fi
            echo "Found: \$file"
        done

        # Copy to create actual files from symlinks
        cp "${cis_per_element_results}" cis_per_element_output.${results_ext}
        cp "${cis_per_guide_results}" cis_per_guide_output.${results_ext}
        cp "${trans_per_element_results}" trans_per_element_output.${results_ext}
        cp "${trans_per_guide_results}" trans_per_guide_output.${results_ext}

        echo "All files copied and ready for publishing"
    """
}
