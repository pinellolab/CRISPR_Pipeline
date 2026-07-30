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
        path local_analysis_per_element_results
        path local_analysis_per_guide_results
        path global_analysis_per_element_results
        path global_analysis_per_guide_results

    output:
        path "local_analysis_per_element_output.tsv.gz"
        path "local_analysis_per_guide_output.tsv.gz"
        path "global_analysis_per_element_output.tsv.gz"
        path "global_analysis_per_guide_output.tsv.gz"

    script:
    """
        # Check all files exist
        for file in "${local_analysis_per_element_results}" "${local_analysis_per_guide_results}" "${global_analysis_per_element_results}" "${global_analysis_per_guide_results}"; do
            if [[ ! -f "\$file" ]]; then
                echo "ERROR: File not found: \$file"
                exit 1
            fi
            echo "Found: \$file"
        done

        # Copy to create actual files from symlinks
        cp "${local_analysis_per_element_results}" local_analysis_per_element_output.tsv.gz
        cp "${local_analysis_per_guide_results}" local_analysis_per_guide_output.tsv.gz
        cp "${global_analysis_per_element_results}" global_analysis_per_element_output.tsv.gz
        cp "${global_analysis_per_guide_results}" global_analysis_per_guide_output.tsv.gz

        echo "All files copied and ready for publishing"
    """
}
