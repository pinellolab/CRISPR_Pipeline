process publishFiles {
    publishDir './pipeline_outputs', mode: 'copy', overwrite: true

    input:
        path cis_element
        path cis_guide
        path trans_element
        path trans_guide

    output:
        path "cis_per_element_output.tsv"
        path "cis_per_guide_output.tsv"
        path "trans_per_element_output.tsv"
        path "trans_per_guide_output.tsv"

    script:
    """
        # Check all files exist
        for file in "${cis_element}" "${cis_guide}" "${trans_element}" "${trans_guide}"; do
            if [[ ! -f "\$file" ]]; then
                echo "ERROR: File not found: \$file"
                exit 1
            fi
            echo "Found: \$file"
        done

        # Copy  to create actual files from symlinks
        cp "${cis_element}" cis_per_element_output.tsv
        cp "${cis_guide}" cis_per_guide_output.tsv
        cp "${trans_element}" trans_per_element_output.tsv
        cp "${trans_guide}" trans_per_guide_output.tsv

        # Remove original files from pipeline_outputs
        rm -f ./pipeline_outputs/per_element_output.tsv
        rm -f ./pipeline_outputs/per_guide_output.tsv

        echo "All files copied and ready for publishing"
    """
}
