
process mergedResults {
    cache 'lenient'

    input:
    path sceptre_per_guide
    path sceptre_per_element
    path perturbo_per_guide
    path perturbo_per_element
    path base_mudata

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "per_element_output.tsv.gz", emit: per_element_output
    path "per_guide_output.tsv.gz", emit: per_guide_output

    script:
        """
        merge_method_results.py \\
            --sceptre_per_guide ${sceptre_per_guide} \\
            --sceptre_per_element ${sceptre_per_element} \\
            --perturbo_per_guide ${perturbo_per_guide} \\
            --perturbo_per_element ${perturbo_per_element} \\
            --base_mudata ${base_mudata}
        """

}
