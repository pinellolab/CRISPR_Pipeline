
process mergedResults {
    cache 'lenient'

    input:
    path sceptre_per_guide
    path sceptre_per_element
    path perturbo_per_guide
    path perturbo_per_element
    path base_mudata
    // Whether this call needs the intermediate MuData. The default workflow does not
    // (mergeMudata builds the published one); the sceptre,perturbo workflow ends here,
    // so it does.
    val write_mudata

    output:
    // Optional: mergeMudata assembles the published MuData from these tables, so the
    // intermediate copy is skipped unless something asks for it.
    path "inference_mudata.h5mu", optional: true, emit: inference_mudata
    path "per_element_output.tsv.gz", emit: per_element_output
    path "per_guide_output.tsv.gz", emit: per_guide_output

    script:
        """
        merge_method_results.py \\
            --sceptre_per_guide ${sceptre_per_guide} \\
            --sceptre_per_element ${sceptre_per_element} \\
            --perturbo_per_guide ${perturbo_per_guide} \\
            --perturbo_per_element ${perturbo_per_element} \\
            --base_mudata ${base_mudata} \\
            ${write_mudata ? '--write_mudata' : ''}
        """

}
