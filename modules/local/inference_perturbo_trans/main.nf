
process inference_perturbo_trans {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }

    input:
    path mudata
    val inference_method
    val efficiency_mode
    val dummy // fake dependency to force this to run after cis analysis

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "perturbo_trans_per_element_output.tsv.gz", emit: per_element_output
    path "perturbo_trans_per_guide_output.tsv.gz", emit: per_guide_output

    script:
        """
        # Run PerTurbo inference for per-element results
        perturbo_inference.py ${mudata} perturbo_trans_per_element_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type element --test_all_pairs
        
        # Run PerTurbo inference for per-guide results  
        perturbo_inference.py ${mudata} perturbo_trans_per_guide_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type guide --test_all_pairs
        
        # Add both results to the base mudata file
        add_perturbo_results_to_mudata.py \\
            --per_guide_tsv perturbo_trans_per_guide_output.tsv.gz \\
            --per_element_tsv perturbo_trans_per_element_output.tsv.gz \\
            --base_mudata ${mudata} \\
            --output inference_mudata.h5mu
        """
}
