
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
    path dependency_token // dependency to force this to run after cis analysis
    val write_mudata

    output:
    path "inference_mudata.h5mu", emit: inference_mudata, optional: true
    path "perturbo_trans_done.token", emit: completion_token
    path "perturbo_trans_per_element_output.tsv.gz", emit: per_element_output
    path "perturbo_trans_per_guide_output.tsv.gz", emit: per_guide_output

    script:
        """
        # Run PerTurbo inference for per-element results
        perturbo_inference.py ${mudata} perturbo_trans_per_element_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type element --test_all_pairs
        
        # Run PerTurbo inference for per-guide results  
        perturbo_inference.py ${mudata} perturbo_trans_per_guide_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type guide --test_all_pairs
        
        if [ "${write_mudata}" = "true" ]; then
            add_perturbo_results_to_mudata.py \\
                --per_guide_tsv perturbo_trans_per_guide_output.tsv.gz \\
                --per_element_tsv perturbo_trans_per_element_output.tsv.gz \\
                --base_mudata ${mudata} \\
                --output inference_mudata.h5mu
        fi
        touch perturbo_trans_done.token
        """
}
