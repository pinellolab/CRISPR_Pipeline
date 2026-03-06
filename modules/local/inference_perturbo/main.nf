
process inference_perturbo {
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
    path pairs_to_test_file
    val run_all_pairs
    val write_mudata
    
    output:
    path "inference_mudata.h5mu", emit: inference_mudata, optional: true
    path "perturbo_done.token", emit: completion_token
    path "perturbo_cis_per_element_output.tsv.gz", emit: per_element_output
    path "perturbo_cis_per_guide_output.tsv.gz", emit: per_guide_output

    script:
        def testAllPairsArg = run_all_pairs ? "--test_all_pairs" : ""
        def pairsArg = (!run_all_pairs && pairs_to_test_file.getName() != 'no_pairs_to_test.tsv') ? "--pairs_to_test_file ${pairs_to_test_file}" : ""
        """
        # Run PerTurbo inference for per-element results
        perturbo_inference.py ${mudata} perturbo_cis_per_element_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type element ${testAllPairsArg} ${pairsArg}
        
        # Run PerTurbo inference for per-guide results  
        perturbo_inference.py ${mudata} perturbo_cis_per_guide_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type guide ${testAllPairsArg} ${pairsArg}
        
        if [ "${write_mudata}" = "true" ]; then
            add_perturbo_results_to_mudata.py \\
                --per_guide_tsv perturbo_cis_per_guide_output.tsv.gz \\
                --per_element_tsv perturbo_cis_per_element_output.tsv.gz \\
                --base_mudata ${mudata} \\
                --output inference_mudata.h5mu
        fi
        touch perturbo_done.token
        """
}
