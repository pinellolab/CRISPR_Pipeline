
process inference_perturbo {
    cache 'lenient'
    publishDir './pipeline_outputs'
    maxForks 1

    input:
    path mudata
    val inference_method
    val efficiency_mode
    val mode // 'cis' or 'trans'
    
    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "perturbo_${mode}_per_element_output.tsv.gz", emit: per_element_output
    path "perturbo_${mode}_per_guide_output.tsv.gz", emit: per_guide_output

    script:
        def test_all_pairs_flag = mode == 'trans' ? '--test_all_pairs' : ''
        """
        # Run PerTurbo inference for per-element results
        perturbo_inference.py ${mudata} perturbo_${mode}_per_element_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type element ${test_all_pairs_flag}
        
        # Run PerTurbo inference for per-guide results  
        perturbo_inference.py ${mudata} perturbo_${mode}_per_guide_output.tsv.gz --efficiency_mode ${efficiency_mode} --inference_type guide ${test_all_pairs_flag}
        
        # Add both results to the base mudata file
        add_perturbo_results_to_mudata.py \\
            --per_guide_tsv perturbo_${mode}_per_guide_output.tsv.gz \\
            --per_element_tsv perturbo_${mode}_per_element_output.tsv.gz \\
            --base_mudata ${mudata} \\
            --output inference_mudata.h5mu
        """
}
