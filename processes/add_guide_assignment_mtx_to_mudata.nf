process add_guide_assignment_mtx_to_mudata {
    cache 'lenient'
    debug true
  
    input:
    tuple val(key), path(mudata_input), path(guide_assignment_mtx)

    output:
    path "${mudata_input.simpleName}_output.h5mu", emit: guide_assignment_mudata_output
    
    script:
      """
        echo "Processing ${key}: add ${guide_assignment_mtx} to ${mudata_input}"
        export NUMBA_CACHE_DIR=/tmp
        add_guide_assignment.py --guide_assignment ${guide_assignment_mtx} --mudata ${mudata_input}
        mv sceptre_assignment_mudata.h5mu ${mudata_input.simpleName}_output.h5mu
      """
}