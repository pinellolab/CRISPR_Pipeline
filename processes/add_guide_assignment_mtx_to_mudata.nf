process add_guide_assignment_mtx_to_mudata {
    cache 'lenient'
    debug true
  
    input:
    each mudata_input
    each guide_assignment_mtx

    output:
    path "${mudata_input.simpleName}_output.h5mu", emit: guide_assignment_mudata_output
    
    script:
    """
      export NUMBA_CACHE_DIR=/tmp
      add_guide_assignment.py --guide_assignment ${guide_assignment_mtx} --mudata ${mudata_input}
      mv sceptre_assignment_mudata.h5mu ${mudata_input.simpleName}_output.h5mu
    """
}