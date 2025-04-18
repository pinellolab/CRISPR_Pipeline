process guide_assignment_sceptre {
    cache 'lenient'
    debug true
  
    input:
    each mudata_input

    output:
    path "guide_assignment.mtx", emit: guide_assignment_mtx_output
    
    script:
    """
      export NUMBA_CACHE_DIR=/tmp
      assign_grnas_sceptre.R ${mudata_input}
    """
}