process guide_assignment_sceptre {
    cache 'lenient'
    debug true
  
    input:
    each path(mudata_input)

    output:
    path "${mudata_input.simpleName}.mtx", emit: guide_assignment_mtx_output
    
    script:
    """
      export NUMBA_CACHE_DIR=/tmp
      assign_grnas_sceptre.R ${mudata_input}
      mv guide_assignment.mtx ${mudata_input.simpleName}.mtx
    """
}