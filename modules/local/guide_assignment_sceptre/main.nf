process guide_assignment_sceptre {
    cache 'lenient'
    debug true

    input:
    path mudata_input
    val probability_threshold
    val SCEPTRE_n_em_rep

    output:
    path "${mudata_input.simpleName}_out.h5mu", emit: guide_assignment_mudata_output

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    assign_grnas_sceptre.R ${mudata_input} ${probability_threshold} ${SCEPTRE_n_em_rep} ${task.cpus}
    """
}
