process guide_assignment_sceptre {
    cache 'lenient'

    input:
    path mudata_input
    val probability_threshold
    val SCEPTRE_n_em_rep

    output:
    path "${mudata_input.simpleName}_output.h5mu", emit: guide_assignment_mudata_output

    script:
    """
    # SCEPTRE's mixture-model assignment does real dense work in OpenBLAS: measured
    # at 16 of 16 allotted cores on CPU with this cap, 64 without it. The global
    # env {} cap of 1 is right for the forked workers in inference_sceptre, not here.
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus}
    export NUMBA_CACHE_DIR=/tmp
    assign_grnas_sceptre.R ${mudata_input} ${probability_threshold} ${SCEPTRE_n_em_rep}
    add_guide_assignment.py --guide_assignment guide_assignment.mtx --mudata ${mudata_input}
    mv sceptre_assignment_mudata.h5mu ${mudata_input.simpleName}_output.h5mu
    """
}
