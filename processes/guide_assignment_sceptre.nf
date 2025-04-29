process guide_assignment_sceptre {
    cache 'lenient'
    debug true

    input:
    path mudata_input

    output:
    path "${mudata_input.simpleName}_output.h5mu", emit: guide_assignment_mudata_output

    script:
    """
    export NUMBA_CACHE_DIR=/tmp
    assign_grnas_sceptre.R ${mudata_input}
    add_guide_assignment.py --guide_assignment guide_assignment.mtx --mudata ${mudata_input}
    mv sceptre_assignment_mudata.h5mu ${mudata_input.simpleName}_output.h5mu
    """
}
