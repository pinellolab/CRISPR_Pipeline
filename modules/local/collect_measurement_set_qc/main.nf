process collect_measurement_set_qc {
    cache 'lenient'

    input:
    path measurement_set_qc_dirs

    output:
    path "figures", emit: figures_dir

    script:
    def inputs = measurement_set_qc_dirs.collect { it.toString() }.sort().join(' ')
    """
    collect_measurement_set_qc.py ${inputs} --output figures
    """
}
