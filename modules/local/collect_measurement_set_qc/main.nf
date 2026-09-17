process collect_measurement_set_qc {
    cache 'lenient'

    // Publish these artifacts as soon as the aggregation task completes so the
    // live W&B monitor can display per-measurement-set knees and QC distributions
    // before the final dashboard process runs.
    publishDir "${params.outdir}/measurement_set_qc", mode: params.publish_dir_mode, overwrite: true

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
