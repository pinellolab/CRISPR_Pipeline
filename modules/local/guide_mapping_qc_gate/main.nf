process guide_mapping_qc_gate {
    tag "${params.GUIDE_MAPPING_QC_action}"
    // Stopping here makes this MuData the run's deliverable; in a full run
    // mergeMudata publishes the final object and this one is an intermediate.
    publishDir "${params.outdir}/guide_mapping_qc", mode: params.publish_dir_mode, pattern: "guide_mapping_qc_passed_mudata.h5mu", enabled: params.stop_after_guide_mapping_qc, overwrite: true
    cache 'lenient'

    input:
        path mudata_input
        path guide_mapping_qc_dir

    output:
        path "guide_mapping_qc_passed_mudata.h5mu", emit: mudata

    script:
        """
        guide_mapping_qc_gate.py \
            ${guide_mapping_qc_dir}/guide_mapping_qc.json \
            --action ${params.GUIDE_MAPPING_QC_action}
        # A symlink is copied to object storage as the link target itself,
        # so the next process stages a dangling link. Copy the file.
        cp -L ${mudata_input} guide_mapping_qc_passed_mudata.h5mu
        """
}
