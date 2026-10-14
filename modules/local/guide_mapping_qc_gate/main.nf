process guide_mapping_qc_gate {
    tag "${params.GUIDE_MAPPING_QC_action}"
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
