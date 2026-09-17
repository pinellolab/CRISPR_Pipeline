process filter_guide_assignment_qc {
    tag "max_assigned_guides_${params.GUIDE_ASSIGNMENT_max_guides_per_cell}"
    cache 'lenient'
    publishDir "${params.outdir}/guide_assignment_qc", mode: params.publish_dir_mode, pattern: "guide_assignment_qc/*", overwrite: true

    input:
        path mudata_input

    output:
        path "guide_assignment_filtered_mudata.h5mu", emit: filtered_mudata
        path "guide_assignment_qc", emit: guide_assignment_qc

    script:
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        filter_guide_assignment_qc.py \
            ${mudata_input} \
            guide_assignment_filtered_mudata.h5mu \
            --outdir guide_assignment_qc \
            --max-guides-per-cell ${params.GUIDE_ASSIGNMENT_max_guides_per_cell} \
            --batch-column ${params.QC_batch_col ?: 'batch'}
        """

    stub:
        """
        mkdir -p guide_assignment_qc
        touch guide_assignment_filtered_mudata.h5mu \
            guide_assignment_qc/guide_assignment_filter_flow.tsv \
            guide_assignment_qc/guide_assignment_filter_steps_all.png
        """
}
