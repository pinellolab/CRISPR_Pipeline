process remove_clonal_cells {
    tag "${params.CLONE_REMOVAL_action}"
    cache 'lenient'
    publishDir "${params.outdir}/clone_removal", mode: params.publish_dir_mode, pattern: "clone_qc/*", overwrite: true

    input:
        path mudata_input

    output:
        path "clone_filtered_mudata.h5mu", emit: filtered_mudata
        path "clone_qc", emit: clone_qc

    script:
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        remove_clonal_cells.py \
            ${mudata_input} \
            clone_filtered_mudata.h5mu \
            --outdir clone_qc \
            --alpha ${params.CLONE_REMOVAL_alpha} \
            --min-clone-size ${params.CLONE_REMOVAL_min_clone_size} \
            --action ${params.CLONE_REMOVAL_action}
        """

    stub:
        """
        mkdir -p clone_qc
        touch clone_filtered_mudata.h5mu clone_qc/clone_metrics.tsv clone_qc/clone_filter_summary.png
        """
}
