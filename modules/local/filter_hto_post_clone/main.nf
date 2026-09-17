process filter_hto_post_clone {
    tag "min_${params.HTO_min_positive_cells}_singlets_${params.HTO_keep_singlets_only}"
    cache 'lenient'
    publishDir "${params.outdir}/hashing_qc", mode: params.publish_dir_mode, pattern: "hto_qc/*", overwrite: true

    input:
        path mudata_input

    output:
        path "hto_filtered_mudata.h5mu", emit: filtered_mudata
        path "post_clone_hashing_filtered.h5ad", emit: filtered_hashing
        path "post_clone_hashing_unfiltered.h5ad", emit: unfiltered_hashing
        path "hto_qc", emit: hto_qc

    script:
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        filter_hto_post_clone.py \
            ${mudata_input} \
            hto_filtered_mudata.h5mu \
            --outdir hto_qc \
            --min-positive-cells ${params.HTO_min_positive_cells} \
            --singlet-only ${params.HTO_keep_singlets_only} \
            --batch-column ${params.QC_batch_col ?: 'batch'} \
            --filtered-hashing-output post_clone_hashing_filtered.h5ad \
            --unfiltered-hashing-output post_clone_hashing_unfiltered.h5ad
        """

    stub:
        """
        mkdir -p hto_qc
        touch hto_filtered_mudata.h5mu post_clone_hashing_filtered.h5ad post_clone_hashing_unfiltered.h5ad \
            hto_qc/hto_filter_flow.tsv hto_qc/hto_filter_steps_stub.png
        """
}
