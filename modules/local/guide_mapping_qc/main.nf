process guide_mapping_qc {
    tag "orientation_${reverse_complement_guides}"
    cache 'lenient'
    publishDir "${params.outdir}/guide_mapping_qc", mode: params.publish_dir_mode, pattern: "guide_mapping_qc", overwrite: true

    input:
        path adata_rna
        path adata_guide
        path mudata_input
        path guide_metadata
        val reverse_complement_guides
        val spacer_tag

    output:
        path "guide_mapping_qc", emit: qc_dir

    script:
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}
        python ${projectDir}/bin/guide_mapping_qc.py \
            --rna ${adata_rna} \
            --guide ${adata_guide} \
            --mudata ${mudata_input} \
            --guide-metadata ${guide_metadata} \
            --reverse-complement-guides ${reverse_complement_guides} \
            --spacer-tag "${spacer_tag ?: ''}" \
            --batch-column "${params.QC_batch_col ?: 'batch'}" \
            --min-overlap-cells-per-set ${params.GUIDE_MAPPING_QC_min_overlap_cells_per_set} \
            --min-guide-to-rna-fraction ${params.GUIDE_MAPPING_QC_min_guide_to_rna_fraction} \
            --min-overlap-to-guide-fraction ${params.GUIDE_MAPPING_QC_min_overlap_to_guide_fraction} \
            --min-recovered-guide-fraction ${params.GUIDE_MAPPING_QC_min_recovered_guide_fraction} \
            --outdir guide_mapping_qc
        """
}
