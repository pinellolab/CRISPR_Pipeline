process sequencing_saturation {
    tag "RNA BUS rarefaction"
    cache 'lenient'
    // Publish the emitted directory itself. A nested `saturation_qc/*` pattern
    // does not match a directory output, leaving the published folder empty.
    publishDir "${params.outdir}", mode: params.publish_dir_mode, overwrite: true

    input:
        path mapping_dirs
        path filtered_anndata
        path transcriptome_t2g
        path parsed_covariates

    output:
        path "sequencing_saturation", emit: saturation_qc

    script:
        def mapping_args = mapping_dirs.collect { it.toString() }.sort().join(' ')
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        sequencing_saturation.py \
            --mapping-dirs ${mapping_args} \
            --filtered-anndata ${filtered_anndata} \
            --t2g ${transcriptome_t2g} \
            --covariates ${parsed_covariates} \
            --fractions '${params.SATURATION_downsample_fractions}' \
            --outdir sequencing_saturation
        """

    stub:
        """
        mkdir -p sequencing_saturation
        touch sequencing_saturation/sequencing_saturation_curve.tsv \
              sequencing_saturation/sequencing_saturation_metrics.tsv \
              sequencing_saturation/sequencing_saturation_curve.png
        """
}
