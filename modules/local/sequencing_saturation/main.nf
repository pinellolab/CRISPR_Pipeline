process sequencing_saturation {
    tag "RNA BUS rarefaction"
    cache 'lenient'
    publishDir "${params.outdir}/sequencing_saturation", mode: params.publish_dir_mode, pattern: "saturation_qc/*", overwrite: true

    input:
        path mapping_dirs
        path filtered_anndata
        path transcriptome_t2g
        path parsed_covariates

    output:
        path "saturation_qc", emit: saturation_qc

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
            --outdir saturation_qc
        """

    stub:
        """
        mkdir -p saturation_qc
        touch saturation_qc/sequencing_saturation_curve.tsv \
              saturation_qc/sequencing_saturation_metrics.tsv \
              saturation_qc/sequencing_saturation_curve.png
        """
}
