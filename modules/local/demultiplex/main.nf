
process demultiplex{

    cache 'lenient'

    input:
    path adata_path

    output:
    path "*_hashing_filtered_demux.h5ad", emit: hashing_demux_anndata
    path "*_hashing_unfiltered_demux.h5ad", emit: hashing_demux_unfiltered_anndata
    path "gmm_demux_qc.json", emit: gmm_demux_qc


    script:
    def rejectDegenerate = params.HTO_GMM_reject_nonzero_positive ? '--reject-nonzero-positive' : ''
    """
    adata_name=\$(basename ${adata_path} .h5ad)
    hto_string=\$(demultiplex_prepare.py --adata ${adata_path} -o demuxfile)

    (
        export OPENBLAS_NUM_THREADS=1
        GMM_DEMUX_PATH=\$(type -p GMM-demux)
        export PATH=\$(dirname \$GMM_DEMUX_PATH):\$PATH
        run_gmm_demux_with_qc.py \
            --matrix-dir demuxfile \
            --hto-names \$hto_string \
            --output-dir FULL \
            --ssd-output-dir SSD_mtx \
            --qc-json gmm_demux_qc.json \
            --seed ${params.HTO_GMM_random_seed} \
            --max-attempts ${params.HTO_GMM_max_seed_attempts} \
            ${rejectDegenerate}
    )

    demultiplex_filter.py --adata ${adata_path} --demux_report FULL/GMM_full.csv --demux_config FULL/GMM_full.config --demux_qc gmm_demux_qc.json --filtered_output \${adata_name}_hashing_filtered_demux.h5ad --unfiltered_output \${adata_name}_hashing_unfiltered_demux.h5ad
    """
}
