
process demultiplex{

    cache 'lenient'
    debug true

    input:
    path adata_path

    output:
    path "*_hashing_filtered_demux.h5ad", emit: hashing_demux_anndata
    path "*_hashing_unfiltered_demux.h5ad", emit: hashing_demux_unfiltered_anndata


    script:
    """
    adata_name=\$(basename ${adata_path} .h5ad)
    hto_string=\$(demultiplex_prepare.py --adata ${adata_path} -o demuxfile)

    (
        export OPENBLAS_NUM_THREADS=1
        GMM_DEMUX_PATH=\$(type -p GMM-demux)
        export PATH=\$(dirname \$GMM_DEMUX_PATH):\$PATH
        GMM-demux demuxfile/ \$hto_string -f FULL
    )

    demultiplex_filter.py --adata ${adata_path} --demux_report FULL/GMM_full.csv --demux_config FULL/GMM_full.config --filtered_output \${adata_name}_hashing_filtered_demux.h5ad --unfiltered_output \${adata_name}_hashing_unfiltered_demux.h5ad
    """
}
