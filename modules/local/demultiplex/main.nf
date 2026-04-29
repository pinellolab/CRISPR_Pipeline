process demultiplex {
    input:
    path hashing_filtered_anndata

    output:
    path 'hashing_demux.h5ad',            emit: hashing_demux_anndata
    path 'hashing_demux_unfiltered.h5ad', emit: hashing_demux_unfiltered_anndata

    script:
    """
    cp ${hashing_filtered_anndata} hashing_demux.h5ad
    cp ${hashing_filtered_anndata} hashing_demux_unfiltered.h5ad
    """
}