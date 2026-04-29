process hashing_concat {
    input:
    path hashing_demux_anndata
    path hashing_unfiltered_demux_anndata

    output:
    path 'concatenated_hashing_demux.h5ad',            emit: concatenated_hashing_demux
    path 'concatenated_unfiltered_hashing_demux.h5ad', emit: concatenated_hashing_unfiltered_demux

    script:
    """
    cp ${hashing_demux_anndata} concatenated_hashing_demux.h5ad
    cp ${hashing_unfiltered_demux_anndata} concatenated_unfiltered_hashing_demux.h5ad
    """
}