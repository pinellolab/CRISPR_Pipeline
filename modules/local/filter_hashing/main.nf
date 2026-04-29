process filter_hashing {
    input:
    path filtered_anndata_rna
    path anndata_hashing

    output:
    path 'hashing_filtered.h5ad',      emit: hashing_filtered_anndata
    path 'hashing_concatenated_adata.h5ad', emit: adata_hashing

    script:
    """
    cp ${anndata_hashing} hashing_filtered.h5ad
    cp ${anndata_hashing} hashing_concatenated_adata.h5ad
    """
}