
process filter_hashing{

    cache 'lenient'
    debug true

    input:
    path filtered_anndata_rna
    path anndata_hashing

    output:
    path "*_filtered.h5ad", emit: hashing_filtered_anndata
    path "hashing_concatenated_adata.h5ad", emit: adata_hashing

    script:
    """
    filter_hashing.py --hash_file ${anndata_hashing} --rna_file ${filtered_anndata_rna}
    mv concatenated_adata.h5ad hashing_concatenated_adata.h5ad
    """
}
