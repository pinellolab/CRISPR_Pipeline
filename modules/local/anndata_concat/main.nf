process anndata_concat {
    cache 'lenient'
    debug true

    input:
    path parsed_covariate_df
    path adata_filepath

    output:
    path "concatenated_adata.h5ad", emit: concat_anndata

    script:
    """
    anndata_concat.py ${adata_filepath} ${parsed_covariate_df} --output concatenated_adata.h5ad
    """
}
