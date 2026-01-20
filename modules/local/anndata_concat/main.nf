process anndata_concat {
    cache 'lenient'
    debug true

    input:
    path adata_filepath
    path parsed_covariate_df
    val modality

    output:
    path "concatenated_adata.h5ad", emit: concat_anndata

    script:
    """
    anndata_concat.py ${adata_filepath} ${parsed_covariate_df} --output concatenated_adata.h5ad --modality ${modality}
    """
}
