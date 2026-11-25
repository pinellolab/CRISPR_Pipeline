process anndata_concat {
    cache 'lenient'

    input:
    path parsed_covariate_df
    path adata_filepath
    // Whether this modality's mapping wrote counts_unfiltered_modified. One
    // module serves RNA, guide and hashing, and only the modalities whose
    // barcodes were translated have that folder.
    val bc_replacement

    output:
    path "concatenated_adata.h5ad", emit: concat_anndata

    script:
    """
    anndata_concat.py ${adata_filepath} ${parsed_covariate_df} --output concatenated_adata.h5ad --bc_replacement ${bc_replacement} --mm ${params.use_multimapping} --normalize-feature-index
    """
}
