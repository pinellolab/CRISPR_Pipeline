process CreateMuData {
    input:
    path adata_rna
    path adata_guide
    path guide_metadata
    path gtf_file
    val  moi
    val  capture_method
    path adata_hashing

    output:
    path 'mudata.h5mu',                 emit: mudata
    path 'guide_concatenated_adata.h5ad', emit: adata_guide

    script:
    """
    touch mudata.h5mu
    cp ${adata_guide} guide_concatenated_adata.h5ad
    """
}