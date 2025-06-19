
process CreateMuData_HASHING {
    cache 'lenient'

    input:

    path adata_rna
    path adata_guide
    path adata_hashing
    path guide_metadata
    path gtf_file
    val moi
    val capture_method

    output:
    path "mudata.h5mu" , emit: mudata
    path "guide_concatenated_adata.h5ad", emit: adata_guide

    script:
        """
        create_mdata_HASHING.py ${adata_rna} ${adata_guide} ${adata_hashing} ${guide_metadata} ${gtf_file} ${moi} ${capture_method}
        mv concatenated_adata.h5ad guide_concatenated_adata.h5ad
        """

}
