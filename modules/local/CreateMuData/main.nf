
process CreateMuData {
    cache 'lenient'
    debug true

    input:

    path adata_rna
    path adata_guide
    path guide_metadata
    path gtf_file
    val moi
    val capture_method
    path adata_hashing

    output:
    path "mudata.h5mu" , emit: mudata
    path "guide_concatenated_adata.h5ad", emit: adata_guide
    path "debug_create_mdata", optional: true, emit: debug_create_mdata

    script:
        def hashing_arg = adata_hashing.name != 'NO_FILE' ? "--adata_hashing ${adata_hashing}" : ""
        def debug_arg = params.DEBUG_VAR ? "--debug_var" : ""
        """
        create_mdata.py ${adata_rna} ${adata_guide} ${guide_metadata} ${gtf_file} ${moi} ${capture_method} ${hashing_arg} ${debug_arg}
        mv concatenated_adata.h5ad guide_concatenated_adata.h5ad
        """

}
