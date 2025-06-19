
process doublets_scrub {

    cache 'lenient'
    
    input:
        path mudata

    output:
        path "mdata_doublets.h5mu", emit: mudata_doublet

    script:
        """
        doublets.py --input ${mudata}
        """
}
