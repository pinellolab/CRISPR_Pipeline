process doublets_scrub {
    input:
    path mudata

    output:
    path 'mudata_doublet.h5mu', emit: mudata_doublet

    script:
    """
    cp ${mudata} mudata_doublet.h5mu
    """
}