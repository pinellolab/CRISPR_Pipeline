process GUIDE_ASSIGNMENT_STUB {
    input:
    path mudata

    output:
    path 'concat_mudata.h5mu', emit: concat_mudata

    script:
    """
    cp ${mudata} concat_mudata.h5mu
    """
}

workflow guide_assignment_pipeline {
    take:
    mudata

    main:
    GUIDE_ASSIGNMENT_STUB(mudata)

    emit:
    concat_mudata = GUIDE_ASSIGNMENT_STUB.out.concat_mudata
}