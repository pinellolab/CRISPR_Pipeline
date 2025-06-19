process prepare_assignment {
    cache 'lenient'
    debug true

    input:
        path mudata_input

    output:
        path "*_mudata.h5mu", emit: prepare_assignment_mudata

    script:
    """
        prepare_assignment.py --input ${mudata_input} 
    """
}