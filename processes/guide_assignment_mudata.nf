process guide_assignment_mudata {
    
    cache 'lenient'
    
    input:
    path guide_assignment_matrix
    path mudata

    output:
        path "sceptre_assignment_mudata.h5mu", emit: guide_assignment_mudata_output

    script:
        """
          add_guide_assignment.py --guide_assignment ${guide_assignment_matrix} --mudata ${mudata} 
        """

}
