
process evaluation_controls {

    cache 'lenient'
    input:

    path mdata


    output:
    path "plots" , emit: evaluation_controls

    script:
            """
            
            evaluate_controls.py ${mdata}

            """
}
