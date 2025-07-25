process mergeMudata {
    cache 'lenient'
    debug true
    publishDir './pipeline_outputs', mode: 'copy', overwrite: true
    
    input:
        path cis_file
        path trans_file

    output:
        path "inference_mudata.h5mu", emit: inference_mudata

    script:
    """
        merge_mudata.py --cis_mudata ${cis_file} --trans_mudata ${trans_file} --output inference_mudata.h5mu
    """
}
