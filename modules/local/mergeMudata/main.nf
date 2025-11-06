process mergeMudata {
    cache 'lenient'
    debug true
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
        path cis_per_guide
        path cis_per_element
        path trans_per_guide
        path trans_per_element
        path base_mudata

    output:
        path "inference_mudata.h5mu", emit: inference_mudata
        path "cis_per_guide_results.tsv.gz", emit: cis_per_guide_results
        path "cis_per_element_results.tsv.gz", emit: cis_per_element_results
        path "trans_per_guide_results.tsv.gz", emit: trans_per_guide_results
        path "trans_per_element_results.tsv.gz", emit: trans_per_element_results

    script:
    """
        merge_cis_trans_results.py \\
            --cis_per_guide ${cis_per_guide} \\
            --cis_per_element ${cis_per_element} \\
            --trans_per_guide ${trans_per_guide} \\
            --trans_per_element ${trans_per_element} \\
            --base_mudata ${base_mudata} \\
            --output inference_mudata.h5mu
    """
}
