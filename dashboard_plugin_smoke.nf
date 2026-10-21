nextflow.enable.dsl=2

process DASHBOARD_PLUGIN_SMOKE {
    output:
    path 'smoke.txt'

    script:
    '''
    sleep 15
    printf 'dashboard plugin smoke test\n' > smoke.txt
    '''
}

workflow {
    DASHBOARD_PLUGIN_SMOKE()
}
