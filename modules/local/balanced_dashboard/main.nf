process balanced_dashboard {
    cache 'lenient'
    cpus 2
    memory '64 GB'
    time '12h'
    publishDir "${params.outdir}/pipeline_dashboard", mode: 'copy', pattern: 'dashboard.html'
    publishDir "${params.outdir}/pipeline_dashboard", mode: 'copy', pattern: 'additional_qc'
    publishDir "${params.outdir}/pipeline_dashboard", mode: 'copy', pattern: 'evaluation_output'
    publishDir "${params.outdir}/pipeline_dashboard", mode: 'copy', pattern: 'control_output'
    publishDir "${params.outdir}/pipeline_dashboard", mode: 'copy', pattern: 'benchmark_output'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'pipeline_qc_metrics.json'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'pipeline_dashboard.tar.gz'
    input:
    path mudata
    path report
    path additional, stageAs: 'inputs/additional_qc'
    path evaluation, stageAs: 'inputs/evaluation_output'
    path controls, stageAs: 'inputs/control_output'
    path benchmark, stageAs: 'inputs/benchmark_output'
    output:
    path 'pipeline_dashboard.tar.gz', emit: archive
    path 'pipeline_qc_metrics.json', emit: metrics
    path 'dashboard.html'
    path 'additional_qc'
    path 'evaluation_output'
    path 'control_output'
    path 'benchmark_output'
    script:
    def resolved_params = groovy.json.JsonOutput.toJson(params)
    """
    cat > resolved_params.json <<'PARAMS_EOF'
${resolved_params}
PARAMS_EOF
    cp -RL inputs/. .
    python ${projectDir}/bin/balanced_dashboard.py --input ${mudata} --report ${report} --params resolved_params.json
    mkdir pipeline_dashboard
    cp -R dashboard.html additional_qc evaluation_output control_output benchmark_output pipeline_dashboard/
    tar -czf pipeline_dashboard.tar.gz pipeline_dashboard
    """
}
