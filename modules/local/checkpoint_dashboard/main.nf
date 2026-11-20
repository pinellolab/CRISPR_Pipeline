process save_checkpoint_dashboard {
    container { params.containers.base }
    publishDir { "${params.outdir}/preinference_checkpoint/${workflow.runName}" }, mode: 'copy', overwrite: false

    input:
    path manifest
    path qc_metrics
    path seqspec_table
    path qc_dirs

    output:
    path 'checkpoint_qc'

    script:
    """
    python ${projectDir}/bin/checkpoint_dashboard.py snapshot --manifest "${manifest}" \\
      --metrics "${qc_metrics}" --seqspec "${seqspec_table}" --output checkpoint_qc ${qc_dirs}
    """
}

process stage_checkpoint_dashboard {
    container { params.containers.base }
    cache 'deep'
    publishDir { params.outdir }, mode: 'copy', overwrite: true

    input:
    path external_report
    path sources, stageAs: 'parent/*'
    path guide_report, stageAs: 'guide/*'

    output:
    path 'checkpoint_dashboard', emit: dashboard

    script:
    """
    python ${projectDir}/bin/checkpoint_dashboard.py restore --report "${external_report}" \\
      --source parent --guide-report guide/guide_metadata.validation.json --output checkpoint_dashboard
    """
}

process render_continuation_dashboard {
    container { params.containers.base }
    publishDir { "${params.outdir}/pipeline_dashboard" }, mode: 'copy', overwrite: true

    input:
    path checkpoint_dashboard
    path inference_mudata
    path evaluation_output
    path controls, stageAs: 'controls/*'
    path calibration, stageAs: 'calibration_input/*'

    output:
    path 'dashboard.html'
    path 'inference_preview.html'
    path 'evaluation_output'

    script:
    """
    mv "${evaluation_output}" staged_evaluation
    mkdir evaluation_output
    cp -RL staged_evaluation/. evaluation_output/
    cp -RL "${controls}/." evaluation_output/
    if [[ "${calibration}" == */calibration ]]; then
        mkdir -p pipeline_outputs/calibration
        cp -RL "${calibration}/." pipeline_outputs/calibration/
    fi

    python ${projectDir}/bin/checkpoint_dashboard.py complete --root . \\
      --mudata "${inference_mudata}" --run-name "${workflow.runName}"
    """
}
