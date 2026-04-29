process PREPAREGUIDECHEMISTRY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    path raw_technology_file
    val  is_10x3v3
    val  spacer_tag

    output:
    path "kb_technology.txt", emit: technology_file
    path "kb_workflow.txt",   emit: workflow_file
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prepare_guide_chemistry.py \\
        --technology-file ${raw_technology_file} \\
        --is-10x3v3 ${is_10x3v3} \\
        --spacer-tag '${spacer_tag}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch kb_technology.txt
    touch kb_workflow.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0"
    END_VERSIONS
    """
}