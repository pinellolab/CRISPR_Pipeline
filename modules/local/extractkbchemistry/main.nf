process EXTRACTKBCHEMISTRY {
    tag "$parsed_seqspec_file"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    path parsed_seqspec_file

    output:
    path "kb_technology.txt", emit: technology_file
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    extract_parsed_seqspec.py --file ${parsed_seqspec_file} > kb_technology.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch kb_technology.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "0"
        pandas: "0"
    END_VERSIONS
    """
}