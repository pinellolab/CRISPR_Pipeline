process PREPARECOVARIATE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    val covariate_list  // val: JSON-serialisable map/list of batch covariates, e.g. [batch: ['A', 'B']]

    output:
    path "parse_covariate.csv", emit: parsed_covariate_file
    path "cov_string.txt",      emit: covariate_string
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def json_string = groovy.json.JsonOutput.toJson(covariate_list)
    """
    parse_covariate.py '${json_string}'
    prepare_formula.py parse_covariate.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch parse_covariate.csv
    touch cov_string.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
