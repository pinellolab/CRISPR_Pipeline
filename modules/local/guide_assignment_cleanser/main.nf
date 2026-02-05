process guide_assignment_cleanser {
    cache 'lenient'
    debug true

    input:
        path mudata_input
        val threshold
        val capture_method

    output:
        path "${mudata_input.simpleName}_output.h5mu", emit: guide_assignment_mudata_output

    script:
        def thresh_opt = threshold ? "-t ${threshold}" : ""
        """
        export CMDSTAN=/root/.cmdstan/cmdstan-2.36.0
        export PATH=\$PATH:\$CMDSTAN/bin
        export TMPDIR='./tmp'
        mkdir -p ./tmp
        cat /usr/local/lib/python3.11/site-packages/cleanser/configuration.py
        cleanser -i ${mudata_input} --posteriors-output ${mudata_input.simpleName}_output.h5mu --modality guide --${capture_method} --output-layer guide_assignment ${thresh_opt}
        """
}
