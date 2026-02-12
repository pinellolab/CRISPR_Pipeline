process additional_qc_plots {
    cache 'lenient'

    input:
        path mudata

    output:
        path "additional_qc", emit: additional_qc

    script:
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        mkdir -p additional_qc/gene additional_qc/guide additional_qc/intended_target additional_qc/trans

        mapping_gene.py \\
            --input ${mudata} \\
            --outdir additional_qc/gene \\
            --batch-col ${params.QC_batch_col ?: 'batch'}

        mapping_guide.py \\
            --input ${mudata} \\
            --outdir additional_qc/guide \\
            --batch-col ${params.QC_batch_col ?: 'batch'}

        HAS_RESULTS=\$(additional_has_results.py --input ${mudata})

        if [[ "\${HAS_RESULTS}" == "YES" ]]; then
            intended_target.py \\
                --input ${mudata} \\
                --outdir additional_qc/intended_target \\
                --results-key auto \\
                --log2fc-col auto \\
                --pvalue-col auto

            trans.py \\
                --input ${mudata} \\
                --outdir additional_qc/trans \\
                --results-key auto \\
                --log2fc-col auto \\
                --pvalue-col auto
        else
            echo "No inference results found in mudata.uns; skipping intended_target/trans QC."
        fi
        """
}
