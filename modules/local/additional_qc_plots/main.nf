process additional_qc_plots {
    cache 'lenient'

    input:
        path mudata
        path clone_qc_input
        path saturation_qc_input

    output:
        path "additional_qc", emit: additional_qc

    script:
        """
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        mkdir -p additional_qc/gene additional_qc/guide additional_qc/intended_target additional_qc/global_analysis

        if find ${clone_qc_input} -mindepth 1 -type f ! -name '.gitkeep' -print -quit | grep -q .; then
            mkdir -p additional_qc/clones
            cp -R ${clone_qc_input}/. additional_qc/clones/
        fi
        if find ${saturation_qc_input} -mindepth 1 -type f ! -name '.gitkeep' -print -quit | grep -q .; then
            mkdir -p additional_qc/sequencing_saturation
            cp -R ${saturation_qc_input}/. additional_qc/sequencing_saturation/
        fi

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
                --outdir additional_qc/global_analysis \\
                --results-key auto \\
                --log2fc-col auto \\
                --pvalue-col auto \\
                --prefix global_analysis
        else
            echo "No inference results found in mudata.uns; skipping local/global analysis QC."
        fi
        """
}
