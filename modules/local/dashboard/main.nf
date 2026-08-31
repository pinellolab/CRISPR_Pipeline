process createDashboard {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        return "${out}/pipeline_dashboard"
    }, mode: 'copy', saveAs: { filename -> filename in ['pipeline_dashboard.tar.gz', 'pipeline_qc_metrics.json'] ? null : filename }
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        return out.replaceAll('/$','')
    }, mode: 'copy', overwrite: true, pattern: 'pipeline_dashboard.tar.gz'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        return out.replaceAll('/$','')
    }, mode: 'copy', overwrite: true, pattern: 'pipeline_qc_metrics.json'

    input:
        path guide_seqSpecCheck_plots
        path guide_fq_tbl
        path mudata
        path gene_ann
        path gene_ann_filtered
        path guide_ann
        path ks_transcripts_out_dir_collected
        path ks_guide_out_dir_collected
        path additional_qc_dir
        path figures_dir
        path evaluation_output_dir
        path css
        path js
        path svg
        path controls_evaluation_output_dir
        path benchmark_output_dir
        

    output:
        tuple path("evaluation_output"), path("figures"), path("guide_seqSpec_plots"), path("additional_qc"), path("dashboard.html"), path("svg"), path("benchmark_output")
        path "pipeline_dashboard.tar.gz", emit: dashboard_archive
        path "pipeline_qc_metrics.json", emit: qc_metrics_json

    script:
        def dashboard_params_json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
        """
        cat > dashboard_params.json <<'EOF_DASHBOARD_PARAMS'
${dashboard_params_json}
EOF_DASHBOARD_PARAMS

        echo "=== RENAMING INPUT DIRECTORIES ==="
        [[ -e guide_seqSpec_plots ]] && mv guide_seqSpec_plots input_guide_seqSpec_plots
        [[ -e figures ]] && mv figures input_figures
        [[ -e evaluation_output ]] && mv evaluation_output input_evaluation_output
        [[ -e svg ]] && mv svg input_svg
        [[ -e additional_qc ]] && mv additional_qc input_additional_qc
        [[ -e benchmark_output ]] && mv benchmark_output input_benchmark_output
        [[ -e ${mudata} ]] && mv ${mudata} input_mudata.h5mu
        [[ -e plots ]] && mv plots input_controls_evaluation_output_dir


        # Create new output directories with actual content
        echo "=== CREATING OUTPUT DIRECTORIES ==="
        mkdir -p guide_seqSpec_plots figures evaluation_output svg additional_qc benchmark_output

        # Copy content from renamed inputs to new outputs
        if [[ -L input_guide_seqSpec_plots ]]; then
            cp -rL input_guide_seqSpec_plots/* guide_seqSpec_plots/ 2>/dev/null || true
        fi
        if [[ -L input_figures ]]; then
            cp -rL input_figures/* figures/ 2>/dev/null || true
        fi
        if [[ -L input_evaluation_output ]]; then
            cp -rL input_evaluation_output/* evaluation_output/ 2>/dev/null || true
        fi
        if [[ -L input_svg ]]; then
            cp -rL input_svg/* svg/ 2>/dev/null || true
        fi
        if [[ -L input_additional_qc ]]; then
            cp -rL input_additional_qc/* additional_qc/ 2>/dev/null || true
        fi
        if [[ -L input_benchmark_output ]]; then
            cp -rL input_benchmark_output/* benchmark_output/ 2>/dev/null || true
        fi
        if [[ -L input_controls_evaluation_output_dir ]]; then
            cp -rL input_controls_evaluation_output_dir/* evaluation_output/ 2>/dev/null || true
        fi
  
        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        # Run scripts using the renamed input files
        process_json.py --output_dir json_dir
        create_dashboard_plots.py --mudata input_mudata.h5mu --output_dir figures
        create_dashboard_df.py --json_dir json_dir --guide_fq_tbl ${guide_fq_tbl} --mudata input_mudata.h5mu --gene_ann ${gene_ann} --gene_ann_filtered ${gene_ann_filtered} --guide_ann ${guide_ann} --params_json dashboard_params.json --additional_qc_dir additional_qc --benchmark_dir benchmark_output --qc_metrics_json pipeline_qc_metrics.json
        create_dashboard.py --input all_df.pkl

        mkdir -p pipeline_dashboard
        cp -R evaluation_output figures guide_seqSpec_plots additional_qc dashboard.html svg benchmark_output pipeline_dashboard/
        tar -czf pipeline_dashboard.tar.gz pipeline_dashboard
        rm -rf pipeline_dashboard

        echo "=== FINAL OUTPUT SIZES ==="
        du -sh guide_seqSpec_plots figures evaluation_output svg additional_qc benchmark_output *.html pipeline_dashboard.tar.gz pipeline_qc_metrics.json 2>/dev/null || true
        """
}
