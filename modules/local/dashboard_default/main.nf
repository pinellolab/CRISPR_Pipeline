process createDashboard_default {
    cache 'lenient'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        return "${out}/pipeline_dashboard"
    }, mode: 'copy'
    publishDir path: {
        def out = params.outdir?.toString() ?: './pipeline_outputs'
        out = out.replaceAll('/$','')
        if (out == 'pipeline_outputs' || out.endsWith('/pipeline_outputs')) {
            return out
        }
        return "${out}/pipeline_outputs"
    }, mode: 'copy', overwrite: true, pattern: 'inference_mudata.h5mu'

    input:
        path guide_seqSpecCheck_plots
        path guide_fq_tbl
        path mudata
        path gene_ann
        path gene_ann_filtered
        path guide_ann
        path ks_transcripts_out_dir_collected
        path ks_guide_out_dir_collected
        path figures_dir
        path evaluation_output_dir
        path css
        path js
        path svg
        path controls_evaluation_output_dir

    output:
        tuple path("evaluation_output"), path("figures"), path("guide_seqSpec_plots"), path("dashboard.html"), path("svg"), path("inference_mudata.h5mu")

    script:
        """
        echo "=== RENAMING INPUT DIRECTORIES ==="
        [[ -e guide_seqSpec_plots ]] && mv guide_seqSpec_plots input_guide_seqSpec_plots
        [[ -e figures ]] && mv figures input_figures
        [[ -e evaluation_output ]] && mv evaluation_output input_evaluation_output
        [[ -e svg ]] && mv svg input_svg
        [[ -e ${mudata} ]] && mv ${mudata} input_mudata.h5mu
        [[ -e plots ]] && mv plots input_controls_evaluation_output_dir


        # Create new output directories with actual content
        echo "=== CREATING OUTPUT DIRECTORIES ==="
        mkdir -p guide_seqSpec_plots figures evaluation_output svg

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
        if [[ -f input_mudata.h5mu ]]; then
            cp -L input_mudata.h5mu inference_mudata.h5mu
        fi

        if [[ -L input_controls_evaluation_output_dir ]]; then
            cp -rL input_controls_evaluation_output_dir/* evaluation_output/ 2>/dev/null || true
        fi

        export MPLCONFIGDIR="./tmp/mplconfigdir"
        mkdir -p \${MPLCONFIGDIR}

        # Run scripts using the renamed input files
        process_json.py --output_dir json_dir
        create_dashboard_plots.py --mudata input_mudata.h5mu --output_dir figures
        create_dashboard_df.py --json_dir json_dir --guide_fq_tbl ${guide_fq_tbl} --mudata input_mudata.h5mu --gene_ann ${gene_ann} --gene_ann_filtered ${gene_ann_filtered} --guide_ann ${guide_ann} --default
        create_dashboard.py --input all_df.pkl

        echo "=== FINAL OUTPUT SIZES ==="
        du -sh guide_seqSpec_plots figures evaluation_output svg *.html *.h5mu 2>/dev/null || true
        """
}
