
process createDashboard_HASHING {
    cache 'lenient'
    publishDir './pipeline_dashboard'

    input:
        path guide_seqSpecCheck_plots
        path guide_fq_tbl
        path hashing_seqSpecCheck_plots
        path hashing_fq_tbl
        path mudata
        path gene_ann
        path gene_ann_filtered
        path guide_ann
        path hashing_ann
        path hashing_demux
        path hashing_unfiltered_demux
        path ks_transcripts_out_dir_collected
        path ks_guide_out_dir_collected
        path ks_hashing_out_dir_collected
        path figures_dir
        path evaluation_output_dir
        path css
        path js
        path svg

    output:
        tuple path("evaluation_output"), path("figures"), path("guide_seqSpec_plots"), path("hashing_seqSpec_plots"), path("dashboard.html"), path("svg"), path("inference_mudata.h5mu")

    script:
        """
        echo "Guide seqSpec plots directory: ${guide_seqSpecCheck_plots}"
        echo "Hashing seqSpec plots directory: ${hashing_seqSpecCheck_plots}"
        echo "Transcriptome output directory: ${ks_transcripts_out_dir_collected}"
        echo "Guide output directory: ${ks_guide_out_dir_collected}"
        echo "Hashing output directory: ${ks_hashing_out_dir_collected}"
        echo "Figures directory: ${figures_dir}"
        echo "Evaluation output directory: ${evaluation_output_dir}"
        echo "css directory: ${css}"
        echo "js directory: ${js}"
        echo "svg directory: ${svg}"

        export MPLCONFIGDIR="./tmp/mplconfigdir"
        process_json_HASHING.py --output_dir json_dir
        create_dashboard_plots_HASHING.py --mudata ${mudata} --hashing_demux ${hashing_demux} --unfiltered_hashing_demux ${hashing_unfiltered_demux} --output_dir figures
        create_dashboard_df_HASHING.py --json_dir json_dir --guide_fq_tbl ${guide_fq_tbl} --hashing_fq_tbl ${hashing_fq_tbl} --mudata ${mudata} --gene_ann ${gene_ann} --gene_ann_filtered ${gene_ann_filtered} --guide_ann ${guide_ann} --hashing_ann ${hashing_ann} --hashing_demux ${hashing_demux} --hashing_unfiltered_demux ${hashing_unfiltered_demux}
        create_dashboard_HASHING.py --input all_df.pkl
        """

}
