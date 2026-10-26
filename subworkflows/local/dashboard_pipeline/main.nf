nextflow.enable.dsl=2

include { createDashboard } from '../../../modules/local/dashboard'
include { createDashboard_default } from '../../../modules/local/dashboard_default'

workflow dashboard_pipeline {
    take:
    guide_seqSpecCheck_plots
    guide_position_table
    concat_anndata_rna
    filtered_anndata_rna
    ks_transcripts_out_dir_collected
    concat_anndata_guide
    ks_guide_out_dir_collected
    mdata
    additional_qc_dir
    figures_dir
    evaluation_output_dir
    controls_evaluation_output_dir
    benchmark_output_dir


    main:

    if (params.INFERENCE_method == 'default') {
        DashboardProcess = createDashboard_default(
            guide_seqSpecCheck_plots,
            guide_position_table,
            mdata,
            concat_anndata_rna,
            filtered_anndata_rna,
            concat_anndata_guide,
            ks_transcripts_out_dir_collected,
            ks_guide_out_dir_collected,
            additional_qc_dir,
            figures_dir,
            evaluation_output_dir,
            file(params.css),
            file(params.js),
            file(params.svg),
            controls_evaluation_output_dir,
            benchmark_output_dir
        )
    } else {
        DashboardProcess = createDashboard(
            guide_seqSpecCheck_plots,
            guide_position_table,
            mdata,
            concat_anndata_rna,
            filtered_anndata_rna,
            concat_anndata_guide,
            ks_transcripts_out_dir_collected,
            ks_guide_out_dir_collected,
            additional_qc_dir,
            figures_dir,
            evaluation_output_dir,
            file(params.css),
            file(params.js),
            file(params.svg),
            controls_evaluation_output_dir,
            benchmark_output_dir
        )
    }

    emit:
    dashboard_archive = DashboardProcess.dashboard_archive

}
