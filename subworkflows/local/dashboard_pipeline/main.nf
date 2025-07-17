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
    figures_dir
    evaluation_output_dir

    main:

    if (params.INFERENCE_method == 'default') {
        createDashboard_default(
            guide_seqSpecCheck_plots,
            guide_position_table,
            mdata,
            concat_anndata_rna,
            filtered_anndata_rna,
            concat_anndata_guide,
            ks_transcripts_out_dir_collected,
            ks_guide_out_dir_collected,
            figures_dir,
            evaluation_output_dir,
            file(params.css),
            file(params.js),
            file(params.svg)
        )
    } else {
        createDashboard(
            guide_seqSpecCheck_plots,
            guide_position_table,
            mdata,
            concat_anndata_rna,
            filtered_anndata_rna,
            concat_anndata_guide,
            ks_transcripts_out_dir_collected,
            ks_guide_out_dir_collected,
            figures_dir,
            evaluation_output_dir,
            file(params.css),
            file(params.js),
            file(params.svg)
        )
    }

}
