include { DASHBOARD as DASHBOARD_STD         } from '../../../modules/local/dashboard'
include { DASHBOARD as DASHBOARD_STD_DEFAULT } from '../../../modules/local/dashboard'

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
    def ch_empty = channel.value([])

    if (params.INFERENCE_method == 'default') {
        DASHBOARD_STD_DEFAULT(
            guide_seqSpecCheck_plots,
            guide_position_table,
            ch_empty,
            ch_empty,
            mdata,
            concat_anndata_rna,
            filtered_anndata_rna,
            concat_anndata_guide,
            ch_empty,
            ch_empty,
            ch_empty,
            ks_transcripts_out_dir_collected,
            ks_guide_out_dir_collected,
            ch_empty,
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
        DASHBOARD_STD(
            guide_seqSpecCheck_plots,
            guide_position_table,
            ch_empty,
            ch_empty,
            mdata,
            concat_anndata_rna,
            filtered_anndata_rna,
            concat_anndata_guide,
            ch_empty,
            ch_empty,
            ch_empty,
            ks_transcripts_out_dir_collected,
            ks_guide_out_dir_collected,
            ch_empty,
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
    dashboard_output = (params.INFERENCE_method == 'default' ? DASHBOARD_STD_DEFAULT.out : DASHBOARD_STD.out)
}