nextflow.enable.dsl=2

include { createDashboard_HASHING } from '../../../modules/local/dashboard_HASHING'

workflow dashboard_pipeline_HASHING {
    take:
    guide_seqSpecCheck_plots
    guide_position_table
    hashing_seqSpecCheck_plots
    hashing_position_table
    concat_anndata_rna
    filtered_anndata_rna
    ks_transcripts_out_dir_collected
    concat_anndata_guide
    ks_guide_out_dir_collected
    concat_anndata_hashing
    ks_hashing_out_dir_collected
    adata_demux
    adata_unfiltered_demux
    mdata
    figures_dir
    evaluation_output_dir

    main:

    createDashboard_HASHING(
        guide_seqSpecCheck_plots,
        guide_position_table,
        hashing_seqSpecCheck_plots,
        hashing_position_table,
        mdata,
        concat_anndata_rna,
        filtered_anndata_rna,
        concat_anndata_guide,
        concat_anndata_hashing,
        adata_demux,
        adata_unfiltered_demux,
        ks_transcripts_out_dir_collected,
        ks_guide_out_dir_collected,
        ks_hashing_out_dir_collected,
        figures_dir,
        evaluation_output_dir,
        file(params.css),
        file(params.js),
        file(params.svg)
    )

}
