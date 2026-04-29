process PREPROCESSING_STUB {
    input:
    path concat_anndata_rna
    path trans_out_dir

    output:
    path 'filtered_anndata_rna.h5ad', emit: filtered_anndata_rna
    path 'adata_rna.h5ad',            emit: adata_rna
    path 'figures',                   emit: figures_dir
    path 'gencode.gtf',               emit: gencode_gtf

    script:
    """
    cp ${concat_anndata_rna} filtered_anndata_rna.h5ad
    cp ${concat_anndata_rna} adata_rna.h5ad
    mkdir -p figures
    touch figures/placeholder.txt
    cat <<'EOF' > gencode.gtf
chr1\tstub\tgene\t1\t1000\t.\t+\t.\tgene_id "GENE0"; gene_name "GENE0";
EOF
    """
}

workflow preprocessing_pipeline {
    take:
    concat_anndata_rna
    trans_out_dir

    main:
    PREPROCESSING_STUB(concat_anndata_rna, trans_out_dir)

    emit:
    filtered_anndata_rna = PREPROCESSING_STUB.out.filtered_anndata_rna
    figures_dir          = PREPROCESSING_STUB.out.figures_dir
    adata_rna            = PREPROCESSING_STUB.out.adata_rna
    gencode_gtf          = PREPROCESSING_STUB.out.gencode_gtf
}