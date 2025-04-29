
process PreprocessAnnData {
  
  cache 'lenient'

  input:
  path adata_rna
  path gname_rna
  val min_genes
  val min_cells
  val pct_mito
  val reference

  output:
  path "filtered_anndata.h5ad" , emit: filtered_anndata_rna
  path "rna_concatenated_adata.h5ad", emit: adata_rna
  path "figures", emit: figures_dir
  
  script:
        """
        preprocess_adata.py ${adata_rna} ${gname_rna} --min_genes ${min_genes} --min_cells ${min_cells} --pct_mito ${pct_mito} --reference ${reference}
        mv concatenated_adata.h5ad rna_concatenated_adata.h5ad
        """

}
