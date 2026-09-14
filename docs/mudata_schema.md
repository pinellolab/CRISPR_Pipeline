# MuData field reference

The final `pipeline_outputs/inference_mudata.h5mu` is a MuData object whose rows are retained cell barcodes and whose principal modalities are `gene` and `guide`. A third `hashing` modality is present when data hashing is enabled. RNA and guide `.X` matrices contain raw UMI counts; the binary guide calls used for inference are stored separately in `guide.layers['guide_assignment']`.

The canonical multi-tab field dictionary is also available as a [Google Sheet](https://docs.google.com/spreadsheets/d/1hwGyxCtwwgnpzEgdt7BlkK1tuJ22BPMpY-BEtCpTKm0/edit). The checked-in source for both representations is [`mudata_schema_catalog.tsv`](mudata_schema_catalog.tsv). Regenerate this page with `python3 bin/render_mudata_schema_docs.py`.

## Scope and stability

Fields marked **Always** are part of the current pipeline contract. Optional fields depend on hashing, Scrublet, guide-assignment method, inference method, or available metadata. The pipeline deliberately preserves additional `guide.var` columns from the input guide metadata; therefore dataset-specific columns may appear beyond the currently observed extensions documented below. Current output keys use `local_analysis_*` and `global_analysis_*`; `cis_*` and `trans_*` are legacy aliases retained for interpretation of older files.

The catalog was audited against `dev` source commit `33d068dad7f8163c313bd6fa62aa45c6b35bffa2` and completed SCEPTRE, CLEANSER, enhancer-screen, Gasperini, and Replogle outputs.

## Structure

| Path | Field | Type | Availability | Description | Producer/source |
|---|---|---|---|---|---|
| MuData | obs_names | string index | Always | Canonical final cell-barcode index shared by all modalities after barcode intersection and filtering. | create_mdata.py |
| MuData | var_names | string index | Always | Union feature index maintained by MuData; use modality-specific var_names for biological identifiers. | MuData |
| MuData.mod | gene | AnnData | Always | RNA count and annotation modality with cells by genes. | create_mdata.py |
| MuData.mod | guide | AnnData | Always | Guide UMI, assignment, and guide-design modality with cells by guides. | create_mdata.py |
| MuData.mod | hashing | AnnData | When ENABLE_DATA_HASHING=true | Hashtag-oligo count and demultiplexing modality aligned to the retained cells. | create_mdata.py |
| MuData.obsm | gene | boolean membership matrix | Always | MuData-generated mapping from global observations to gene-modality observations. | MuData |
| MuData.obsm | guide | boolean membership matrix | Always | MuData-generated mapping from global observations to guide-modality observations. | MuData |
| MuData.obsm | hashing | boolean membership matrix | Hashing only | MuData-generated mapping from global observations to hashing-modality observations. | MuData |
| MuData.varm | gene | boolean membership matrix | Always | MuData-generated mapping from global variables to gene features. | MuData |
| MuData.varm | guide | boolean membership matrix | Always | MuData-generated mapping from global variables to guide features. | MuData |
| MuData.varm | hashing | boolean membership matrix | Hashing only | MuData-generated mapping from global variables to hashtag features. | MuData |

## Matrices

| Path | Field | Type | Availability | Description | Producer/source |
|---|---|---|---|---|---|
| gene.X | RNA counts | sparse integer matrix | Always | Filtered cell-by-gene UMI count matrix; values remain counts rather than normalized expression. | RNA mapping and preprocessing |
| guide.X | Guide counts | sparse integer matrix | Always | Raw cell-by-guide UMI count matrix before binary guide assignment. | Guide mapping |
| guide.layers | guide_assignment | sparse binary matrix | After guide assignment | Final cell-by-guide assignment used by inference and QC; nonzero means that guide is assigned to the cell. | add_guide_assignment.py or CLEANSER integration |
| guide.layers | guide_assignment_posteriors | sparse/dense float matrix | CLEANSER when emitted | Per-cell, per-guide assignment posterior/probability retained separately from the binary assignment. | CLEANSER integration |
| hashing.X | Hash counts | sparse integer matrix | Hashing only | Cell-by-hashtag UMI count matrix carried through hashing demultiplexing. | Hash mapping/demultiplexing |

## Cell fields

| Path | Field | Type | Availability | Description | Producer/source |
|---|---|---|---|---|---|
| MuData.obs | batch | categorical | When common to all modalities | Shared batch label promoted from modality observation tables. | create_mdata.py |
| MuData.obs | concat_batch | categorical | After multi-measurement concatenation | Stable concatenation batch/source label common to modalities. | concatenation |
| gene.obs | batch | categorical | Always | Input or mapping batch label for each cell. | RNA mapping |
| gene.obs | concat_batch | categorical | After concatenation | Measurement-set/source label added during concatenation. | concatenation |
| gene.obs | batch_number | integer | Always | One-based numeric encoding of batch used as an inference/QC covariate. | preprocess_adata.py |
| gene.obs | n_counts | integer | Always | Number of detected genes; renamed from n_genes_by_counts during MuData creation. | Scanpy QC/create_mdata.py |
| gene.obs | num_expressed_genes | integer | Version/input dependent | Number of genes with a nonzero count in the cell; legacy/input alias retained when supplied. | RNA preprocessing/create_mdata.py |
| gene.obs | total_gene_umis | integer | Always | Total RNA UMI count for the cell; renamed from total_counts. | Scanpy QC/create_mdata.py |
| gene.obs | log1p_n_genes_by_counts | float | Always | Natural log of one plus the number of detected genes. | Scanpy QC |
| gene.obs | log1p_total_counts | float | Always | Natural log of one plus total RNA counts. | Scanpy QC |
| gene.obs | pct_counts_in_top_50_genes | float percent | Always | Percent of cell RNA counts contributed by its 50 most abundant genes. | Scanpy QC |
| gene.obs | pct_counts_in_top_100_genes | float percent | Always | Percent of cell RNA counts contributed by its 100 most abundant genes. | Scanpy QC |
| gene.obs | pct_counts_in_top_200_genes | float percent | Always | Percent of cell RNA counts contributed by its 200 most abundant genes. | Scanpy QC |
| gene.obs | pct_counts_in_top_500_genes | float percent | Always | Percent of cell RNA counts contributed by its 500 most abundant genes. | Scanpy QC |
| gene.obs | total_counts_mt | integer | Always | Total RNA UMIs assigned to mitochondrial genes. | Scanpy QC |
| gene.obs | log1p_total_counts_mt | float | Always | Natural log of one plus mitochondrial RNA UMIs. | Scanpy QC |
| gene.obs | percent_mito | float percent | Always | Percent of RNA counts from mitochondrial genes; renamed from pct_counts_mt. | preprocess_adata.py/create_mdata.py |
| gene.obs | total_counts_ribo | integer | Always | Total RNA UMIs assigned to ribosomal genes. | Scanpy QC |
| gene.obs | log1p_total_counts_ribo | float | Always | Natural log of one plus ribosomal RNA UMIs. | Scanpy QC |
| gene.obs | pct_counts_ribo | float percent | Always | Percent of RNA counts from ribosomal genes. | Scanpy QC |
| gene.obs | doublet_scores | float | When ENABLE_SCRUBLET=true | Scrublet doublet score before predicted doublets are removed. | doublets.py |
| gene.obs | predicted_doublets | boolean | When ENABLE_SCRUBLET=true | Scrublet binary doublet prediction used for filtering. | doublets.py |
| gene.obs | doublet_info | string | When ENABLE_SCRUBLET=true | String representation of the Scrublet prediction retained for reporting. | doublets.py |
| guide.obs | batch | categorical | Always | Input or mapping batch label for each guide-count row. | Guide mapping |
| guide.obs | concat_batch | categorical | After concatenation | Measurement-set/source label added during concatenation. | concatenation |
| guide.obs | batch_number | integer | Hashing workflows | One-based guide batch encoding added when hashing data are present. | create_mdata.py |
| guide.obs | num_expressed_guides | integer | Always | Number of guides with at least one raw guide UMI in the cell; this is not the number of assigned guides. | create_mdata.py |
| guide.obs | total_guide_umis | integer | Always | Sum of raw guide UMIs in the cell; this is not derived from guide_assignment. | create_mdata.py |
| hashing.obs | batch | categorical | Hashing only | Input batch label for each hashtag-count row. | Hash mapping |
| hashing.obs | concat_batch | categorical | Hashing after concatenation | Measurement-set/source label added during concatenation. | hashing_concat.py |
| hashing.obs | cluster_id | integer/string | Hashing after demultiplexing | Demultiplexing cluster assigned from the HTO report. | demultiplex_filter.py |
| hashing.obs | hto_type | string/categorical | Hashing after demultiplexing | HTO identity or demultiplexing class associated with cluster_id. | demultiplex_filter.py |
| hashing.obs | hto_type_split | string/categorical | Hashing after demultiplexing | Normalized HTO class; multi-HTO labels are collapsed to multiplets and negatives remain negative. | demultiplex_filter.py |

## Gene fields

| Path | Field | Type | Availability | Description | Producer/source |
|---|---|---|---|---|---|
| gene.var_names | gene_id | string index | Always | Version-stripped Ensembl gene identifier used as the canonical gene feature key. | RNA reference/preprocess_adata.py |
| gene.var | symbol | string | Always | Human-readable gene symbol derived from the RNA reference. | preprocess_adata.py |
| gene.var | mt | boolean | Always | True when the gene symbol begins with MT-. | preprocess_adata.py |
| gene.var | ribo | boolean | Always | True for ribosomal gene-symbol prefixes used by pipeline QC. | preprocess_adata.py |
| gene.var | n_cells_by_counts | integer | Always | Number of retained cells with a nonzero count for the gene. | Scanpy QC |
| gene.var | n_cells | integer | After gene filtering | Number of retained cells supporting the gene; used by the fractional minimum-cell filter. | preprocessing/aggregation |
| gene.var | mean_counts | float | Always | Mean raw RNA UMI count across retained cells. | Scanpy QC |
| gene.var | log1p_mean_counts | float | Always | Natural log of one plus mean_counts. | Scanpy QC |
| gene.var | pct_dropout_by_counts | float percent | Always | Percent of retained cells with zero counts for the gene. | Scanpy QC |
| gene.var | total_counts | integer | Always | Total raw RNA UMIs for the gene across retained cells. | Scanpy QC |
| gene.var | log1p_total_counts | float | Always | Natural log of one plus the gene total count. | Scanpy QC |
| gene.var | gene_chr | string/categorical | After GTF annotation | Chromosome/contig of the gene from the configured GTF. | create_mdata.py |
| gene.var | gene_start | integer/float nullable | After GTF annotation | Gene start coordinate from the configured GTF. | create_mdata.py |
| gene.var | gene_end | integer/float nullable | After GTF annotation | Gene end coordinate from the configured GTF. | create_mdata.py |

## Guide fields

| Path | Field | Type | Availability | Description | Producer/source |
|---|---|---|---|---|---|
| guide.var_names | guide_id | string index | Always | Canonical unique guide identifier; also retained as guide.var['guide_id']. | create_mdata.py |
| guide.var | guide_id | string | Always | Canonical guide identifier used to join metadata, assignments, and per-guide results. | Guide metadata/create_mdata.py |
| guide.var | spacer | DNA string | Always | Captured guide spacer sequence used to construct and validate the guide reference. | Guide metadata |
| guide.var | targeting | boolean | Always | True for targeting guides and false for non-targeting/control guides. | Guide metadata/create_mdata.py |
| guide.var | type | categorical | Always | Guide class such as targeting, non-targeting, safe-targeting, or negative control; filled from targeting when missing. | Guide metadata/create_mdata.py |
| guide.var | guide_chr | string/categorical | Targeting guides when known | Chromosome of the guide/protospacer genomic mapping. | Guide metadata |
| guide.var | guide_start | integer/float nullable | Targeting guides when known | Guide/protospacer genomic start coordinate. | Guide metadata |
| guide.var | guide_end | integer/float nullable | Targeting guides when known | Guide/protospacer genomic end coordinate. | Guide metadata |
| guide.var | strand | string/categorical | When supplied | Guide/protospacer strand. | Guide metadata |
| guide.var | pam | string/categorical | When supplied | Protospacer-adjacent motif associated with the guide. | Guide metadata |
| guide.var | genomic_element | string/categorical | When supplied | Source label for the targeted genomic element. | Guide metadata |
| guide.var | intended_target_name | string/categorical | Always after validation | Canonical element or target-gene group used to aggregate guides and construct inference pairs. | Guide metadata/intended_target_key_utils.py |
| guide.var | intended_target_chr | string/categorical | Targeting elements when known | Chromosome of the intended target element. | Guide metadata |
| guide.var | intended_target_start | nullable integer | Targeting elements when known | Start coordinate of the intended target element. | Guide metadata |
| guide.var | intended_target_end | nullable integer | Targeting elements when known | End coordinate of the intended target element. | Guide metadata |
| guide.var | intended_target_key | string/categorical | Always after canonicalization | Collision-safe element grouping key; controls are bucketed as non-targeting\|N and targeting guides use canonical target identity. | intended_target_key_utils.py |
| guide.var | putative_target_genes | string/categorical | When supplied | Source-provided candidate gene or genes for the targeted element. | Guide metadata passthrough |
| guide.var | reporter | boolean/numeric nullable | When supplied | Source annotation indicating reporter-associated guides/elements. | Guide metadata passthrough |
| guide.var | imperfect | boolean/numeric nullable | When supplied | Source annotation flagging imperfect guide matches/designs. | Guide metadata passthrough |
| guide.var | gene_name | string | When supplied | Source gene symbol annotation associated with the guide or target. | Guide metadata passthrough |
| guide.var | label | string | When supplied | Source display or group label. | Guide metadata passthrough |
| guide.var | description | string/categorical | When supplied | Free-text source description of the guide or element. | Guide metadata passthrough |
| guide.var | guide_class | string/categorical | When supplied | Source-specific guide classification. | Guide metadata passthrough |
| guide.var | source_group_id | string/categorical | When supplied | Original group/element identifier retained for provenance. | Guide metadata passthrough |
| guide.var | element_chr | string/categorical | When supplied | Canonical or repaired element chromosome retained separately from intended_target_chr. | Guide metadata passthrough |
| guide.var | element_start | nullable integer/string | When supplied | Canonical or repaired element start coordinate. | Guide metadata passthrough |
| guide.var | element_end | nullable integer/string | When supplied | Canonical or repaired element end coordinate. | Guide metadata passthrough |
| guide.var | general_group | string/categorical | When supplied | Source-defined broader guide/element grouping. | Guide metadata passthrough |
| guide.var | mapping_status | string/categorical | When supplied | Outcome of guide-to-genome mapping, such as unique, multi-hit, or unmapped. | Metadata preparation |
| guide.var | mapping_hit_count | integer | When supplied | Number of genomic mappings found for the spacer. | Metadata preparation |
| guide.var | mapping_scope | string/categorical | When supplied | Reference/search scope used for guide mapping. | Metadata preparation |
| guide.var | pam_model | string/categorical | When supplied | PAM rule/model used during sequence validation or mapping. | Metadata preparation |
| guide.var | spacer_encoding | string/categorical | When supplied | Description of how the stored spacer sequence was transformed or encoded. | Metadata preparation |
| guide.var | genomic_protospacer | DNA string | When supplied | Reference-matched genomic protospacer sequence. | Metadata preparation |
| guide.var | biological_spacer | DNA string | When supplied | Original biological spacer before capture-specific transformation. | Metadata preparation |
| guide.var | intended_target_gene_id | string/categorical | When supplied | Canonical gene identifier linked to the intended target. | Metadata preparation |
| guide.var | gene_annotation_chr | string/categorical | When supplied | Chromosome of the target-gene annotation used during metadata repair. | Metadata preparation |
| guide.var | gene_annotation_start | nullable numeric | When supplied | Start of the target-gene annotation used during metadata repair. | Metadata preparation |
| guide.var | gene_annotation_end | nullable numeric | When supplied | End of the target-gene annotation used during metadata repair. | Metadata preparation |
| guide.var | gene_annotation_status | string/categorical | When supplied | Status of matching the intended target to the gene annotation. | Metadata preparation |
| guide.var | coordinate_system | string/categorical | When supplied | Coordinate convention, for example zero-based half-open or one-based inclusive. | Metadata preparation |
| guide.var | reference_build | string/categorical | When supplied | Genome build used by the canonical coordinates. | Metadata preparation |
| guide.var | reference_release | string/categorical | When supplied | Genome/GTF reference release used during annotation. | Metadata preparation |
| guide.var | source | string/categorical | When supplied | Source dataset or metadata table identifier. | Guide metadata passthrough |
| guide.var | source_reference_build | string/categorical | When supplied | Genome build of the original source coordinates. | Metadata preparation |
| guide.var | source_mapping_status | string/categorical | When supplied | Original mapping status before canonical repair. | Metadata preparation |
| guide.var | source_guide_chr | string/categorical | When supplied | Original guide chromosome retained before coordinate conversion. | Metadata preparation |
| guide.var | source_guide_start | nullable numeric | When supplied | Original guide start coordinate retained before conversion. | Metadata preparation |
| guide.var | source_guide_end | nullable numeric | When supplied | Original guide end coordinate retained before conversion. | Metadata preparation |
| guide.var | source_guide_strand | string/categorical | When supplied | Original guide strand retained before conversion. | Metadata preparation |
| guide.var | guide_liftover_status | string/categorical | When supplied | Status of guide-coordinate liftover between genome builds. | Metadata preparation |
| guide.var | element_liftover_status | string/categorical | When supplied | Status of element-coordinate liftover. | Metadata preparation |
| guide.var | intended_target_liftover_status | string/categorical | When supplied | Status of intended-target-coordinate liftover. | Metadata preparation |
| guide.var | gene_annotation_liftover_status | string/categorical | When supplied | Status of gene-annotation coordinate liftover. | Metadata preparation |
| guide.var | sequence_validation_hg38 | string/categorical | When supplied | Result of validating the spacer/protospacer against hg38. | Metadata preparation |
| guide.var | coordinate_transfer_method | string/categorical | When supplied | Method used to transfer coordinates between builds. | Metadata preparation |
| guide.var | target_coordinate_repair_method | string/categorical | When supplied | Method used to repair or infer intended-target coordinates. | Metadata preparation |
| guide.var | source_guide_id | string | When supplied | Original guide identifier retained after alias consolidation. | Metadata preparation |
| guide.var | library_id | string | When supplied | Source guide-library identifier. | Guide metadata passthrough |
| guide.var | source_element_hg19 | string | When supplied | Original hg19 element label. | Metadata preparation |
| guide.var | element_hg19_chr | string | When supplied | Original hg19 element chromosome. | Metadata preparation |
| guide.var | element_hg19_start | integer | When supplied | Original hg19 element start. | Metadata preparation |
| guide.var | element_hg19_end | integer | When supplied | Original hg19 element end. | Metadata preparation |
| guide.var | element_hg38_chr | string | When supplied | Lifted/canonical hg38 element chromosome. | Metadata preparation |
| guide.var | element_hg38_start | integer | When supplied | Lifted/canonical hg38 element start. | Metadata preparation |
| guide.var | element_hg38_end | integer | When supplied | Lifted/canonical hg38 element end. | Metadata preparation |
| guide.var | alias_guide_ids | string | When supplied | Delimited guide aliases collapsed into the canonical guide. | Metadata preparation |
| guide.var | alias_source_guide_ids | string | When supplied | Delimited original-source guide identifiers represented by the canonical guide. | Metadata preparation |
| guide.var | capture_length | integer | When supplied | Length of the capture sequence used for guide mapping. | Metadata preparation |
| guide.var | capture_rule | string | When supplied | Rule used to derive the captured sequence from the biological spacer. | Metadata preparation |
| guide.var | vector_primer_motif | DNA string | When supplied | Vector/primer motif used to identify or transform the captured guide sequence. | Metadata preparation |

## Hashing fields

| Path | Field | Type | Availability | Description | Producer/source |
|---|---|---|---|---|---|
| hashing.var_names | hash_id | string index | Hashing only | Canonical hashtag/HTO feature identifier from the mapping reference. | Hash metadata/mapping |
| hashing.var | sequence | DNA string | Hashing when retained by mapping | Hashtag oligo sequence supplied in hash metadata. | Hash metadata |
| hashing.var | Additional columns | input dependent | Optional; input dependent | Hashing metadata columns carried by the mapped AnnData; the pipeline does not impose a richer normalized hashing-var schema. | Hash mapping passthrough |

## Unstructured

| Path | Field | Type | Availability | Description | Producer/source |
|---|---|---|---|---|---|
| guide.uns | capture_method | one-element string array | Always | Guide capture design, for example crop-seq or direct-capture. | create_mdata.py |
| guide.uns | moi | one-element string array | Always | Recorded or inferred multiplicity-of-infection class: high or low. | create_mdata.py |
| MuData.uns | pairs_to_test | DataFrame | Inference preparation only | Requested guide-gene pairs with guide_id and gene_id; used by SCEPTRE and native PerTurbo local/cis inference. | prepare_inference.py |
| MuData.uns | local_analysis_per_guide_results | DataFrame | Current final default inference | Merged local/cis guide-gene results and annotations. | merge_local_global_results.py |
| MuData.uns | local_analysis_per_element_results | DataFrame | Current final default inference | Merged local/cis target-element-gene results aggregated across guides. | merge_local_global_results.py |
| MuData.uns | global_analysis_per_guide_results | DataFrame | Current final default inference | Global/trans all-by-all PerTurbo guide-gene results and annotations. | merge_local_global_results.py |
| MuData.uns | global_analysis_per_element_results | DataFrame | Current final default inference | Global/trans all-by-all PerTurbo element-gene results aggregated across guides. | merge_local_global_results.py |
| MuData.uns | cis_per_guide_results | DataFrame | Legacy output alias | Older name for local per-guide results; read for backward compatibility. | legacy pipeline |
| MuData.uns | cis_per_element_results | DataFrame | Legacy output alias | Older name for local per-element results; read for backward compatibility. | legacy pipeline |
| MuData.uns | trans_per_guide_results | DataFrame | Legacy output alias | Older name for global per-guide results; read for backward compatibility. | legacy pipeline |
| MuData.uns | trans_per_element_results | DataFrame | Legacy output alias | Older name for global per-element results; read for backward compatibility. | legacy pipeline |

## Result columns

| Path | Field | Type | Availability | Description | Producer/source |
|---|---|---|---|---|---|
| per-guide results | gene_id | string | All per-guide tables | Canonical tested-gene identifier. | analysis_output_formatting.py |
| per-guide results | guide_id | string | All per-guide tables | Canonical tested-guide identifier. | analysis_output_formatting.py |
| per-element results | intended_target_name | string | All per-element tables | Canonical element/group identifier used to aggregate guides. | analysis_output_formatting.py |
| per-element results | intended_target_chr | string nullable | All per-element tables | Chromosome of the intended target element. | analysis_output_formatting.py |
| per-element results | intended_target_start | integer nullable | All per-element tables | Start coordinate of the intended target element. | analysis_output_formatting.py |
| per-element results | intended_target_end | integer nullable | All per-element tables | End coordinate of the intended target element. | analysis_output_formatting.py |
| result metrics | sceptre_log2_fc | float | Local results when SCEPTRE ran | SCEPTRE estimated log2 fold change for the perturbation test. | SCEPTRE formatter |
| result metrics | sceptre_p_value | float | Local results when SCEPTRE ran | Raw SCEPTRE p-value. | SCEPTRE formatter |
| result metrics | sceptre_q_value | float | Local results when SCEPTRE ran | Benjamini-Hochberg-adjusted SCEPTRE p-value within the exported test family. | merge_local_global_results.py |
| result metrics | sceptre_fc_se | float | Local results when available | Standard error of the SCEPTRE fold-change estimate. | SCEPTRE formatter |
| result metrics | sceptre_negLog10p | float | Local results when SCEPTRE ran | Negative base-10 logarithm of the SCEPTRE p-value with numerical clipping. | merge_local_global_results.py |
| result metrics | perturbo_log2_fc | float | When PerTurbo ran | PerTurbo estimated log2 fold change. | PerTurbo adapter |
| result metrics | perturbo_p_value | float | When PerTurbo ran | Native/empirical PerTurbo p-value exposed by the pipeline adapter. | PerTurbo adapter |
| result metrics | perturbo_q_value | float | When PerTurbo ran | PerTurbo q-value; current native output is used when provided. | PerTurbo adapter/merger |
| result metrics | perturbo_fc_se | float | When PerTurbo ran | Standard error of the PerTurbo fold-change estimate. | PerTurbo adapter |
| result metrics | perturbo_negLog10p | float | When PerTurbo ran | Negative base-10 logarithm of perturbo_p_value with numerical clipping. | merge_local_global_results.py |
| legacy result metrics | sceptre_log10_p_value | float | Legacy tables only | Older signed/legacy log10 p-value field; superseded by sceptre_negLog10p. | legacy pipeline |
| legacy result metrics | perturbo_log10_p_value | float | Legacy tables only | Older signed/legacy log10 p-value field; superseded by perturbo_negLog10p. | legacy pipeline |
| per-guide annotations | guide_sequence | DNA string | Current annotated tables | Guide spacer/capture sequence copied from guide.var['spacer']. | analysis_output_formatting.py |
| per-guide annotations | guide_type | string | Current annotated tables | Guide type copied from guide.var['type']. | analysis_output_formatting.py |
| per-guide annotations | targeting | boolean | Current annotated tables | Whether the guide is targeting rather than a control. | analysis_output_formatting.py |
| per-guide annotations | guide_chr | string nullable | Current annotated tables | Guide genomic chromosome. | analysis_output_formatting.py |
| per-guide annotations | guide_start | integer nullable | Current annotated tables | Guide genomic start coordinate. | analysis_output_formatting.py |
| per-guide annotations | guide_end | integer nullable | Current annotated tables | Guide genomic end coordinate. | analysis_output_formatting.py |
| per-guide annotations | guide_strand | string nullable | Current annotated tables | Guide strand copied from guide.var['strand']. | analysis_output_formatting.py |
| per-guide annotations | pam | string nullable | Current annotated tables | Guide PAM annotation. | analysis_output_formatting.py |
| per-guide annotations | intended_target_name | string | Current annotated tables | Canonical target element/gene group for the guide. | analysis_output_formatting.py |
| per-guide annotations | intended_target_chr | string nullable | Current annotated tables | Intended-target chromosome. | analysis_output_formatting.py |
| per-guide annotations | intended_target_start | integer nullable | Current annotated tables | Intended-target start coordinate. | analysis_output_formatting.py |
| per-guide annotations | intended_target_end | integer nullable | Current annotated tables | Intended-target end coordinate. | analysis_output_formatting.py |
| shared annotations | gene_name | string nullable | Current annotated tables | Human-readable tested-gene symbol. | analysis_output_formatting.py |
| shared annotations | nPerturbedCells | integer nullable | Current annotated tables | Number of cells assigned to the tested guide or element under guide_assignment. | analysis_output_formatting.py |
| per-element annotations | element_id | string | Current annotated element tables | Catalog-facing canonical element identifier. | build_catalog_per_element_output.py |
| per-element annotations | element_type | string | Current annotated element tables | Catalog-facing element class. | build_catalog_per_element_output.py |
| per-element annotations | element_chr | string nullable | Current annotated element tables | Catalog-facing element chromosome. | build_catalog_per_element_output.py |
| per-element annotations | element_start | integer nullable | Current annotated element tables | Catalog-facing element start coordinate. | build_catalog_per_element_output.py |
| per-element annotations | element_end | integer nullable | Current annotated element tables | Catalog-facing element end coordinate. | build_catalog_per_element_output.py |
| per-element annotations | element_name | string | Current annotated element tables | Human-readable/catalog element name. | build_catalog_per_element_output.py |
| per-element annotations | guide_ids | delimited string | Current annotated element tables | Guide identifiers aggregated into the element. | analysis_output_formatting.py |
| per-element annotations | num_guides | integer | Current annotated element tables | Number of guides represented by the element. | analysis_output_formatting.py |
