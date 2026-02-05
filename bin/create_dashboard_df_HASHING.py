#!/usr/bin/env python

import os
import pandas as pd
import argparse
import muon as mu
import anndata as ad
import numpy as np
import pickle
from scipy import sparse


def human_format(num):
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])

def new_block(modality, description, subject, value_display, highlighted=False, table=None, table_description=None, image='', image_description=''):
    if table is None:
        table = pd.DataFrame()
    if table_description is None:
        table_description = '' * len(table)

    data = {
        'modality': [modality],
        'description': [description],
        'subject': [subject],
        'value_display': [value_display],
        'highlighted': [highlighted],
        'table': [table],
        'table_description': [table_description],
        'image': [image],
        'image_description': [image_description]
    }
    return pd.DataFrame(data)


def _safe_read_tsv(path):
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()


def _get_all_row(metrics_df):
    if metrics_df.empty:
        return None
    if "batch" in metrics_df.columns:
        all_rows = metrics_df[metrics_df["batch"].astype(str) == "all"]
        if not all_rows.empty:
            return all_rows.iloc[0]
    return metrics_df.iloc[0]


def _collect_images(base_dir, image_specs):
    images = []
    descs = []
    for fname, desc in image_specs:
        fpath = os.path.join(base_dir, fname)
        if os.path.exists(fpath):
            images.append(fpath)
            descs.append(desc)
    return images, descs


def collect_additional_qc_blocks(additional_qc_dir):
    blocks = []
    if not additional_qc_dir or not os.path.exists(additional_qc_dir):
        return blocks

    # Gene QC (scRNA)
    gene_dir = os.path.join(additional_qc_dir, "gene")
    gene_metrics = _safe_read_tsv(os.path.join(gene_dir, "gene_metrics.tsv"))
    gene_row = _get_all_row(gene_metrics)
    gene_highlight = ""
    if gene_row is not None:
        parts = []
        if "n_cells" in gene_row:
            parts.append(f"Cells: {human_format(gene_row['n_cells'])}")
        if "umi_median" in gene_row:
            parts.append(f"Median UMIs: {gene_row['umi_median']:.0f}")
        if "genes_median" in gene_row:
            parts.append(f"Median genes: {gene_row['genes_median']:.0f}")
        if "mito_median" in gene_row:
            parts.append(f"Median mito%: {gene_row['mito_median']:.2f}%")
        gene_highlight = ", ".join(parts)
    gene_images, gene_descs = _collect_images(
        gene_dir,
        [
            ("gene_knee_plot.png", "Knee plot of gene UMI counts."),
            ("gene_histograms.png", "Distributions of gene QC metrics (all cells)."),
            ("gene_histograms_by_batch.png", "Distributions of gene QC metrics by batch/lane."),
            ("gene_cells_per_batch.png", "Number of cells per batch/lane."),
        ],
    )
    if not gene_metrics.empty or gene_images:
        blocks.append(
            new_block(
                "scRNA",
                "Gene mapping QC",
                "Gene Mapping QC",
                gene_highlight,
                bool(gene_highlight),
                table=gene_metrics,
                table_description="Gene mapping QC metrics (overall + per batch)",
                image=gene_images,
                image_description=gene_descs,
            )
        )

    # Guide QC
    guide_dir = os.path.join(additional_qc_dir, "guide")
    guide_metrics = _safe_read_tsv(os.path.join(guide_dir, "guide_metrics.tsv"))
    guide_row = _get_all_row(guide_metrics)
    guide_highlight = ""
    if guide_row is not None:
        parts = []
        if "n_cells" in guide_row:
            parts.append(f"Cells: {human_format(guide_row['n_cells'])}")
        if "frac_cells_with_guide" in guide_row:
            parts.append(f"% cells w/ guide: {guide_row['frac_cells_with_guide']*100:.1f}%")
        if "guides_per_cell_mean" in guide_row:
            parts.append(f"Mean guides/cell: {guide_row['guides_per_cell_mean']:.2f}")
        if "guide_umi_median" in guide_row:
            parts.append(f"Median guide UMIs: {guide_row['guide_umi_median']:.0f}")
        guide_highlight = ", ".join(parts)
    guide_images, guide_descs = _collect_images(
        guide_dir,
        [
            ("guide_knee_plot.png", "Knee plot of guide UMI counts."),
            ("guide_histograms.png", "Distributions of guide QC metrics (all cells)."),
            ("guide_histograms_by_batch.png", "Distributions of guide QC metrics by batch/lane."),
        ],
    )
    if not guide_metrics.empty or guide_images:
        blocks.append(
            new_block(
                "Guide",
                "Guide mapping QC",
                "Guide Mapping QC",
                guide_highlight,
                bool(guide_highlight),
                table=guide_metrics,
                table_description="Guide mapping QC metrics (overall + per batch)",
                image=guide_images,
                image_description=guide_descs,
            )
        )

    # Intended target QC
    intended_dir = os.path.join(additional_qc_dir, "intended_target")
    intended_metrics = _safe_read_tsv(os.path.join(intended_dir, "intended_target_metrics.tsv"))
    intended_row = _get_all_row(intended_metrics)
    intended_highlight = ""
    if intended_row is not None:
        parts = []
        if "n_guides_tested" in intended_row:
            parts.append(f"Guides tested: {human_format(intended_row['n_guides_tested'])}")
        if "frac_strong_knockdowns" in intended_row:
            parts.append(f"% strong KD: {intended_row['frac_strong_knockdowns']*100:.1f}%")
        if "frac_significant" in intended_row:
            parts.append(f"% significant: {intended_row['frac_significant']*100:.1f}%")
        if "auroc" in intended_row and pd.notna(intended_row["auroc"]):
            parts.append(f"AUROC: {intended_row['auroc']:.2f}")
        intended_highlight = ", ".join(parts)
    intended_images, intended_descs = _collect_images(
        intended_dir,
        [
            ("intended_target_volcano.png", "Volcano plot of intended target effects."),
            ("intended_target_log2fc_distribution.png", "Distribution of intended target log2FC."),
            ("intended_target_roc_pr_curves.png", "ROC/PR curves for intended vs non-targeting."),
        ],
    )
    if not intended_metrics.empty or intended_images:
        blocks.append(
            new_block(
                "Inference",
                "Intended target QC",
                "Intended Target QC",
                intended_highlight,
                bool(intended_highlight),
                table=intended_metrics,
                table_description="Intended target QC metrics",
                image=intended_images,
                image_description=intended_descs,
            )
        )

    # Trans QC
    trans_dir = os.path.join(additional_qc_dir, "trans")
    trans_metrics = _safe_read_tsv(os.path.join(trans_dir, "trans_metrics.tsv"))
    trans_row = _get_all_row(trans_metrics)
    trans_highlight = ""
    if trans_row is not None:
        parts = []
        if "n_guides_tested" in trans_row:
            parts.append(f"Guides tested: {human_format(trans_row['n_guides_tested'])}")
        if "median_significant_per_guide_targeting" in trans_row:
            parts.append(f"Median sig/guide: {trans_row['median_significant_per_guide_targeting']:.1f}")
        if "total_significant_tests" in trans_row:
            parts.append(f"Total sig tests: {human_format(trans_row['total_significant_tests'])}")
        if "auroc" in trans_row and pd.notna(trans_row["auroc"]):
            parts.append(f"AUROC: {trans_row['auroc']:.2f}")
        trans_highlight = ", ".join(parts)
    trans_images, trans_descs = _collect_images(
        trans_dir,
        [
            ("trans_volcano.png", "Volcano plot of trans effects."),
            ("trans_per_guide_distribution.png", "Distribution of trans effects per guide."),
            ("trans_roc_pr_curves.png", "ROC/PR curves for validated trans links."),
        ],
    )
    if not trans_metrics.empty or trans_images:
        blocks.append(
            new_block(
                "Inference",
                "Trans-regulatory QC",
                "Trans QC",
                trans_highlight,
                bool(trans_highlight),
                table=trans_metrics,
                table_description="Trans-regulatory QC metrics",
                image=trans_images,
                image_description=trans_descs,
            )
        )

    return blocks

def create_json_df(json_dir):
    ## prepare for json files

    list_of_params = []
    json_files = [f for f in os.listdir(json_dir) if f.endswith('.json')]
    file_groups = {prefix: [] for prefix in set('-'.join(f.split('-')[:2]) for f in json_files)}

    for file_name in json_files:
        prefix = '-'.join(file_name.split('-')[:2])
        file_groups[prefix].append(file_name)

    for prefix, files in file_groups.items():
        inspect_file = next((f for f in files if 'inspect' in f), None)
        run_info_file = next((f for f in files if 'run_info' in f), None)

        if not inspect_file or not run_info_file:
            continue

        if prefix.startswith('trans-'):
            modality, subject = ('scRNA', 'Mapping scRNA')
        elif prefix.startswith('hashing-'):
            modality, subject = ('Hashing', 'Mapping Hashing')
        else:
            modality, subject = ('Guide', 'Mapping Guide')

        description = prefix.split('-', 1)[1]

        combined_data = pd.concat([
            pd.read_json(os.path.join(json_dir, inspect_file), typ="series"),
            pd.read_json(os.path.join(json_dir, run_info_file), typ="series")
        ])
        table = combined_data.to_frame().reset_index().rename(columns={'index': 'parameter', 0: 'value'})

        # Exclude specific parameters
        table = table[~table['parameter'].isin(['start_time', 'call'])]

        # Human format certain parameters
        for param in ['n_targets','n_processed', 'n_unique', 'n_pseudoaligned',
                        'numRecords', 'numReads', 'numBarcodes', 'numUMIs', 'numBarcodeUMIs', 'gtRecords', 'numBarcodesOnOnlist', 'numReadsOnOnlist']:
            bm = table['parameter'] == param
            table.loc[bm, 'value'] = table.loc[bm, 'value'].astype(float).apply(human_format)

        # Round certain parameters
        for param in ['meanReadsPerBarcode', 'meanUMIsPerBarcode']:
            bm = table['parameter'] == param
            table.loc[bm, 'value'] = table.loc[bm, 'value'].astype(float).apply(lambda x: f"{x:.1f}")

        # Percentage conversion
        for param in ['p_pseudoaligned', 'p_unique',
                    'percentageBarcodesOnOnlist', 'percentageReadsOnOnlist']:
            bm = table['parameter'] == param
            table.loc[bm, 'value'] = table.loc[bm, 'value'].astype(float).apply(lambda x: f"{x:.1f}%")

        # Generate value_display
        n_processed = table.loc[table['parameter'] == 'n_processed', 'value'].values[0]
        n_mapped = table.loc[table['parameter'] == 'n_pseudoaligned', 'value'].values[0]
        p_mapped = table.loc[table['parameter'] == 'p_pseudoaligned', 'value'].values[0]
        numBarcodes = table.loc[table['parameter'] == 'numBarcodes', 'value'].values[0]
        value_display = f"Total Reads for {modality}: {n_processed}, Paired Reads Mapped: {n_mapped}, Alignment Percentage: {p_mapped}, Total Detected {modality} Barcodes (Unfiltered): {numBarcodes}"

        # Set highlighted to True if value_display is not empty
        highlighted = bool(value_display.strip())

        ### Create json_df
        list_of_params.append(new_block(modality, description, subject, value_display, highlighted=highlighted, table=table, table_description='Mapping summary'))
        json_df = pd.concat(list_of_params, ignore_index=True)
    return json_df


#  base_mdata.uns['per_guide_results'] = merged_guide_df
#     base_mdata.uns['per_element_results'] = merged_element_df


def create_inference_blocks(mudata, use_default=False):
    """Create inference blocks for either single test_results or separate cis/trans results"""
    inference_blocks = []

    if use_default:
        # Process both cis and trans results
        for analysis_type in ['cis', 'trans']:
            results_key = f"{analysis_type}_per_guide_results"

            if results_key in mudata.uns:
                inference_table = pd.DataFrame({k: v for k, v in mudata.uns[results_key].items()})
                n=10000

                if analysis_type == 'cis':
                    inference_table =   inference_table.sort_values(by='sceptre_p_value').head(n)
                else:
                    inference_table =   inference_table.sort_values(by='p_value').head(n)



                gi_highlight = f"Top lowest pvalues  {str(n)} tested sgRNA-gene pairs ({analysis_type}): {inference_table.shape[0]}"

                gi_df = new_block('Inference', f'{analysis_type.capitalize()} Analysis', 'Guide Inference', gi_highlight, True, inference_table,
                        table_description=f'{analysis_type.capitalize()} inference table (gene, guide, target name, lfc2, p-value, pair-type)')
                inference_blocks.append(gi_df)
            else:
                print(f"Warning: {results_key} not found in mudata.uns")
    else:
        # Process single test_results
        if 'cis_per_guide_results' in mudata.uns:
            inference_table = pd.DataFrame({k: v for k, v in mudata.uns['cis_per_guide_results'].items()})
            
            gi_highlight = f"Total tested sgRNA-gene pairs: {inference_table.shape[0]}"

            gi_df = new_block('Inference', '', 'Guide Inference', gi_highlight, True, inference_table,
                    table_description='Inference table (gene, guide, target name, lfc2, p-value, pair-type)')
            inference_blocks.append(gi_df)
        else:
            print("Warning: cis_per_guide_results not found in mudata.uns")

    return inference_blocks

def collect_evaluation_plots(use_default=False):
    """Collect evaluation plots based on whether default mode is used"""
    network_plots = []
    network_descs = []
    volcano_plots = []
    volcano_descs = []
    precission_plots = []
    precission_desc = []
    bar_plot_direct_x_control_plots = []
    bar_plot_direct_x_control_desc = []

    

    if use_default:
        # Look for cis and trans specific plots
        plot_configs = [
            ('cis_sceptre', 'Cis Sceptre'),
            ('cis_perturbo', 'Cis Perturbo'),
            ('trans_sceptre', 'Trans Sceptre'),
            ('trans_perturbo', 'Trans Perturbo')
        ]

        for plot_prefix, plot_desc in plot_configs:
            # Network plots
            network_path = f'evaluation_output/{plot_prefix}_network_plot.png'
            if os.path.exists(network_path):
                network_plots.append(network_path)
                network_descs.append(f'{plot_desc} network plot')



            #original names plot names barplot_direct_vs_control.png , precision_recall_roc.png, volcano_plot.png
            # Volcano plots
            volcano_path = f'evaluation_output/{plot_prefix}_volcano_plot.png'
            if os.path.exists(volcano_path):
                volcano_plots.append(volcano_path)
                volcano_descs.append(f'{plot_desc} direct targets x controls volcano plot')

            precission_recall_path = f'evaluation_output/{plot_prefix}_precision_recall_roc.png'
            if os.path.exists(precission_recall_path):
                precission_plots.append(precission_recall_path)
                precission_desc.append(f'{plot_desc} precision recall roc plot of direct targets x controls')



            bar_plot_direct_x_control_path = f'evaluation_output/{plot_prefix}_barplot_direct_vs_control.png'
            if os.path.exists(bar_plot_direct_x_control_path):
                bar_plot_direct_x_control_plots.append(bar_plot_direct_x_control_path)
                bar_plot_direct_x_control_desc.append(f'{plot_desc} bar plot direct targets vs control plot')
            
    else:
        # Look for standard plots
        plot_configs = [
            ('sceptre', 'Sceptre'),
            ('perturbo', 'Perturbo')
        ]

        for plot_prefix, plot_desc in plot_configs:
            # Network plots
            network_path = f'evaluation_output/{plot_prefix}_network_plot.png'
            if os.path.exists(network_path):
                network_plots.append(network_path)
                network_descs.append(f'{plot_desc} network plot')

            # Volcano plots
            volcano_path = f'evaluation_output/{plot_prefix}_volcano_plot.png'
            if os.path.exists(volcano_path):
                volcano_plots.append(volcano_path)
                volcano_descs.append(f'{plot_desc} volcano plot')
            # Volcano plots

            

    return network_plots, network_descs, volcano_plots, volcano_descs,precission_plots ,precission_desc, bar_plot_direct_x_control_plots ,bar_plot_direct_x_control_desc

def create_dashboard_df(guide_fq_tbl, hashing_fq_tbl, mudata_path, gene_ann_path, filtered_ann_path, guide_ann_path, hashing_ann_path, hashing_demux_path, hashing_unfiltered_demux_path, additional_qc_dir=None, use_default=False):
    ### Create df for cell statistics
    guide_fq_table = pd.read_csv(guide_fq_tbl)
    hashing_fq_table = pd.read_csv(hashing_fq_tbl)
    mudata = mu.read(mudata_path)
    guide_ann = ad.read_h5ad(guide_ann_path)
    gene_ann = ad.read_h5ad(gene_ann_path)
    gene_filtered_ann =ad.read_h5ad(filtered_ann_path)
    hashing_ann = ad.read_h5ad(hashing_ann_path)
    hashing_demux = ad.read_h5ad(hashing_demux_path)
    hashing_unfiltered_demux = ad.read_h5ad(hashing_unfiltered_demux_path)

    intersection_guides_and_scrna_unfitered = set(gene_ann.obs.index).intersection(guide_ann.obs.index)
    intersection_guidebc_scrnabc = len(intersection_guides_and_scrna_unfitered)

    intersection_guides_and_scrna_and_hashing_unfitered = set(gene_ann.obs.index).intersection(guide_ann.obs.index).intersection(hashing_ann.obs.index)
    intersection_guidebc_scrnabc_hashingbc = len(intersection_guides_and_scrna_and_hashing_unfitered)

    cn_highlight=f"Number of guide barcodes (unfiltered) intersecting with scRNA barcodes (unfiltered): {human_format(intersection_guidebc_scrnabc)}, Number of guide barcodes(unfiltered) intersecting with both scRNA and HTO barcodes(unfiltered): {human_format(intersection_guidebc_scrnabc_hashingbc)}, Number of cells after filtering by the minimal number of genes to consider a cell usable: {human_format(gene_filtered_ann.shape[0])},  Number of cells after filtering negative HTOs: {human_format(sum(hashing_unfiltered_demux.obs['hto_type_split'] != 'negative'))}, Number of cells after filtering negative and multiplet HTOs: {human_format(mudata.shape[0])}"

    gn_highlight=f"Number of genes detected after filtering: {human_format(mudata.mod['gene'].var.shape[0])}, Mean UMI counts per cell after filtering: {human_format(mudata.mod['gene'].X.sum(axis=1).mean())}"
    cell_stats = new_block('Filtering Summary', '', 'Filter to select high quality cells', cn_highlight, True)
    gene_stats = new_block('Filtering Summary', '', 'Gene Statistics', gn_highlight, True)

    ### Create image_df for scRNA preprocessing
    rna_img_df = new_block('scRNA', 'scRNA preprocessing', 'Visualization','', False,
        image = ['figures/knee_plot_scRNA.png', 'figures/scatterplot_scrna.png', 'figures/violinplot_scrna.png', 'figures/scRNA_barcodes_UMI_thresholds.png'],
        image_description= ['Knee plot of UMI counts vs. barcode index.', 'Scatterplot of total counts vs. genes detected, colored by mitochondrial content.','Distribution of gene counts, total counts, and mitochondrial content.', 'Number of scRNA barcodes using different\nTotal UMI thresholds.'])

    ### Create image_df for guide
    guide_assignment_matrix = mudata.mod['guide'].layers['guide_assignment']
    guide_highlight = f"Number of guide barcodes (unfiltered) intersecting with scRNA barcodes (unfiltered): {human_format(intersection_guidebc_scrnabc)}, % of guides barcodes (unfiltered) intersecting with the scRNA barcode (unfiltered): {str(np.round((intersection_guidebc_scrnabc / guide_ann.obs.shape[0]) * 100, 2))}%"
    guide_img_df = new_block('Guide', '', 'Visualization', guide_highlight, True,
                    image = ['figures/guides_per_cell_histogram.png', 'figures/cells_per_guide_histogram.png', 'figures/guides_UMI_thresholds.png'],
                    image_description=['Histogram of guides per cell.', 'Histogram of cells per guide.', 'Simulating the final number of cells with assigned guides using different minimal number thresholds (at least one guide > threshold value). (Use it to inspect how many cells would have assigned guides. This can be used to check if the final number of cells with guides fit with your expected number of cells)'])

    ### Create guide inference blocks (handles both single and cis/trans)
    inference_blocks = create_inference_blocks(mudata, use_default)

    ### Create guide assignment df
    cell_ids = mudata.mod['guide'].obs.index
    guide_ids = mudata.mod['guide'].var.index

    if sparse.issparse(mudata.mod['guide'].layers['guide_assignment']):
        df_guide_assignment = pd.DataFrame.sparse.from_spmatrix(guide_assignment_matrix, index=cell_ids, columns=guide_ids)
    else:
        df_guide_assignment = pd.DataFrame(guide_assignment_matrix, index=cell_ids, columns=guide_ids)

    sgRNA_frequencies = df_guide_assignment.sum(axis=0)
    df_sgRNA_frequencies = sgRNA_frequencies.reset_index()
    df_sgRNA_frequencies.columns = ['sgRNA', 'Frequency']
    freq_series = df_sgRNA_frequencies['Frequency']
    if isinstance(freq_series.dtype, pd.SparseDtype):
        freq_series = freq_series.sparse.to_dense()
    df_sgRNA_frequencies['Frequency'] = pd.to_numeric(freq_series, errors='coerce').fillna(0.0)
    total_sgrna_assignments = df_sgRNA_frequencies['Frequency'].sum()
    median_sgrna_frequency = df_sgRNA_frequencies['Frequency'].median()
    gs_highlight=f"Total sgRNA assignment values across all guides: {human_format(total_sgrna_assignments)}, Median assignment value per sgRNA: {human_format(median_sgrna_frequency)}"

    df_sgRNA_table = df_sgRNA_frequencies.copy()
    df_sgRNA_table.columns = ['sgRNA', 'Assignment sum']
    gs_img_df = new_block('Guide', '', 'Guide Assignment', gs_highlight, True, table=df_sgRNA_table,
                    table_description='Sum of sgRNA assignment values per sgRNA',
                    image = ['figures/guides_hist_num_sgRNA.png'],
                    image_description=['Histogram of summed sgRNA assignment values per sgRNA'])

    ### Create inference visualization df
    ##mean guides/cell
    guides_per_cell = np.sum(mudata.mod['guide'].X, axis=1)
    mean_guides_per_cell = np.mean(guides_per_cell)
    ##mean cell/guides
    cells_per_guide = np.sum(mudata.mod['guide'].X, axis=0)
    mean_cells_per_guide = np.mean(cells_per_guide)

    iv_highlight = f"Mean guides per cell: {human_format(mean_guides_per_cell)}, Mean cells per guide: {human_format(mean_cells_per_guide)}"

    # Collect evaluation plots based on use_default flag
    (
        network_plots, 
        network_descs, 
        volcano_plots, 
        volcano_descs, 
        precission_plots,
        precission_desc, 
        bar_plot_direct_x_control_plots,
        bar_plot_direct_x_control_desc
    ) = collect_evaluation_plots(use_default)
    # Combine network + volcano
    all_plots = network_plots + volcano_plots + precission_plots + bar_plot_direct_x_control_plots
    all_descs = network_descs + volcano_descs + precission_desc + bar_plot_direct_x_control_desc


    inf_img_df = new_block('Inference', '', 'Visualization', iv_highlight, True, image=all_plots, image_description=all_descs)

    ### Create hashing demultiplex df
    hs_highlight = f"% of cells identified as negative(no signals to any hashtags): {(hashing_unfiltered_demux.obs['hto_type_split'] == 'negative').mean() *100:.2f}, % of cells identified as singlet positives (positive signal to one hashtag) after demultiplex: {hashing_demux.shape[0]/hashing_unfiltered_demux.shape[0]*100:.2f}"
    ### adding barplot

    hs_demux_df = new_block('Hashing', '', 'Demultiplex', hs_highlight, True,
                        image = ['figures/cells_per_hto_barplot.png', 'figures/umap_hto.png', 'figures/umap_hto_singlets.png'],
                        image_description = ['Number of Cells across Different HTOs', 'UMAP Clustering of Cells Based on HTOs (The dimensions represent the distribution of HTOs in each cell)', 'UMAP Clustering of Cells Based on HTOs (multiplets removed)'])

    ### Collect Additional QC blocks (optional)
    qc_blocks = collect_additional_qc_blocks(additional_qc_dir)

    ### check guide seqspec check df
    guide_check_df = new_block("Guide", '', 'Fastq Overview', '', False, table = guide_fq_table,
                        table_description='Summary of Sequence Index: A summary of the positions where the Guide starts are mapped on the reads (Use to inspect or calibrate the position where the guide is supposed to be found in your SeqSpec File)',
                        image = ['guide_seqSpec_plots/seqSpec_check_plots.png'],
                        image_description= ['The frequency of each nucleotides along the Read 1 (Use to inspect the expected read parts with their expected signature) and Read 2 (Use to inspect the expected read parts with their expected signature)'])

    ### check hashing seqspec check df
    hashing_check_df = new_block("Hashing", '', 'Fastq Overview', '', False, table = hashing_fq_table,
                        table_description='Summary of Sequence Index: A summary of the positions where the Hashtag starts are mapped on the reads (Use to inspect or calibrate the position where the hashtag is supposed to be found in your SeqSpec File)',
                        image = ['hashing_seqSpec_plots/seqSpec_check_plots.png'],
                        image_description= ['The frequency of each nucleotides along the Read 1 (Use to inspect the expected read parts with their expected signature )and Read 2 (Use to inspect the expected read parts with their expected signature)'])

    return guide_check_df, hashing_check_df, cell_stats, gene_stats, rna_img_df, guide_img_df, inference_blocks, gs_img_df, inf_img_df, hs_demux_df, qc_blocks

def main():
    parser = argparse.ArgumentParser(description="Process JSON files and generate dashboard dataframes.")
    parser.add_argument('--json_dir', type=str, help="Directory containing JSON files to process")
    parser.add_argument('--guide_fq_tbl', required=True, help='Path to the guide fastq position table')
    parser.add_argument('--hashing_fq_tbl', required=True, help='Path to the hashing fastq position table')
    parser.add_argument('--mudata', required=True, help='Path to the mudata object')
    parser.add_argument('--gene_ann', required=True, help='Path to the gene anndata file')
    parser.add_argument('--gene_ann_filtered', required=True, help='Path to the gene filtered anndata file')
    parser.add_argument('--guide_ann', required=True, help='Path to the guide anndata file')
    parser.add_argument('--hashing_ann', required=True, help='Path to the hashing anndata file')
    parser.add_argument('--hashing_demux', required=True, help='Path to the hashing demux anndata file')
    parser.add_argument('--hashing_unfiltered_demux', required=True, help='Path to the hashing unfiltered demux anndata file')
    parser.add_argument('--additional_qc_dir', default=None, help='Path to Additional QC output directory')
    parser.add_argument('--default', action="store_true",
                      help="Process mudata with cis_per_guide_results and trans_per_guide_results instead of single test_results")
    parser.add_argument('--output', type=str, default='all_df.pkl', help='Path to output pickle file')

    args = parser.parse_args()

    json_df = create_json_df(args.json_dir)
    ## adding plots
    guide_check_df, hashing_check_df, cell_stats, gene_stats, rna_img_df, guide_img_df, inference_blocks, gs_img_df, inf_img_df, hs_demux_df, qc_blocks = create_dashboard_df(
        args.guide_fq_tbl,
        args.hashing_fq_tbl,
        args.mudata,
        args.gene_ann,
        args.gene_ann_filtered,
        args.guide_ann,
        args.hashing_ann,
        args.hashing_demux,
        args.hashing_unfiltered_demux,
        args.additional_qc_dir,
        args.default
    )

    ## consider the order of modules
    json_df_sorted = json_df.sort_values(by='description', ascending=True)

    # Combine all dataframes, with inference_blocks being a list now
    df_list = [guide_check_df, hashing_check_df, cell_stats, gene_stats, json_df_sorted, hs_demux_df, rna_img_df, guide_img_df, inf_img_df, gs_img_df]

    # Add QC blocks (if any)
    df_list.extend(qc_blocks)

    # Add inference blocks to the list
    df_list.extend(inference_blocks)

    all_df = pd.concat(df_list, ignore_index=True)

    with open(args.output, 'wb') as f:
        pickle.dump(all_df, f)
    print(f"DataFrame saved to {args.output}")

if __name__ == "__main__":
    main()
