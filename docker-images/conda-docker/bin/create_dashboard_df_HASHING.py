#!/usr/bin/env python

import os
import pandas as pd
import argparse
import muon as mu
import anndata as ad
import numpy as np
import pickle

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

def create_dashboard_df(guide_fq_tbl, hashing_fq_tbl, mudata_path, gene_ann_path, filtered_ann_path, guide_ann_path, hashing_ann_path, hashing_demux_path, hashing_unfiltered_demux_path):
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
    guide_assignment_matrix = mudata.mod['guide'].layers['guide_assignment'].toarray()
    guide_highlight = f"Number of guide barcodes (unfiltered) intersecting with scRNA barcodes (unfiltered): {human_format(intersection_guidebc_scrnabc)}, % of guides barcodes (unfiltered) intersecting with the scRNA barcode (unfiltered): {str(np.round((intersection_guidebc_scrnabc / guide_ann.obs.shape[0]) * 100, 2))}%"
    guide_img_df = new_block('Guide', '', 'Visualization', guide_highlight, True,
                    image = ['figures/guides_per_cell_histogram.png', 'figures/cells_per_guide_histogram.png', 'figures/guides_UMI_thresholds.png'],
                    image_description=['Histogram of guides per cell.', 'Histogram of cells per guide.', 'Simulating the final number of cells with assigned guides using different minimal number thresholds (at least one guide > threshold value). (Use it to inspect how many cells would have assigned guides. This can be used to check if the final number of cells with guides fit with your expected number of cells)'])

    ### Create guide inference df
    inference_table = pd.DataFrame({k: v for k, v in mudata.uns['test_results'].items()})
    targets = mudata.mod['guide'].var['intended_target_name'].values

    gi_highlight = f"Total tested sgRNA-gene pairs: {inference_table.shape[0]}"

    gi_df = new_block('Inference', '', 'Guide Inference', gi_highlight, True, inference_table,
            table_description='Inference table gene, guide, target name, lfc2, p-value, pair-type)')

    ### Create guide assignment df
    positive_calls = guide_assignment_matrix > 0
    number_of_guide_barcodes_with_positive_call = (positive_calls.sum(axis=1) > 0).sum()
    cell_ids = mudata.mod['guide'].obs.index 
    guide_ids = mudata.mod['guide'].var.index 
    df_guide_assignment = pd.DataFrame(guide_assignment_matrix, index=cell_ids, columns=guide_ids)
    sgRNA_frequencies = (df_guide_assignment > 0).sum(axis=0)
    df_sgRNA_frequencies = sgRNA_frequencies.reset_index()
    df_sgRNA_frequencies.columns = ['sgRNA', 'Frequency']
    median_frequency = df_sgRNA_frequencies['Frequency'].median()
    gs_highlight=f"Number of Cells with at least one sgRNA assigned: {human_format(number_of_guide_barcodes_with_positive_call)}, The median of cells with a positive sgRNA call is: {median_frequency}"
    
    df_sgRNA_table = df_sgRNA_frequencies.copy()
    df_sgRNA_table.columns = ['sgRNA', '# guide barcodes']
    gs_img_df = new_block('Guide', '', 'Guide Assignment', gs_highlight, True, table=df_sgRNA_table,
                    table_description='Number of guide barcodes with a positive sgRNA call',
                    image = ['figures/guides_hist_num_sgRNA.png'],
                    image_description=['Histogram of the number of sgRNA represented per cell'])
    ### Create inference df
    ##mean guides/cell
    guides_per_cell = np.sum(mudata.mod['guide'].X, axis=1)
    mean_guides_per_cell = np.mean(guides_per_cell)
    ##mean cell/guides
    cells_per_guide = np.sum(mudata.mod['guide'].X, axis=0)
    mean_cells_per_guide = np.mean(cells_per_guide)

    iv_highlight = f"Mean guides per cell: {human_format(mean_guides_per_cell)}, Mean cells per guide: {human_format(mean_cells_per_guide)}"
    
    # Collect existing network plots
    network_plots = []
    network_descs = []

    if os.path.exists('evaluation_output/sceptre_network_plot.png'):
        network_plots.append('evaluation_output/sceptre_network_plot.png')
        network_descs.append('Sceptre network plot')

    if os.path.exists('evaluation_output/perturbo_network_plot.png'):
        network_plots.append('evaluation_output/perturbo_network_plot.png')
        network_descs.append('Perturbo network plot')

    # Collect existing volcano plots
    volcano_plots = []
    volcano_descs = []

    if os.path.exists('evaluation_output/sceptre_volcano_plot.png'):
        volcano_plots.append('evaluation_output/sceptre_volcano_plot.png')
        volcano_descs.append('Sceptre volcano plot')

    if os.path.exists('evaluation_output/perturbo_volcano_plot.png'):
        volcano_plots.append('evaluation_output/perturbo_volcano_plot.png')
        volcano_descs.append('Perturbo volcano plot')

    # Combine network + volcano
    all_plots = network_plots + volcano_plots
    all_descs = network_descs + volcano_descs

    inf_img_df = new_block('Inference', '', 'Visualization', iv_highlight, True, image=all_plots, image_description=all_descs)
    
    ### Create hashing demultiplex df
    hs_highlight = f"% of cells identified as negative(no signals to any hashtags): {(hashing_unfiltered_demux.obs['hto_type_split'] == 'negative').mean() *100:.2f}, % of cells identified as singlet positives (positive signal to one hashtag) after demultiplex: {hashing_demux.shape[0]/hashing_unfiltered_demux.shape[0]*100:.2f}"
    ### adding barplot

    hs_demux_df = new_block('Hashing', '', 'Demultiplex', hs_highlight, True, 
                        image = ['figures/cells_per_hto_barplot.png', 'figures/umap_hto.png', 'figures/umap_hto_singlets.png'], 
                        image_description = ['Number of Cells across Different HTOs', 'UMAP Clustering of Cells Based on HTOs (The dimensions represent the distribution of HTOs in each cell)', 'UMAP Clustering of Cells Based on HTOs (multiplets removed)'])
    
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

    
    return guide_check_df, hashing_check_df, cell_stats, gene_stats, rna_img_df, guide_img_df, gi_df, gs_img_df, inf_img_df, hs_demux_df

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
    
    parser.add_argument('--output', type=str, default='all_df.pkl', help='Path to output pickle file')
    
    args = parser.parse_args()

    json_df = create_json_df(args.json_dir)
    ## adding plots
    guide_check_df, hashing_check_df, cell_stats, gene_stats, rna_img_df, guide_img_df, gi_df, gs_img_df, inf_img_df, hs_demux_df = create_dashboard_df(args.guide_fq_tbl, args.hashing_fq_tbl, args.mudata, args.gene_ann, args.gene_ann_filtered, args.guide_ann, args.hashing_ann, args.hashing_demux, args.hashing_unfiltered_demux)
    
    ## consider the order of modules
    json_df_sorted = json_df.sort_values(by='description', ascending=True)
    all_df = pd.concat([guide_check_df, hashing_check_df, cell_stats, gene_stats, json_df_sorted, hs_demux_df, rna_img_df, guide_img_df, inf_img_df, gs_img_df, gi_df])

    with open(args.output, 'wb') as f:
        pickle.dump(all_df, f)
    print(f"DataFrame saved to {args.output}")

if __name__ == "__main__":
    main()
