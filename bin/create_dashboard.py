#!/usr/bin/env python

import os
import re
import argparse
import pandas as pd
import numpy as np
import pickle

def generate_html_content(df, svg_mapping):
    html = '<div class="dashboard-container">'
    html += '<div class="sidebar" role="navigation" aria-label="Dashboard navigation">'

    # Create tab buttons
    tab_order = ['Filtering Summary', 'scRNA', 'Guide', 'Inference']
    modality_set = sorted(set(df['modality']), key=lambda x: tab_order.index(x))
    
    for i, modality in enumerate(modality_set):
        active_class = 'active' if i == 0 else ''
        html += f'<button class="tablinks {active_class}" onclick="openTab(event, \'{modality}\')" aria-controls="{modality}" role="tab" aria-selected="{str(i == 0).lower()}">{modality}</button>'

    html += '</div>'
    html += '<div class="content" role="main">'
    html += '''
        <div id="settings-panel" class="settings-panel">
            <h3>Dashboard Settings</h3>
            <div class="settings-options">
                <label class="setting-item">
                    <input type="checkbox" id="show-images" checked> Show Images
                </label>
                <label class="setting-item">
                    Cards per row:
                    <select id="cards-per-row">
                        <option subject="2">2</option>
                        <option subject="3" selected>3</option>
                        <option subject="4">4</option>
                    </select>
                </label>
            </div>
        </div>
        <div class="filter-sort-container">
            <input type="text" id="filter-input" placeholder="Filter by dataset name...">
            <select id="sort-select">
                <option subject="name-asc">Name (A-Z)</option>
                <option subject="name-desc">Name (Z-A)</option>
            </select>
        </div>
    '''

    # Create content for each tab
    for i, modality in enumerate(modality_set):
        rows = df[df['modality'] == modality]
        display_style = 'block' if i == 0 else 'none'
        html += f'<div id="{modality}" class="tabcontent" style="display: {display_style};" role="tabpanel" aria-labelledby="tab-{modality}">'
        html += '<div class="flex-container">'

        for j, (_, row) in enumerate(rows.iterrows()):
            highlight_class = 'highlight' if row['highlighted'] else ''
            table_html = row['table'].to_html(index=False) if not row['table'].empty else ""

            # image handling
            image_html = ""
            if isinstance(row['image'], list) and row['image']:
                image_html += '<div class="figure-container">'
                for idx, img_path in enumerate(row['image']):
                    image_name = os.path.basename(img_path)
                    svg_path = svg_mapping.get(image_name, f'svg/default_icon.svg')
                    image_description = row['image_description'][idx] if isinstance(row['image_description'], list) and len(row['image_description']) > idx else "Image description"

                    image_html += f'''
                    <div class="button-description-container">
                        <button class="toggle-button" onclick="toggleElement(this, 'figure')" data-imgsrc="{img_path}" aria-label="Show figure">
                            <img src="{svg_path}" alt="Figure {idx + 1}" class="button-icon">
                        </button>
                        <span class="description">{image_description}</span>
                    </div>
                    '''
                image_html += '</div>'
            elif isinstance(row['image'], str) and row['image']:
                # Use a generic or default SVG path for single images
                image_name = os.path.basename(row['image'])
                svg_path = svg_mapping.get(image_name, 'svg/default_icon.svg')
                image_description = row['image_description'] if isinstance(row['image_description'], str) else "Image description"

                image_html += f'''
                <div class="figure-container">
                    <div class="figure-wrapper">
                        <button class="toggle-button" onclick="toggleElement(this, 'figure')" data-imgsrc="{row['image']}">
                            <img src="{svg_path}" alt="Figure 1">
                        </button>
                        <span class="image-description">{image_description}</span>
                    </div>
                </div>
                '''
            # Highlight values in value_display if highlighted is True
            if row['highlighted']:
                # Highlight numeric values and percentages
                value_display_lines = re.sub(r'(\d+\.?\d*\w?%)', r'<span class="highlighted-value">\1</span>', row['value_display'])
                value_display_lines = re.sub(r'(\d+\.?\d*\w?)', r'<span class="highlighted-value">\1</span>', value_display_lines)
            else:
                value_display_lines = row['value_display']

            value_display_lines = value_display_lines.replace(', ', '<br>')

            table_button_html = ''
            search_bar_html = ''
            if not row['table'].empty:
                table_button_html = f'''
                <div class="button-description-container">
                    <button class="toggle-button" onclick="toggleElement(this, 'table')" aria-label="Toggle table visibility">
                        <img src="svg/table.svg" alt="Table icon" class="button-icon">
                    </button>
                    <span class="description">{row.get('table_description', '')}</span>
                </div>
                '''
                # Add search bar for "guide assignment" table
                if "guide assignment" in row['subject'].lower():
                    search_bar_html = f'''
                    <div class="table-search-container" style="display: none;">
                        <input type="text" id="table-search-{j}" class="table-search" placeholder="Search sgRNA...">
                    </div>
                    '''
                # Add search bar for "guide inference" table
                elif "guide inference" in row['subject'].lower():
                    search_bar_html = f'''
                    <div class="table-search-container" style="display: none;">
                        <input type="text" id="table-search-{j}" class="table-search" placeholder="Search by gene_id, guide_id, or intended_target_name...">
                    </div>
                    '''
                    table_html = re.sub(
                        r'<th>(log2_fc|p_value)</th>',
                        r'''<th>\1 
                        <div class="sort-icon-wrapper">
                            <span class="sort-icon" onclick="sortTableByColumn(this.closest('.table-container'), '\1', true)">▲</span>
                            <span class="sort-icon" onclick="sortTableByColumn(this.closest('.table-container'), '\1', false)">▼</span>
                        </div>
                        </th>''',
                        table_html
                    )

                    # Add the reset button after the table
                    reset_button_html = f'''
                    <div class="reset-button-container">
                        <button class="reset-button" onclick="resetTableSort(this.closest('.table-container'))">Reset Sort</button>
                    </div>
                    '''
                    # Append the reset button HTML to the table HTML
                    table_html += reset_button_html

            html += f'''
                <div class="card" role="article" aria-labelledby="card-title-{j}">
                    <h3 id="card-title-{j}">{row['subject']}</h3>
                    <p>{row['description']}</p>
                    <p class="value-display {highlight_class}" aria-live="polite">{value_display_lines}</p>
                    {table_button_html}
                    {search_bar_html}
                    <div class="table-container" style="display: none;">
                        <div class="table-content">{table_html}</div>
                        <div class="pagination-wrapper">
                            <div class="table-pagination">
                                <button class="active">1</button>
                                <button>2</button>
                                <button>3</button>
                            </div>
                        </div>
                    </div>
                    {image_html}
                </div>
                '''

        html += '</div>'
        html += f'<div id="pagination-{modality}" class="pagination"></div>'
        html += '</div>'

    html += '</div></div>'

    # Add modal HTML
    html += '''
    <div id="figureModal" class="modal">
        <div class="modal-content">
            <span class="close">&times;</span>
            <img id="modalImage" src="" alt="Figure">
        </div>
    </div>
    '''

    return html

def df_to_html_side_tabs(df, svg_mapping):
    try:
        if df.empty:
            raise ValueError("DataFrame is empty")

        required_columns = ['modality', 'description', 'subject', 'value_display', 'highlighted', 'table', 'table_description', 'image', 'image_description']
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f"Required column '{col}' is missing from the DataFrame")

        # Generate HTML content
        html_content = generate_html_content(df, svg_mapping)

        # Read CSS and JS files
        with open('css/styles.css', 'r') as css_file:
            css_content = css_file.read()
        with open('js/script.js', 'r') as js_file:
            js_content = js_file.read()

        # Combine all parts
        full_html = f'''
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Dashboard</title>
            <style>
            {css_content}
            </style>
        </head>
        <body>
            {html_content}
            <script>
            {js_content}
            </script>
        </body>
        </html>
        '''

        return full_html

    except Exception as e:
        error_html = f'''
        <html>
        <body>
            <h1>Error generating dashboard</h1>
            <p>An error occurred while generating the dashboard: {str(e)}</p>
        </body>
        </html>
        '''
        return error_html

def main():
    parser = argparse.ArgumentParser(description="Generate an HTML dashboard from a DataFrame.")
    parser.add_argument('--input', type=str, help='Paths to the input containing the dataframes.')
    parser.add_argument('--output', type=str, default='dashboard.html', help='Path to the output HTML file.')
    parser.add_argument('--css-file', type=str, default='css/styles.css', help='Path to the CSS file.')
    parser.add_argument('--js-file', type=str, default='js/script.js', help='Path to the JavaScript file.')

    args = parser.parse_args()

    with open(args.input, 'rb') as f:
        all_df = pickle.load(f)

    # SVG mapping
    svg_mapping = {
        "seqSpec_check_plots.png": "svg/line.svg",
        "knee_plot_scRNA.png": "svg/knee_plot.svg",
        "scatterplot_scrna.png": "svg/scatter.svg",
        "violinplot_scrna.png": "svg/violin.svg",
        "scRNA_barcodes_UMI_thresholds.png" : "svg/barplot.svg",
        "guides_UMI_thresholds.png": "svg/barplot.svg",
        "guides_hist_num_sgRNA.png": "svg/barplot.svg",
        "sceptre_network_plot.png": "svg/network.svg",
        "perturbo_network_plot.png": "svg/network.svg",
        "sceptre_volcano_plot.png": "svg/volcano.svg",
        "perturbo_volcano_plot.png": "svg/volcano.svg",
        "guides_per_cell_histogram.png": "svg/barplot.svg",
        "cells_per_guide_histogram.png": "svg/barplot.svg",
        "cells_per_hto_barplot.png": "svg/barplot.svg",
        "table_icon": "svg/table.svg",
        "default_icon": "svg/table.svg"
    }

    # Generate the HTML
    dashboard_html = df_to_html_side_tabs(all_df, svg_mapping)

    # Write the HTML to the output file
    with open(args.output, 'w') as f:
        f.write(dashboard_html)
        print(f"Dashboard generated successfully. Open {args.output} in your web browser.")

if __name__ == "__main__":
    main()
