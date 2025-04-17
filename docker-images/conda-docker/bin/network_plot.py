#!/usr/bin/env python
import os
import argparse
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData
import networkx as nx
import matplotlib.pyplot as plt
from typing import Optional, Literal
from math import ceil

def plot_network(mdata: MuData, 
                central_node: str, 
                method: Literal['sceptre', 'perturbo'],
                source_column: str = "source",
                target_column: str = "target", 
                min_weight: Optional[float] = None,
                results_key: Optional[str] = "test_results", 
                ax: Optional[plt.Axes] = None):
    
    if ax is None:
        ax = plt.gca()

    # Set method-specific column names
    weight_column = f"{method}_log2_fc"
    node_size_column = f"{method}_p_value"
    
    results_df = pd.DataFrame({k: v for k, v in mdata.uns[results_key].items()})
    results_df = results_df.drop_duplicates()
    
    if min_weight is not None:
        results_df = results_df[results_df[weight_column].abs() >= min_weight]
    
    if central_node is not None:
        filtered_df = results_df[results_df[source_column] == central_node]
        
        if filtered_df.empty:
            # Create empty plot with message instead of returning False
            ax.text(0.5, 0.5, f"No network found for\n{central_node}", 
                   ha='center', va='center')
            ax.set_title(f"{method.capitalize()} Network - {central_node}")
            ax.axis('off')
            return False
        results_df = filtered_df

    # Create the network 
    G = nx.DiGraph()
    for i, row in results_df.iterrows():
        G.add_edge(row[source_column], row[target_column], 
                  weight=row[weight_column],
                  pvalue=row[node_size_column])
    
    pos = nx.circular_layout(G)

    # Calculate node sizes based on -log10(p-value)
    node_sizes = {}
    for node in G.nodes():
        edges = G.in_edges(node, data=True)
        if edges:
            min_pvalue = min(d['pvalue'] for _, _, d in edges)
            node_sizes[node] = -np.log10(min_pvalue) * 100
        else:
            node_sizes[node] = 100

    # Draw nodes
    sizes = [node_sizes.get(n, 100) for n in G.nodes()]
    nx.draw(G, pos, with_labels=True, node_color="skyblue", 
            node_size=sizes, edge_cmap=plt.cm.coolwarm, arrowsize=20, ax=ax)

    # Draw labels
    label_pos = {k: (v[0]+0.2, v[1] + 0.05) for k, v in pos.items()}
    nx.draw_networkx_labels(G, label_pos, font_size=8)
    
    # Draw edges with weights
    edge_labels = {(u, v): f"FC:{d['weight']:.2f}\np:{d['pvalue']:.2e}" 
                  for u, v, d in G.edges(data=True)}
    edge_colors = [d["weight"] for u, v, d in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, 
                          edge_cmap=plt.cm.coolwarm, arrowsize=5)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)

    ax.set_title(f"{method.capitalize()} Network - {central_node}")
    return True 

def create_method_figure(mdata: MuData,
                        central_node: str,
                        method: str,
                        source_column: str,
                        target_column: str,
                        min_weight: float,
                        results_key: str):
    """Create a figure for a specific method and central node."""
    
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)
    plot_network(mdata, central_node, method,
                source_column=source_column,
                target_column=target_column,
                min_weight=min_weight,
                results_key=results_key,
                ax=ax)
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot networks from MuData for available analysis methods")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--central_nodes", type=str, nargs='+', required=True, 
                      help="Central nodes to plot")
    parser.add_argument("--source_column", type=str, default="intended_target_name",
                      help="Source column name in test results")
    parser.add_argument("--target_column", type=str, default="gene_id",
                      help="Target column name in test results")
    parser.add_argument("--min_weight", type=float, default=0.1,
                      help="Minimum absolute log2 fold change to filter edges")
    parser.add_argument("--results_key", type=str, default="test_results",
                      help="Key for test results in mdata.uns")
    
    args = parser.parse_args()
    
    # Load MuData
    print("Loading MuData file...")
    mdata = mu.read(args.mdata_path)
    
    # Detect available methods
    results_df = pd.DataFrame({k: v for k, v in mdata.uns[args.results_key].items()})
    cols = results_df.columns
    
    available_methods = []
    if 'sceptre_log2_fc' in cols and not results_df['sceptre_log2_fc'].isna().all():
        available_methods.append('sceptre')
    if 'perturbo_log2_fc' in cols and not results_df['perturbo_log2_fc'].isna().all():
        available_methods.append('perturbo')
        
    print(f"Available methods: {available_methods}")
    
    if not available_methods:
        print("No methods with valid data found")
        exit()
    
    # Create output directory
    output_dir = "evaluation_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create plots for each method and central node
    for method in available_methods:
        print(f"\nGenerating plots for {method}...")
        
        # Calculate number of subplots needed
        n_plots = len(args.central_nodes)
        n_cols = min(3, n_plots)  # Max 3 columns
        n_rows = ceil(n_plots / n_cols)
        
        # Create figure for this method
        fig = plt.figure(figsize=(7 * n_cols, 7 * n_rows))
        fig.suptitle(f"{method.capitalize()} Networks", fontsize=16, y=0.95)
        
        # Create subplots for each central node
        for i, central_node in enumerate(args.central_nodes):
            ax = fig.add_subplot(n_rows, n_cols, i + 1)
            plot_network(mdata, central_node, method,
                        source_column=args.source_column,
                        target_column=args.target_column,
                        min_weight=args.min_weight,
                        results_key=args.results_key,
                        ax=ax)
        
        plt.tight_layout()
        
        # Save the plot
        output_file = os.path.join(output_dir, f"{method}_network_plot.png")
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {method} plots to {output_file}")

    print("All plots generated successfully")