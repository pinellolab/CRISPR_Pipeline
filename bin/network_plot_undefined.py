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

def plot_network(mdata: MuData, central_node: str, method: Optional[Literal['sceptre', 'perturbo']],
               source_column: str = "source", target_column: str = "target",
               min_weight: Optional[float] = None, results_key: Optional[str] = "test_results",
               ax: Optional[plt.Axes] = None):
   
   if ax is None:
       fig, ax = plt.subplots(1, 1, figsize=(15, 15))
   
   # Set column names based on whether method is specified
   if not method:  # Generic columns
       weight_column = "log2_fc"
       node_size_column = "p_value"
   else:  # Method-specific columns
       weight_column = f"{method}_log2_fc"
       node_size_column = f"{method}_p_value"
   
   results_df = pd.DataFrame({k: v for k, v in mdata.uns[results_key].items()})
   results_df = results_df.drop_duplicates()
   if min_weight is not None:
       results_df = results_df[results_df[weight_column].abs() >= min_weight]

   # Filter to rows related to the selected central node
   if central_node is not None:
       results_df = results_df[results_df[source_column] == central_node]
       
   if results_df.empty:
       print(f"No data found for {central_node} with method {method}")
       return False

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
           node_size=sizes, edge_cmap=plt.cm.coolwarm, arrowsize=20)

   label_pos = {k: (v[0]+0.2, v[1] + 0.05) for k, v in pos.items()}
   nx.draw_networkx_labels(G, label_pos, font_size=8)
   
   edge_labels = {(u, v): f"FC:{d['weight']:.2f}\np:{d['pvalue']:.2e}"
                 for u, v, d in G.edges(data=True)}
   edge_colors = [d["weight"] for u, v, d in G.edges(data=True)]
   nx.draw_networkx_edges(G, pos, edge_color=edge_colors,
                         edge_cmap=plt.cm.coolwarm, arrowsize=5)
   nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)
   
   method_title = method.capitalize() if method else "Differential Expression"
   ax.set_title(f"{method_title} Network - {central_node}")
   return True

def select_central_nodes(mdata: MuData, num_nodes: int, source_column: str,
                       target_column: str, method: Optional[Literal['sceptre', 'perturbo']],
                       min_weight: Optional[float] = None):  # Added min_weight parameter
    # Set column name based on whether method is specified
    if not method:
        weight_column = "log2_fc"
    else:
        weight_column = f"{method}_log2_fc"
    
    results_df = pd.DataFrame({k: v for k, v in mdata.uns["test_results"].items()})
    results_df = results_df.drop_duplicates()
    
    # Apply min_weight filter if specified
    if min_weight is not None:
        results_df = results_df[results_df[weight_column].abs() >= min_weight]
    
    # Create a graph and calculate the degree of each node
    G = nx.DiGraph()
    for i, row in results_df.iterrows():
        if pd.notna(row.get(weight_column)):  # Only add edges if weight exists
            G.add_edge(row[source_column], row[target_column], weight=abs(row[weight_column]))
    
    # Sort nodes by weighted degree
    degrees = dict(G.degree(weight='weight'))
    sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)

    intended_target_names = set(results_df['intended_target_name'])
    filtered_nodes = [node for node in sorted_nodes if node in intended_target_names]
    
    # Verify each node has data
    valid_nodes = []
    for node in filtered_nodes:
        node_data = results_df[results_df[source_column] == node]
        if not node_data.empty:
            valid_nodes.append(node)
            if len(valid_nodes) == num_nodes:
                break
    
    method_name = method if method else "Analysis"
    if valid_nodes:
        print(f"Top {len(valid_nodes)} nodes for {method_name}:", valid_nodes)
    else:
        print(f"No valid nodes found for {method_name} with the current filters")
    
    return valid_nodes

def main():
   print("Starting program...")
   parser = argparse.ArgumentParser(description="Select nodes and plot networks comparing two analysis methods")
   parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
   parser.add_argument("--num_nodes", type=int, required=True, help="Number of central nodes to select")
   parser.add_argument("--source_column", type=str, default="intended_target_name",
                     help="Source column name in test results")
   parser.add_argument("--target_column", type=str, default="gene_id",
                     help="Target column name in test results")
   parser.add_argument("--min_weight", type=float, default=0.1,
                     help="Minimum absolute log2 fold change to filter edges")
   parser.add_argument("--results_key", type=str, default="test_results",
                     help="Key for test results in mdata.uns")
   
   args = parser.parse_args()
   
   print("Loading MuData file...")
   mdata = mu.read(args.mdata_path)
   
   # Check available methods
   results_df = pd.DataFrame({k: v for k, v in mdata.uns[args.results_key].items()})
   cols = results_df.columns
   print("Available columns:", cols)
   
   output_dir = "evaluation_output"
   os.makedirs(output_dir, exist_ok=True)
   
   # Check for generic columns first
   if 'log2_fc' in cols and 'p_value' in cols:
       print("Using generic log2_fc and p_value columns")
       central_nodes = select_central_nodes(
           mdata, args.num_nodes, args.source_column, args.target_column, None
       )
       
       # Create one figure with subplots for all nodes
       num_rows = (len(central_nodes) + 1) // 2
       fig = plt.figure(figsize=(30, 15 * num_rows))
       
       for i, central_node in enumerate(central_nodes):
           ax = fig.add_subplot(num_rows, 2, i + 1)
           plot_network(
               mdata,
               central_node=central_node,
               method=None,
               source_column=args.source_column,
               target_column=args.target_column,
               min_weight=args.min_weight,
               results_key=args.results_key,
               ax=ax
           )
       
       plt.tight_layout()
       output_file = os.path.join(output_dir, "network_plot.png")
       plt.savefig(output_file, dpi=300, bbox_inches='tight')
       plt.close()
       print(f"Saved combined plot to {output_file}")
           
   else:
       # Check for method-specific columns
       available_methods = []
       if 'sceptre_log2_fc' in cols and not results_df['sceptre_log2_fc'].isna().all():
           available_methods.append('sceptre')
       if 'perturbo_log2_fc' in cols and not results_df['perturbo_log2_fc'].isna().all():
           available_methods.append('perturbo')
           
       print(f"Available methods: {available_methods}")
       
       if not available_methods:
           print("No methods with valid data found")
           return
       
       # Create separate plots for each method
       for method in available_methods:
           print(f"\nGenerating plot for {method}...")
           
           method_nodes = select_central_nodes(mdata, args.num_nodes, args.source_column, args.target_column, method, args.min_weight)
           
           fig = plt.figure(figsize=(30, 15))
           
           if method_nodes:
               num_rows = (len(method_nodes) + 1) // 2
               for i, central_node in enumerate(method_nodes):
                   ax = fig.add_subplot(num_rows, 2, i + 1)
                   plot_network(
                       mdata,
                       central_node=central_node,
                       method=method,
                       source_column=args.source_column,
                       target_column=args.target_column,
                       min_weight=args.min_weight,
                       results_key=args.results_key,
                       ax=ax
                   )
           else:
               # Create a single blank subplot with a message
               ax = fig.add_subplot(1, 1, 1)
               ax.text(0.5, 0.5, f"No nodes available for {method}\nwith min_weight={args.min_weight}",
                      horizontalalignment='center',
                      verticalalignment='center',
                      transform=ax.transAxes,
                      fontsize=14)
               ax.set_axis_off()
           
           plt.tight_layout()
           output_file = os.path.join(output_dir, f"{method}_network_plot.png")
           plt.savefig(output_file, dpi=300, bbox_inches='tight')
           plt.close()
           print(f"Saved {method} plots to {output_file}")

   print("All plots generated successfully")

if __name__ == "__main__":
   main()