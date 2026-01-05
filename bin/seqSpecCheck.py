#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import gzip
from collections import Counter
import numpy as np  
import os

# Set style
plt.style.use('ggplot')

def is_gzipped(filename):
    with open(filename, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def open_file(filename, mode='rt'):
    if is_gzipped(filename):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def calculate_gini(array):
    array = np.sort(array)
    index = np.arange(1, array.shape[0] + 1)
    n = array.shape[0]
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) if np.sum(array) > 0 else 0

def readFastq(filename, max_reads):
    sequences = []
    if filename.endswith('.gz'):
        filetype =  filename.split('.')[-2]
    else:
        filetype = filename.split('.')[-1]

    with open_file(filename) as fh:
        if filetype == 'fastq':
            while len(sequences) < max_reads:
                header = fh.readline()
                if not header: break
                seq = fh.readline().rstrip()
                fh.readline()
                fh.readline()
                sequences.append(seq)
        elif filetype == 'fasta':
            seq = ''
            for line in fh:
                line = line.strip()
                if line.startswith('>'):
                    if seq:
                        sequences.append(seq)
                        if len(sequences) >= max_reads: break
                        seq = ''
                else:
                    seq += line
            if seq and len(sequences) < max_reads:
                sequences.append(seq)
        else: # Fallback
            try:
                fh.seek(0)
                if fh.readline().startswith('>'):
                    fh.seek(0)
                    seq = ''
                    for line in fh:
                        line = line.strip()
                        if line.startswith('>'):
                            if seq:
                                sequences.append(seq)
                                if len(sequences) >= max_reads: break
                                seq = ''
                        else:
                            seq += line
                    if seq and len(sequences) < max_reads:
                        sequences.append(seq)
                else:
                    fh.seek(0)
                    while len(sequences) < max_reads:
                        header = fh.readline()
                        if not header: break
                        seq = fh.readline().rstrip()
                        fh.readline()
                        fh.readline()
                        sequences.append(seq)
            except:
                 raise ValueError("Unsupported filetype")
    return sequences

def fastq_sequence_plot(seqs, file_name, ax):
    if not seqs: return
    df_read = pd.DataFrame([list(n) for n in seqs]).fillna('N')
    freq = []
    for n in ['A', 'C', 'G', 'T']:
        freq.append((df_read == n).sum() / len(df_read))
    freq_df = pd.DataFrame(freq, index=['A', 'C', 'G', 'T']).T
    freq_df.plot(color=['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00'], ax=ax, alpha=0.8, linewidth=2)
    ax.set_title(f'Nuc Freq: {file_name}', fontsize=10, fontweight='bold')
    ax.set_ylabel('Freq', fontsize=8)
    ax.tick_params(axis='both', labelsize=8)

def analyze_guides_in_reads(reads, guide_list):
    positions = []
    upstream_map = {}
    guide_hits = Counter({g: 0 for g in guide_list})
    
    if not reads or not guide_list:
        return positions, upstream_map, guide_hits

    for seq in reads:
        for guide in guide_list:
            idx = seq.find(guide)
            if idx != -1:
                positions.append(idx)
                guide_hits[guide] += 1
                start_cut = max(0, idx - 12)
                upstream_fragment = seq[start_cut:idx]
                if len(upstream_fragment) < 12:
                    padding = '-' * (12 - len(upstream_fragment))
                    upstream_fragment = padding + upstream_fragment
                if idx not in upstream_map: upstream_map[idx] = []
                upstream_map[idx].append(upstream_fragment)
                break 
    return positions, upstream_map, guide_hits

def plot_positions(positions, title, ax, highlight=False):
    if not positions:
        ax.text(0.5, 0.5, "No Guides", ha='center', va='center')
        ax.set_title(title, fontsize=10)
        return
    counts = Counter(positions)
    df = pd.DataFrame.from_dict(counts, orient='index', columns=['Count']).sort_index()
    ax.bar(df.index, df['Count'], color='#984ea3', alpha=0.7)
    
    t_str = title
    if highlight:
        t_str = "★ " + title + " ★"
        ax.set_facecolor('#e6f5c9') # Light green background
        
    ax.set_title(t_str, fontsize=10, fontweight='bold' if highlight else 'normal')
    ax.set_xlabel('Index', fontsize=8)
    ax.tick_params(axis='both', labelsize=8)

def plot_top_upstream(upstream_map, positions, title, ax, highlight=False):
    if not positions:
        ax.axis('off')
        return
    best_pos = Counter(positions).most_common(1)[0][0]
    seqs_at_best = upstream_map.get(best_pos, [])
    if not seqs_at_best:
        ax.text(0.5, 0.5, "No Seq", ha='center')
        return
    top_10 = Counter(seqs_at_best).most_common(10)
    seqs, counts = zip(*top_10)
    y_pos = np.arange(len(seqs))
    colors = plt.cm.Spectral(np.linspace(0, 1, len(seqs)))
    ax.barh(y_pos, counts, color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(seqs, fontdict={'family': 'monospace', 'size': 8})
    ax.invert_yaxis() 
    
    if highlight:
        ax.set_facecolor('#e6f5c9')
        
    ax.set_title(f"{title}\n(at Pos {best_pos})", fontsize=9)
    ax.tick_params(axis='x', labelsize=7)

def plot_guide_diversity(guide_hits, title, ax, highlight=False):
    counts = np.array(list(guide_hits.values()))
    total_hits = np.sum(counts)
    if total_hits == 0:
        ax.text(0.5, 0.5, "No Hits", ha='center')
        ax.set_title(title, fontsize=10)
        return

    gini = calculate_gini(counts)
    sorted_counts = np.sort(counts)[::-1]
    
    color = '#377eb8'
    if gini > 0.6: color = '#e41a1c'
    
    ax.bar(range(len(sorted_counts)), sorted_counts, color=color, width=1.0)
    
    title_str = f"{title}\nGini: {gini:.2f}"
    if highlight:
        ax.set_facecolor('#e6f5c9')
        title_str = "★ " + title_str
        
    ax.set_title(title_str, fontsize=9, color='black' if gini <= 0.6 else 'red')
    ax.set_xlabel('Rank', fontsize=8)
    
    zeros = np.sum(counts == 0)
    ax.text(0.95, 0.95, f"Zeros: {zeros}", transform=ax.transAxes, ha='right', va='top', fontsize=8)

def calculate_score(total_reads, positions, upstream_map, guide_hits):
    """
    Equation: 
    Score = (3 * HitRatio) + (2 * PosPurity) + (1 * FlankPurity) + (1 * (1-Gini))
    """
    total_hits = len(positions)
    if total_hits == 0: return 0.0, 0.0, 0.0, 0.0, 0.0
    
    # 1. Hit Ratio
    hit_ratio = total_hits / total_reads
    
    # 2. Position Purity (Mode % of hits)
    pos_counts = Counter(positions)
    mode_pos, mode_count = pos_counts.most_common(1)[0]
    pos_purity = mode_count / total_hits
    
    # 3. Flank Purity (Top sequence % at mode position)
    flanks = upstream_map.get(mode_pos, [])
    if flanks:
        top_flank_count = Counter(flanks).most_common(1)[0][1]
        flank_purity = top_flank_count / len(flanks)
    else:
        flank_purity = 0.0
        
    # 4. Diversity (1 - Gini)
    counts_arr = np.array(list(guide_hits.values()))
    gini = calculate_gini(counts_arr)
    diversity_score = 1.0 - gini
    
    # WEIGHTED SUM
    score = (3.0 * hit_ratio) + (2.0 * pos_purity) + (1.0 * flank_purity) + (1.0 * diversity_score)
    
    return score, hit_ratio, pos_purity, flank_purity, gini

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read1', nargs='+', required=True)
    parser.add_argument('--read2', nargs='+', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--max_reads', type=int, default=10000)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()

    if not os.path.exists('seqSpec_plots'): os.makedirs('seqSpec_plots')
    
    # Load Metadata
    sep = ',' if args.metadata.endswith('.csv') else '\t'
    meta = pd.read_csv(args.metadata, sep=sep)
    col = next((c for c in ['spacer','sequence','guide','seq'] if c in meta.columns), None)
    if not col: raise ValueError("No guide column found in metadata")
    
    guides_fwd = [str(x).upper() for x in meta[col].dropna().unique()]
    guides_rev = [get_reverse_complement(g) for g in guides_fwd]
    
    stats_out = []
    best_configs = []
    
    nrows = len(args.read1) * 2
    fig, axs = plt.subplots(nrows=nrows, ncols=7, figsize=(26, 3.5 * nrows))
    if nrows == 1: axs = axs.reshape(1, -1)
    
    for i, (r1, r2) in enumerate(zip(args.read1, args.read2)):
        bn1 = os.path.basename(r1).split('.')[0]
        bn2 = os.path.basename(r2).split('.')[0]
        
        print(f"Processing {bn1} / {bn2}...")
        s1 = readFastq(r1, args.max_reads)
        s2 = readFastq(r2, args.max_reads)
        
        # We need to collect data for 4 configs, score them, THEN plot.
        
        # Config 1: R1 Fwd
        p1f, u1f, h1f = analyze_guides_in_reads(s1, guides_fwd)
        sc1f, _, _, _, _ = calculate_score(len(s1), p1f, u1f, h1f)
        
        # Config 2: R1 Rev
        p1r, u1r, h1r = analyze_guides_in_reads(s1, guides_rev)
        sc1r, _, _, _, _ = calculate_score(len(s1), p1r, u1r, h1r)
        
        # Config 3: R2 Fwd
        p2f, u2f, h2f = analyze_guides_in_reads(s2, guides_fwd)
        sc2f, _, _, _, _ = calculate_score(len(s2), p2f, u2f, h2f)
        
        # Config 4: R2 Rev
        p2r, u2r, h2r = analyze_guides_in_reads(s2, guides_rev)
        sc2r, _, _, _, _ = calculate_score(len(s2), p2r, u2r, h2r)

        # Determine Winner for this sample pair
        scores = [sc1f, sc1r, sc2f, sc2r]
        labels = ['R1_Fwd', 'R1_Rev', 'R2_Fwd', 'R2_Rev']
        best_idx = np.argmax(scores)
        winner_label = labels[best_idx]
        
        print(f"  Scores: R1F={sc1f:.2f}, R1R={sc1r:.2f}, R2F={sc2f:.2f}, R2R={sc2r:.2f} -> Winner: {winner_label}")
        
        best_configs.append({
            'Sample_R1': bn1, 
            'Sample_R2': bn2, 
            'Best_Config': winner_label, 
            'Score': max(scores),
            'Equation': '3*HitRatio + 2*PosPurity + 1*FlankPurity + 1*(1-Gini)'
        })

        # --- PLOTTING ---
        # Helper to check if this specific subplot is part of the winning combo
        def is_winner(curr_label): return curr_label == winner_label

        # Row 1 (Read 1)
        r_idx = i * 2
        fastq_sequence_plot(s1, bn1, axs[r_idx, 0])
        plot_positions(p1f, "R1 Pos (Fwd)", axs[r_idx, 1], highlight=is_winner('R1_Fwd'))
        plot_top_upstream(u1f, p1f, "Upstream (Fwd)", axs[r_idx, 2], highlight=is_winner('R1_Fwd'))
        plot_guide_diversity(h1f, "Div (Fwd)", axs[r_idx, 3], highlight=is_winner('R1_Fwd'))
        
        plot_positions(p1r, "R1 Pos (Rev)", axs[r_idx, 4], highlight=is_winner('R1_Rev'))
        plot_top_upstream(u1r, p1r, "Upstream (Rev)", axs[r_idx, 5], highlight=is_winner('R1_Rev'))
        plot_guide_diversity(h1r, "Div (Rev)", axs[r_idx, 6], highlight=is_winner('R1_Rev'))

        # Row 2 (Read 2)
        r_idx = i * 2 + 1
        fastq_sequence_plot(s2, bn2, axs[r_idx, 0])
        plot_positions(p2f, "R2 Pos (Fwd)", axs[r_idx, 1], highlight=is_winner('R2_Fwd'))
        plot_top_upstream(u2f, p2f, "Upstream (Fwd)", axs[r_idx, 2], highlight=is_winner('R2_Fwd'))
        plot_guide_diversity(h2f, "Div (Fwd)", axs[r_idx, 3], highlight=is_winner('R2_Fwd'))
        
        plot_positions(p2r, "R2 Pos (Rev)", axs[r_idx, 4], highlight=is_winner('R2_Rev'))
        plot_top_upstream(u2r, p2r, "Upstream (Rev)", axs[r_idx, 5], highlight=is_winner('R2_Rev'))
        plot_guide_diversity(h2r, "Div (Rev)", axs[r_idx, 6], highlight=is_winner('R2_Rev'))
        
        # Save detailed stats
        def add_stat(label, p, u, h, s):
            sc, hr, pp, fp, g = calculate_score(len(s1) if 'R1' in label else len(s2), p, u, h)
            stats_out.append({
                'Sample': bn1, 'Config': label, 'TotalHits': len(p), 
                'HitRatio': hr, 'PosPurity': pp, 'FlankPurity': fp, 
                'Gini': g, 'FinalScore': sc, 'IsWinner': (label == winner_label)
            })
            
        add_stat('R1_Fwd', p1f, u1f, h1f, s1)
        add_stat('R1_Rev', p1r, u1r, h1r, s1)
        add_stat('R2_Fwd', p2f, u2f, h2f, s2)
        add_stat('R2_Rev', p2r, u2r, h2r, s2)

    plt.tight_layout()
    if args.plot:
        plt.savefig("seqSpec_plots/seqSpec_check_plots.png", dpi=300)

    if stats_out:
        pd.DataFrame(stats_out).to_csv("position_table.csv", index=False)
    
    if best_configs:
        pd.DataFrame(best_configs).to_csv("best_config.csv", index=False)
        print("Saved best_config.csv")

if __name__ == "__main__":
    main()