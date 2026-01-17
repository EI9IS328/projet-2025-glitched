#!/usr/bin/env python3
"""
Visualize compression benchmark results comparing RLE vs no compression.
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd

# Color palette for compression comparison
COMPRESSION_STYLES = {
    'none': {'color': '#c0392b', 'marker': 's', 'linestyle': '-', 'label': 'No Compression'},
    'rle': {'color': '#8e44ad', 'marker': 'h', 'linestyle': '-.', 'label': 'RLE'},
    'quant1': {'color': '#e67e22', 'marker': 'v', 'linestyle': '--', 'label': 'Quant L1'},
    'quant2': {'color': '#d35400', 'marker': '<', 'linestyle': '--', 'label': 'Quant L2'},
    'quant1_rle': {'color': '#9b59b6', 'marker': '>', 'linestyle': ':', 'label': 'Quant L1 + RLE'},
    'quant2_rle': {'color': '#6c3483', 'marker': '*', 'linestyle': ':', 'label': 'Quant L2 + RLE'},
}


def main():
    parser = argparse.ArgumentParser(
        description="Plot compression benchmark results from CSV"
    )
    parser.add_argument(
        "csv_file",
        help="Path to benchmark CSV file"
    )
    parser.add_argument(
        "--output-prefix", "-o",
        default="compression",
        help="Prefix for output files (default: compression)"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)

    # Group by compression type and grid size, averaging multiple runs
    grouped = df.groupby(['compression', 'grid_total']).mean(numeric_only=True).reset_index()

    compressions = grouped['compression'].unique()

    # Plot 1: Execution time comparison
    fig, ax = plt.subplots(figsize=(8, 6))
    for comp in compressions:
        data = grouped[grouped['compression'] == comp]
        style = COMPRESSION_STYLES.get(comp, COMPRESSION_STYLES['none'])
        ax.plot(data['grid_total'], data['time_kernel_total'],
                marker=style['marker'], linestyle=style['linestyle'],
                color=style['color'], label=style['label'], linewidth=2, markersize=8)
    ax.set_xlabel('Grid Total (elements)', fontsize=12)
    ax.set_ylabel('Total Time (seconds)', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = f'{args.output_prefix}_time_comparison.png'
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close()

    # Plot 2: Data output comparison (log scale)
    fig, ax = plt.subplots(figsize=(8, 6))
    for comp in compressions:
        data = grouped[grouped['compression'] == comp]
        style = COMPRESSION_STYLES.get(comp, COMPRESSION_STYLES['none'])
        ax.plot(data['grid_total'], data['total_mb'],
                marker=style['marker'], linestyle=style['linestyle'],
                color=style['color'], label=style['label'], linewidth=2, markersize=8)
    ax.set_xlabel('Grid Total (elements)', fontsize=12)
    ax.set_ylabel('Data Written (MB)', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = f'{args.output_prefix}_data_comparison.png'
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close()

    # Plot 3: Compression ratio (data reduction) - all modes vs baseline
    fig, ax = plt.subplots(figsize=(8, 6))
    none_data = grouped[grouped['compression'] == 'none'].sort_values('grid_total')

    if len(none_data) > 0:
        for comp in compressions:
            if comp == 'none':
                continue
            comp_data = grouped[grouped['compression'] == comp].sort_values('grid_total')

            if len(comp_data) > 0:
                # Ensure both datasets have the same grid sizes
                common_grids = set(none_data['grid_total']).intersection(set(comp_data['grid_total']))
                none_subset = none_data[none_data['grid_total'].isin(common_grids)].sort_values('grid_total')
                comp_subset = comp_data[comp_data['grid_total'].isin(common_grids)].sort_values('grid_total')

                # Compression ratio = uncompressed / compressed (higher is better)
                compression_ratio = none_subset['total_mb'].values / comp_subset['total_mb'].values

                style = COMPRESSION_STYLES.get(comp, COMPRESSION_STYLES['none'])
                ax.plot(none_subset['grid_total'], compression_ratio,
                        marker=style['marker'], linestyle=style['linestyle'],
                        color=style['color'], label=style['label'], linewidth=2, markersize=8)

        ax.axhline(y=1.0, color='#95a5a6', linestyle='--', linewidth=1, label='No compression (1.0x)')
        ax.set_xlabel('Grid Total (elements)', fontsize=12)
        ax.set_ylabel('Compression Ratio', fontsize=12)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        filename = f'{args.output_prefix}_ratio.png'
        plt.savefig(filename, dpi=150)
        print(f"Saved {filename}")
        plt.close()

    # Plot 4: Snapshot time breakdown
    fig, ax = plt.subplots(figsize=(8, 6))
    for comp in compressions:
        data = grouped[grouped['compression'] == comp]
        style = COMPRESSION_STYLES.get(comp, COMPRESSION_STYLES['none'])
        ax.plot(data['grid_total'], data['time_snapshots'],
                marker=style['marker'], linestyle=style['linestyle'],
                color=style['color'], label=style['label'], linewidth=2, markersize=8)
    ax.set_xlabel('Grid Total (elements)', fontsize=12)
    ax.set_ylabel('Snapshot Time (seconds)', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = f'{args.output_prefix}_snapshot_time.png'
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")
    plt.close()

    # Plot 5: Side-by-side comparison for largest grid
    largest_grid = grouped['grid_total'].max()
    largest_data = grouped[grouped['grid_total'] == largest_grid]

    if len(largest_data) >= 2:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Time comparison
        compressions_list = largest_data['compression'].tolist()
        times = largest_data['time_kernel_total'].tolist()
        colors = [COMPRESSION_STYLES.get(c, COMPRESSION_STYLES['none'])['color'] for c in compressions_list]
        labels = [COMPRESSION_STYLES.get(c, COMPRESSION_STYLES['none'])['label'] for c in compressions_list]

        ax1.bar(labels, times, color=colors, alpha=0.7, edgecolor='black')
        ax1.set_ylabel('Total Time (seconds)', fontsize=12)
        ax1.tick_params(axis='x', rotation=45)
        ax1.grid(True, alpha=0.3, axis='y')

        # Data size comparison
        data_sizes = largest_data['total_mb'].tolist()

        ax2.bar(labels, data_sizes, color=colors, alpha=0.7, edgecolor='black')
        ax2.set_ylabel('Data Written (MB)', fontsize=12)
        ax2.tick_params(axis='x', rotation=45)
        ax2.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        filename = f'{args.output_prefix}_comparison_grid_{largest_grid}.png'
        plt.savefig(filename, dpi=150)
        print(f"Saved {filename}")
        plt.close()

    print("\nAll plots generated successfully!")


if __name__ == "__main__":
    main()
