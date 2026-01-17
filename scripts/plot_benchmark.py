#!/usr/bin/env python3
"""
Visualize benchmark results from benchmark_insitu_vs_adhoc.py
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd

# Scientific color palette for HPC/simulation context
MODE_STYLES = {
    'base': {'color': '#2c3e50', 'marker': 'o', 'linestyle': '-'},
    'adhoc-plain': {'color': '#e74c3c', 'marker': 's', 'linestyle': '-'},
    'adhoc-bin': {'color': '#c0392b', 'marker': 'p', 'linestyle': '-.'},
    'adhoc-bin-rle': {'color': '#8e44ad', 'marker': 'h', 'linestyle': '-.'},
    'adhoc-bin-quant1': {'color': '#e67e22', 'marker': 'v', 'linestyle': '--'},
    'adhoc-bin-quant2': {'color': '#d35400', 'marker': '<', 'linestyle': '--'},
    'adhoc-bin-quant1-rle': {'color': '#9b59b6', 'marker': '>', 'linestyle': ':'},
    'adhoc-bin-quant2-rle': {'color': '#6c3483', 'marker': '*', 'linestyle': ':'},
    'insitu': {'color': '#3498db', 'marker': '^', 'linestyle': '--'},
    'rgb': {'color': '#27ae60', 'marker': 'D', 'linestyle': ':'},
}

MODE_LABELS = {
    'base': 'Base (no export)',
    'adhoc-plain': 'Ad-hoc (plain text)',
    'adhoc-bin': 'Ad-hoc (binary)',
    'adhoc-bin-rle': 'Ad-hoc (binary RLE)',
    'adhoc-bin-quant1': 'Ad-hoc (binary quant L1)',
    'adhoc-bin-quant2': 'Ad-hoc (binary quant L2)',
    'adhoc-bin-quant1-rle': 'Ad-hoc (quant L1 + RLE)',
    'adhoc-bin-quant2-rle': 'Ad-hoc (quant L2 + RLE)',
    'insitu': 'In-situ (grayscale)',
    'rgb': 'In-situ (RGB colormap)',
}


def get_style(mode):
    return MODE_STYLES.get(mode, {'color': '#7f8c8d', 'marker': 'o', 'linestyle': '-'})


def get_label(mode):
    return MODE_LABELS.get(mode, mode)


def main():
    parser = argparse.ArgumentParser(
        description="Plot benchmark results from CSV"
    )
    parser.add_argument(
        "csv_file",
        help="Path to benchmark CSV file"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)

    # Group by mode and grid size, averaging multiple runs
    grouped = df.groupby(['mode', 'grid_total']).mean(numeric_only=True).reset_index()

    modes = grouped['mode'].unique()

    # Plot 1: Snapshot time vs grid size
    fig, ax = plt.subplots(figsize=(6, 5))
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        style = get_style(mode)
        ax.plot(data['grid_total'], data['time_snapshots'],
                marker=style['marker'], linestyle=style['linestyle'],
                color=style['color'], label=get_label(mode), linewidth=2, markersize=8)
    ax.set_xlabel('Grid Total (elements)')
    ax.set_ylabel('Time (seconds)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('snapshot_time_vs_grid_size.png', dpi=150)
    print("Saved snapshot_time_vs_grid_size.png")
    plt.close()

    # Plot 2: Total bytes vs grid size (log scale)
    fig, ax = plt.subplots(figsize=(6, 5))
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        style = get_style(mode)
        ax.plot(data['grid_total'], data['total_bytes'] / 1e6,
                marker=style['marker'], linestyle=style['linestyle'],
                color=style['color'], label=get_label(mode), linewidth=2, markersize=8)
    ax.set_xlabel('Grid Total (elements)')
    ax.set_ylabel('Data Written (MB)')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('data_output_vs_grid_size.png', dpi=150)
    print("Saved data_output_vs_grid_size.png")
    plt.close()

    # Plot 3: Time breakdown (stacked bar) for largest grid
    fig, ax = plt.subplots(figsize=(6, 5))
    largest_grid = grouped['grid_total'].max()
    breakdown = grouped[grouped['grid_total'] == largest_grid][['mode', 'time_simulating', 'time_snapshots', 'time_sismos']]
    breakdown = breakdown.set_index('mode')
    # Reorder index to match our preferred order
    order = [m for m in ['base', 'adhoc-plain', 'adhoc-bin', 'adhoc-bin-rle', 'adhoc-bin-quant1', 'adhoc-bin-quant2', 'adhoc-bin-quant1-rle', 'adhoc-bin-quant2-rle', 'insitu', 'rgb'] if m in breakdown.index]
    breakdown = breakdown.reindex(order)
    breakdown.index = [get_label(m) for m in breakdown.index]
    breakdown.plot(kind='bar', stacked=True, ax=ax, color=['#3498db', '#e74c3c', '#f39c12'])
    ax.set_xlabel('Mode')
    ax.set_ylabel('Time (seconds)')
    ax.legend(['Simulating', 'Snapshots', 'Sismos'])
    ax.tick_params(axis='x', rotation=15)
    plt.tight_layout()
    plt.savefig(f'time_breakdown_grid_{largest_grid}.png', dpi=150)
    print(f"Saved time_breakdown_grid_{largest_grid}.png")
    plt.close()

    # Plot 4: Total time vs grid size
    fig, ax = plt.subplots(figsize=(6, 5))
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        style = get_style(mode)
        ax.plot(data['grid_total'], data['time_kernel_total'],
                marker=style['marker'], linestyle=style['linestyle'],
                color=style['color'], label=get_label(mode), linewidth=2, markersize=8)
    ax.set_xlabel('Grid Total (elements)')
    ax.set_ylabel('Time (seconds)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('total_time_vs_grid_size.png', dpi=150)
    print("Saved total_time_vs_grid_size.png")
    plt.close()


if __name__ == "__main__":
    main()
