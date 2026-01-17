#!/usr/bin/env python3
"""
Visualize benchmark results from benchmark_insitu_vs_adhoc.py
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd

# Scientific color palette for HPC/simulation context
MODE_COLORS = {
    'base': '#2c3e50',    # Dark slate - baseline
    'adhoc': '#e74c3c',   # Red - raw data export
    'insitu': '#3498db',  # Blue - in-situ processing
    'rgb': '#27ae60',     # Green - colormap export
}

MODE_LABELS = {
    'base': 'Base (no export)',
    'adhoc': 'Ad-hoc (raw data)',
    'insitu': 'In-situ (grayscale)',
    'rgb': 'In-situ (RGB colormap)',
}


def get_color(mode):
    return MODE_COLORS.get(mode, '#7f8c8d')


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
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output image file (default: display interactively)"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)

    # Group by mode and grid size, averaging multiple runs
    grouped = df.groupby(['mode', 'grid_total']).mean(numeric_only=True).reset_index()

    modes = grouped['mode'].unique()
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('FUnTiDES Benchmark: Export Mode Comparison', fontsize=14, fontweight='bold')

    # Plot 1: Snapshot time vs grid size
    ax = axes[0, 0]
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        ax.plot(data['grid_total'], data['time_snapshots'], 'o-',
                color=get_color(mode), label=get_label(mode), linewidth=2, markersize=6)
    ax.set_xlabel('Grid Total (elements)')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Snapshot Time vs Grid Size')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Total bytes vs grid size (log scale)
    ax = axes[0, 1]
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        ax.plot(data['grid_total'], data['total_bytes'] / 1e6, 'o-',
                color=get_color(mode), label=get_label(mode), linewidth=2, markersize=6)
    ax.set_xlabel('Grid Total (elements)')
    ax.set_ylabel('Data Written (MB)')
    ax.set_yscale('log')
    ax.set_title('Data Output vs Grid Size')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Time breakdown (stacked bar) for largest grid
    ax = axes[1, 0]
    largest_grid = grouped['grid_total'].max()
    breakdown = grouped[grouped['grid_total'] == largest_grid][['mode', 'time_simulating', 'time_snapshots', 'time_sismos']]
    breakdown = breakdown.set_index('mode')
    # Reorder index to match our preferred order
    order = [m for m in ['base', 'adhoc', 'insitu', 'rgb'] if m in breakdown.index]
    breakdown = breakdown.reindex(order)
    breakdown.index = [get_label(m) for m in breakdown.index]
    breakdown.plot(kind='bar', stacked=True, ax=ax, color=['#3498db', '#e74c3c', '#f39c12'])
    ax.set_xlabel('Mode')
    ax.set_ylabel('Time (seconds)')
    ax.set_title(f'Time Breakdown (grid={largest_grid})')
    ax.legend(['Simulating', 'Snapshots', 'Sismos'])
    ax.tick_params(axis='x', rotation=15)

    # Plot 4: Kernel time vs grid size
    ax = axes[1, 1]
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        ax.plot(data['grid_total'], data['time_kernel_total'], 'o-',
                color=get_color(mode), label=get_label(mode), linewidth=2, markersize=6)
    ax.set_xlabel('Grid Total (elements)')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Kernel Time vs Grid Size')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"Saved to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
