#!/usr/bin/env python3
"""
Visualize benchmark results from benchmark_insitu_vs_adhoc.py
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd


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

    # Plot 1: Snapshot time vs grid size
    ax = axes[0, 0]
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        ax.plot(data['grid_total'], data['time_snapshots'], 'o-', label=mode)
    ax.set_xlabel('Grid Total (elements)')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Snapshot Time vs Grid Size')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Total bytes vs grid size (log scale)
    ax = axes[0, 1]
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        ax.plot(data['grid_total'], data['total_bytes'] / 1e6, 'o-', label=mode)
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
    breakdown.plot(kind='bar', stacked=True, ax=ax)
    ax.set_xlabel('Mode')
    ax.set_ylabel('Time (seconds)')
    ax.set_title(f'Time Breakdown (grid={largest_grid})')
    ax.legend(['Simulating', 'Snapshots', 'Sismos'])
    ax.tick_params(axis='x', rotation=0)

    # Plot 4: Kernel time vs grid size
    ax = axes[1, 1]
    for mode in modes:
        data = grouped[grouped['mode'] == mode]
        ax.plot(data['grid_total'], data['time_kernel_total'], 'o-', label=mode)
    ax.set_xlabel('Grid Total (elements)')
    ax.set_ylabel('Time (seconds)')
    ax.set_title('Kernel Time vs Grid Size')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"Saved to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
