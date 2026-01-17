#!/usr/bin/env python3
"""
Quick plot script to compare adhoc-plain vs adhoc-bin from benchmark results.
Creates two plots: execution time and data amount saved.
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Plot adhoc-plain vs adhoc-bin comparison"
    )
    parser.add_argument(
        "csv_file",
        help="Path to benchmark CSV file from benchmark_insitu_vs_adhoc"
    )
    args = parser.parse_args()

    # Read data
    df = pd.read_csv(args.csv_file)

    # Filter for adhoc-plain and adhoc-bin only
    filtered = df[df['mode'].isin(['adhoc-plain', 'adhoc-bin'])]

    if len(filtered) == 0:
        print("Error: No adhoc-plain or adhoc-bin data found in CSV")
        return

    # Group by mode and grid size, averaging multiple runs
    grouped = filtered.groupby(['mode', 'grid_total']).mean(numeric_only=True).reset_index()

    # Split data
    plain_data = grouped[grouped['mode'] == 'adhoc-plain'].sort_values('grid_total')
    binary_data = grouped[grouped['mode'] == 'adhoc-bin'].sort_values('grid_total')

    # Plot 1: Execution time comparison
    fig, ax = plt.subplots(figsize=(8, 6))

    if len(plain_data) > 0:
        ax.plot(plain_data['grid_total'], plain_data['time_kernel_total'],
                marker='s', linestyle='-', color='#e74c3c',
                label='Ad-hoc Plain Text', linewidth=2, markersize=8)

    if len(binary_data) > 0:
        ax.plot(binary_data['grid_total'], binary_data['time_kernel_total'],
                marker='p', linestyle='-.', color='#c0392b',
                label='Ad-hoc Binary', linewidth=2, markersize=8)

    ax.set_xlabel('Grid Total (elements)', fontsize=12)
    ax.set_ylabel('Total Execution Time (seconds)', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('plain_vs_binary_time.png', dpi=150)
    print("Saved plain_vs_binary_time.png")
    plt.close()

    # Plot 2: Data amount comparison
    fig, ax = plt.subplots(figsize=(8, 6))

    if len(plain_data) > 0:
        ax.plot(plain_data['grid_total'], plain_data['total_bytes'] / 1e6,
                marker='s', linestyle='-', color='#e74c3c',
                label='Ad-hoc Plain Text', linewidth=2, markersize=8)

    if len(binary_data) > 0:
        ax.plot(binary_data['grid_total'], binary_data['total_bytes'] / 1e6,
                marker='p', linestyle='-.', color='#c0392b',
                label='Ad-hoc Binary', linewidth=2, markersize=8)

    ax.set_xlabel('Grid Total (elements)', fontsize=12)
    ax.set_ylabel('Data Written (MB)', fontsize=12)
    ax.set_yscale('log')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('plain_vs_binary_data.png', dpi=150)
    print("Saved plain_vs_binary_data.png")
    plt.close()

    # Print summary
    print("\nSummary:")
    if len(plain_data) > 0 and len(binary_data) > 0:
        # Get data for largest grid size
        common_grids = set(plain_data['grid_total']).intersection(set(binary_data['grid_total']))
        if common_grids:
            largest = max(common_grids)
            plain_largest = plain_data[plain_data['grid_total'] == largest].iloc[0]
            binary_largest = binary_data[binary_data['grid_total'] == largest].iloc[0]

            time_speedup = plain_largest['time_kernel_total'] / binary_largest['time_kernel_total']
            data_ratio = binary_largest['total_bytes'] / plain_largest['total_bytes']
            data_reduction = (1 - data_ratio) * 100

            print(f"  Grid size: {largest} elements")
            print(f"  Binary is {time_speedup:.2f}x faster than plain text")
            print(f"  Binary uses {data_ratio:.2f}x the data ({data_reduction:.1f}% reduction)")


if __name__ == "__main__":
    main()
