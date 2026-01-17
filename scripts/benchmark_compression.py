#!/usr/bin/env python3
"""
Benchmark script to compare RLE compression vs no compression for binary snapshots.

This script focuses on comparing:
- adhoc-bin: Binary snapshots without compression
- adhoc-bin-rle: Binary snapshots with RLE compression

Metrics collected:
- Execution time (simulation, snapshots, total)
- Total data exported (bytes)
"""

import argparse
import csv
import itertools
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from typing import Optional

# Regex patterns to extract metrics from output
NUMBER_REGEX = r"(\d+(?:\.\d+)?(?:e[+-]?\d+)?)"
TIME_SIMULATING_REGEX = rf"Time Spent Simulating : {NUMBER_REGEX} seconds\."
TIME_SNAPSHOTS_REGEX = rf"Time Making and Saving Snapshots : {NUMBER_REGEX} seconds\."
TIME_SISMOS_REGEX = rf"Time Making Sismos : {NUMBER_REGEX} seconds\."
TOTAL_BYTES_REGEX = rf"Total written data: {NUMBER_REGEX} Bytes\."
TIME_KERNEL_TOTAL_REGEX = rf"Time Kernel Total : {NUMBER_REGEX} seconds\."


@dataclass
class BenchmarkResult:
    """Holds the results of a single benchmark run."""
    compression: str  # 'none' or 'rle'
    ex: int
    ey: int
    ez: int
    time_simulating: float
    time_snapshots: float
    time_sismos: float
    time_kernel_total: float
    total_bytes: float


def parse_output(output: str) -> dict:
    """Extract timing metrics from simulation output."""
    metrics = {}

    patterns = {
        'time_simulating': TIME_SIMULATING_REGEX,
        'time_snapshots': TIME_SNAPSHOTS_REGEX,
        'time_sismos': TIME_SISMOS_REGEX,
        'time_kernel_total': TIME_KERNEL_TOTAL_REGEX,
        'total_bytes': TOTAL_BYTES_REGEX,
    }

    for name, pattern in patterns.items():
        match = re.search(pattern, output)
        if match:
            metrics[name] = float(match.group(1))
        else:
            metrics[name] = 0.0

    return metrics


def run_benchmark(
    bin_path: str,
    compression: str,
    ex: int,
    ey: int,
    ez: int,
    timemax: float,
    output_dir: str,
) -> Optional[BenchmarkResult]:
    """
    Run a single benchmark with the given configuration.

    Args:
        bin_path: Path to the semproxy executable
        compression: Compression mode ('none', 'rle', 'quant1', 'quant2', 'quant1_rle', 'quant2_rle')
        ex, ey, ez: Grid dimensions
        timemax: Maximum simulation time
        output_dir: Directory for temporary output files

    Returns:
        BenchmarkResult if successful, None otherwise
    """
    snapshot_dir = os.path.join(output_dir, f"snapshots_{compression}_{ex}x{ey}x{ez}")
    os.makedirs(snapshot_dir, exist_ok=True)

    # Build command
    cmd = [
        bin_path,
        "--ex", str(ex),
        "--ey", str(ey),
        "--ez", str(ez),
        "--timemax", str(timemax),
        "--snapshot-folder-path", snapshot_dir,
        "--snapshot-format", "bin",
    ]

    if compression == 'rle':
        cmd.append("--snapshot-rle")
    elif compression == 'quant1':
        cmd.extend(["--compression", "quant", "--quant-level", "1"])
    elif compression == 'quant2':
        cmd.extend(["--compression", "quant", "--quant-level", "2"])
    elif compression == 'quant1_rle':
        cmd.extend(["--compression", "quant_rle", "--quant-level", "1"])
    elif compression == 'quant2_rle':
        cmd.extend(["--compression", "quant_rle", "--quant-level", "2"])

    print(f"  Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 minute timeout
        )

        if result.returncode != 0:
            print(f"  ERROR: Process returned {result.returncode}", file=sys.stderr)
            print(f"  stderr: {result.stderr}", file=sys.stderr)
            return None

        output = result.stdout + result.stderr
        metrics = parse_output(output)

        return BenchmarkResult(
            compression=compression,
            ex=ex,
            ey=ey,
            ez=ez,
            time_simulating=metrics['time_simulating'],
            time_snapshots=metrics['time_snapshots'],
            time_sismos=metrics['time_sismos'],
            time_kernel_total=metrics['time_kernel_total'],
            total_bytes=metrics['total_bytes'],
        )

    except subprocess.TimeoutExpired:
        print(f"  ERROR: Process timed out", file=sys.stderr)
        return None
    except Exception as e:
        print(f"  ERROR: {e}", file=sys.stderr)
        return None
    finally:
        # Cleanup snapshot directory
        if os.path.exists(snapshot_dir):
            shutil.rmtree(snapshot_dir)


def run_benchmarks(
    bin_path: str,
    compressions: list[str],
    grid_sizes: list[tuple[int, int, int]],
    timemax: float,
    output_dir: str,
    num_runs: int = 1,
) -> list[BenchmarkResult]:
    """
    Run benchmarks for all combinations of compression modes and grid sizes.

    Args:
        bin_path: Path to the semproxy executable
        compressions: List of compression modes to test
        grid_sizes: List of (ex, ey, ez) tuples
        timemax: Maximum simulation time
        output_dir: Directory for temporary output files
        num_runs: Number of runs per configuration (for averaging)

    Returns:
        List of BenchmarkResult objects
    """
    results = []
    total_configs = len(compressions) * len(grid_sizes) * num_runs
    current = 0

    for compression, (ex, ey, ez) in itertools.product(compressions, grid_sizes):
        for run_idx in range(num_runs):
            current += 1
            print(f"[{current}/{total_configs}] Compression: {compression}, Grid: {ex}x{ey}x{ez}, Run: {run_idx + 1}/{num_runs}")

            result = run_benchmark(
                bin_path=bin_path,
                compression=compression,
                ex=ex,
                ey=ey,
                ez=ez,
                timemax=timemax,
                output_dir=output_dir,
            )

            if result:
                results.append(result)
                print(f"    Simulating: {result.time_simulating:.4f}s, "
                      f"Snapshots: {result.time_snapshots:.4f}s, "
                      f"Total: {result.time_kernel_total:.4f}s, "
                      f"Data: {result.total_bytes / 1024 / 1024:.2f} MB")
            else:
                print(f"    FAILED")

    return results


def save_results_csv(results: list[BenchmarkResult], output_path: str):
    """Save benchmark results to a CSV file."""
    fieldnames = [
        'compression', 'ex', 'ey', 'ez', 'grid_total',
        'time_simulating', 'time_snapshots', 'time_sismos',
        'time_kernel_total', 'total_bytes', 'total_mb',
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for r in results:
            writer.writerow({
                'compression': r.compression,
                'ex': r.ex,
                'ey': r.ey,
                'ez': r.ez,
                'grid_total': r.ex * r.ey * r.ez,
                'time_simulating': r.time_simulating,
                'time_snapshots': r.time_snapshots,
                'time_sismos': r.time_sismos,
                'time_kernel_total': r.time_kernel_total,
                'total_bytes': r.total_bytes,
                'total_mb': r.total_bytes / 1024 / 1024,
            })

    print(f"\nResults saved to {output_path}")


def print_summary(results: list[BenchmarkResult]):
    """Print a summary comparison of all compression modes."""
    if not results:
        return

    print("\n" + "="*80)
    print("SUMMARY: Compression Comparison")
    print("="*80)

    # Group by grid size
    grid_sizes = sorted(set((r.ex, r.ey, r.ez) for r in results))
    compression_modes = sorted(set(r.compression for r in results))

    for ex, ey, ez in grid_sizes:
        grid_results = [r for r in results if (r.ex, r.ey, r.ez) == (ex, ey, ez)]

        # Get baseline (no compression)
        none_results = [r for r in grid_results if r.compression == 'none']
        if not none_results:
            continue

        none_time = sum(r.time_kernel_total for r in none_results) / len(none_results)
        none_bytes = sum(r.total_bytes for r in none_results) / len(none_results)

        print(f"\nGrid {ex}x{ey}x{ez} ({ex*ey*ez} elements):")
        print(f"  {'Mode':<20} {'Time (s)':<12} {'Ratio':<10} {'Data (MB)':<12} {'Compression':<12} {'Reduction':<10}")
        print(f"  {'-'*20} {'-'*12} {'-'*10} {'-'*12} {'-'*12} {'-'*10}")

        for mode in compression_modes:
            mode_results = [r for r in grid_results if r.compression == mode]
            if not mode_results:
                continue

            mode_time = sum(r.time_kernel_total for r in mode_results) / len(mode_results)
            mode_bytes = sum(r.total_bytes for r in mode_results) / len(mode_results)

            time_ratio = mode_time / none_time if none_time > 0 else 0
            compression_ratio = mode_bytes / none_bytes if none_bytes > 0 else 0
            reduction = (1 - compression_ratio) * 100

            print(f"  {mode:<20} {mode_time:<12.4f} {time_ratio:<10.2f} {mode_bytes/1024/1024:<12.2f} {compression_ratio:<12.2f} {reduction:<10.1f}%")

    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark compression methods for binary snapshots"
    )
    parser.add_argument(
        "--bin", "-b",
        default="./build/bin/semproxy",
        help="Path to semproxy executable (default: ./build/bin/semproxy)"
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output CSV file path (default: benchmark_compression_TIMESTAMP.csv)"
    )
    parser.add_argument(
        "--output-dir",
        default="./benchmark_compression_tmp",
        help="Directory for temporary output files (default: ./benchmark_compression_tmp)"
    )
    parser.add_argument(
        "--compressions", "-c",
        nargs="+",
        choices=["none", "rle", "quant1", "quant2", "quant1_rle", "quant2_rle"],
        default=["none", "rle", "quant1", "quant2", "quant1_rle", "quant2_rle"],
        help="Compression modes to test (default: all)"
    )
    parser.add_argument(
        "--sizes", "-s",
        nargs="+",
        type=int,
        default=[10, 20, 30, 40, 50, 75, 100, 125, 150],
        help="Grid sizes to test (same for ex, ey, ez) (default: 10 20 30 40 50 75 100 125 150)"
    )
    parser.add_argument(
        "--timemax", "-t",
        type=float,
        default=0.5,
        help="Maximum simulation time in seconds (default: 0.5)"
    )
    parser.add_argument(
        "--runs", "-r",
        type=int,
        default=3,
        help="Number of runs per configuration for averaging (default: 3)"
    )

    args = parser.parse_args()

    # Check binary exists
    if not os.path.exists(args.bin):
        print(f"Error: Binary not found at {args.bin}", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate grid sizes (cubic grids)
    grid_sizes = [(s, s, s) for s in args.sizes]

    # Set output filename
    if args.output is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.output = f"benchmark_compression_{timestamp}.csv"

    print(f"\nCompression Benchmark Configuration:")
    print(f"  Binary: {args.bin}")
    print(f"  Compression modes: {args.compressions}")
    print(f"  Grid sizes: {grid_sizes}")
    print(f"  Timemax: {args.timemax}s")
    print(f"  Runs per config: {args.runs}")
    print(f"  Output: {args.output}")
    print()

    # Run benchmarks
    results = run_benchmarks(
        bin_path=args.bin,
        compressions=args.compressions,
        grid_sizes=grid_sizes,
        timemax=args.timemax,
        output_dir=args.output_dir,
        num_runs=args.runs,
    )

    if not results:
        print("No results collected!", file=sys.stderr)
        sys.exit(1)

    # Save results
    save_results_csv(results, args.output)

    # Print summary
    print_summary(results)

    # Cleanup
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)

    print(f"\nBenchmark complete! {len(results)} results saved to {args.output}")


if __name__ == "__main__":
    main()
