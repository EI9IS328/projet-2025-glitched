#!/usr/bin/env python3
"""
Benchmark script to compare ad-hoc vs in-situ export approaches for HPC simulation.

This script runs the SEM simulation with different configurations and grid sizes,
collecting timing metrics and data output statistics for analysis.
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
    mode: str
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
    mode: str,
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
        mode: Export mode - 'insitu' or 'adhoc'
        ex, ey, ez: Grid dimensions
        timemax: Maximum simulation time
        output_dir: Directory for temporary output files

    Returns:
        BenchmarkResult if successful, None otherwise
    """
    snapshot_dir = os.path.join(output_dir, f"snapshots_{mode}_{ex}x{ey}x{ez}")

    # Build command based on mode
    cmd = [
        bin_path,
        "--ex", str(ex),
        "--ey", str(ey),
        "--ez", str(ez),
        "--timemax", str(timemax),
    ]

    if mode == 'base':
        pass  # No snapshot arguments
    elif mode == 'adhoc-plain':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-format", "plain"])
    elif mode == 'adhoc-bin':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-format", "bin"])
    elif mode == 'adhoc-bin-rle':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-format", "bin", "--snapshot-rle"])
    elif mode == 'adhoc-bin-quant1':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-format", "bin", "--compression", "quant", "--quant-level", "1"])
    elif mode == 'adhoc-bin-quant2':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-format", "bin", "--compression", "quant", "--quant-level", "2"])
    elif mode == 'adhoc-bin-quant1-rle':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-format", "bin", "--compression", "quant_rle", "--quant-level", "1"])
    elif mode == 'adhoc-bin-quant2-rle':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-format", "bin", "--compression", "quant_rle", "--quant-level", "2"])
    elif mode == 'insitu':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-in-situ"])
    elif mode == 'rgb':
        os.makedirs(snapshot_dir, exist_ok=True)
        cmd.extend(["--snapshot-folder-path", snapshot_dir, "--snapshot-in-situ", "--snapshot-colormap", "viridis"])

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
            mode=mode,
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
    modes: list[str],
    grid_sizes: list[tuple[int, int, int]],
    timemax: float,
    output_dir: str,
    num_runs: int = 1,
) -> list[BenchmarkResult]:
    """
    Run benchmarks for all combinations of modes and grid sizes.

    Args:
        bin_path: Path to the semproxy executable
        modes: List of export modes to test
        grid_sizes: List of (ex, ey, ez) tuples
        timemax: Maximum simulation time
        output_dir: Directory for temporary output files
        num_runs: Number of runs per configuration (for averaging)

    Returns:
        List of BenchmarkResult objects
    """
    results = []
    total_configs = len(modes) * len(grid_sizes) * num_runs
    current = 0

    for mode, (ex, ey, ez) in itertools.product(modes, grid_sizes):
        for run_idx in range(num_runs):
            current += 1
            print(f"[{current}/{total_configs}] Mode: {mode}, Grid: {ex}x{ey}x{ez}, Run: {run_idx + 1}/{num_runs}")

            result = run_benchmark(
                bin_path=bin_path,
                mode=mode,
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
                      f"Bytes: {result.total_bytes:.0f}")
            else:
                print(f"    FAILED")

    return results


def save_results_csv(results: list[BenchmarkResult], output_path: str):
    """Save benchmark results to a CSV file."""
    fieldnames = [
        'mode', 'ex', 'ey', 'ez', 'grid_total',
        'time_simulating', 'time_snapshots', 'time_sismos',
        'time_kernel_total', 'total_bytes',
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for r in results:
            writer.writerow({
                'mode': r.mode,
                'ex': r.ex,
                'ey': r.ey,
                'ez': r.ez,
                'grid_total': r.ex * r.ey * r.ez,
                'time_simulating': r.time_simulating,
                'time_snapshots': r.time_snapshots,
                'time_sismos': r.time_sismos,
                'time_kernel_total': r.time_kernel_total,
                'total_bytes': r.total_bytes,
            })

    print(f"Results saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark in-situ vs ad-hoc export approaches for HPC simulation"
    )
    parser.add_argument(
        "--bin", "-b",
        default="./build/bin/semproxy",
        help="Path to semproxy executable (default: ./build/bin/semproxy)"
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output CSV file path (default: benchmark_TIMESTAMP.csv)"
    )
    parser.add_argument(
        "--output-dir",
        default="./benchmark_tmp",
        help="Directory for temporary output files (default: ./benchmark_tmp)"
    )
    parser.add_argument(
        "--modes", "-m",
        nargs="+",
        choices=["base", "adhoc-plain", "adhoc-bin", "adhoc-bin-rle", "adhoc-bin-quant1", "adhoc-bin-quant2", "adhoc-bin-quant1-rle", "adhoc-bin-quant2-rle", "insitu", "rgb"],
        default=["base", "adhoc-plain", "adhoc-bin", "adhoc-bin-rle", "adhoc-bin-quant1", "adhoc-bin-quant2", "adhoc-bin-quant1-rle", "adhoc-bin-quant2-rle", "insitu", "rgb"],
        help="Export modes to benchmark: base (no snapshots), adhoc-plain (plain text), adhoc-bin (binary), adhoc-bin-rle (binary with RLE), adhoc-bin-quant1/2 (binary with quantization), adhoc-bin-quant1/2-rle (quantization + RLE), insitu (image), rgb (colormap) (default: all)"
    )
    parser.add_argument(
        "--sizes", "-s",
        nargs="+",
        type=int,
        default=[10, 20, 30, 40, 50],
        help="Grid sizes to test (same for ex, ey, ez) (default: 10 20 30 40 50)"
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
        default=1,
        help="Number of runs per configuration (default: 1)"
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
        args.output = f"benchmark_{timestamp}.csv"

    print(f"\nBenchmark Configuration:")
    print(f"  Binary: {args.bin}")
    print(f"  Modes: {args.modes}")
    print(f"  Grid sizes: {grid_sizes}")
    print(f"  Timemax: {args.timemax}s")
    print(f"  Runs per config: {args.runs}")
    print(f"  Output: {args.output}")
    print()

    # Run benchmarks
    results = run_benchmarks(
        bin_path=args.bin,
        modes=args.modes,
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

    # Cleanup
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)

    print(f"\nBenchmark complete! {len(results)} results saved to {args.output}")


if __name__ == "__main__":
    main()
