# snaps  -> évolution pression globale sur le domaine
# x sismos -> évolution pression par receiver suivi
# stats (min, max, mean, stdvar…),
# snaps  -> par étape
#        -> en global
# x sismos -> par receiver
# x       -> en global
# fourier transform on sismic signals
# x sismos -> par receiver
# propagation speed (global)

import argparse
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.stats import linregress
from datetime import datetime
import glob
import re

def load_params(path):
    params = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or ':' not in line:
                continue
            keys_part, vals_part = line.split(':', 1)
            keys = [k.strip() for k in keys_part.split(',')]
            vals = []
            for v in vals_part.split(','):
                try:
                    vals.append(float(v.strip()))
                except ValueError:
                    vals.append(v.strip())
            
            if len(keys) == len(vals):
                for k, v in zip(keys, vals):
                    params[k] = v
            elif len(vals) == 1:
                 for k in keys:
                    params[k] = vals[0]
            else:
                print(f"Warning: mismatch in keys/values count in line: {line}")

    return params

def load_sismos(path):
    with open(path, "r") as f:
        header = f.readline().strip().split(";")
        if len(header) < 2:
             return 0, 0, []
        n_receivers = int(header[0])
        n_measures = int(header[1])

        receivers_data = []
        for _ in range(n_receivers):
            coords_line = f.readline().strip()
            if not coords_line: break
            coords = list(map(float, coords_line.split(";")))
            
            measures_line = f.readline().strip()
            if not measures_line: break
            raw_measures = list(map(float, measures_line.split(";")))
            measures = np.abs(np.array(raw_measures)).tolist()
            receivers_data.append({"coords": coords, "measures": measures})
    return n_receivers, n_measures, receivers_data

def load_snapshot(path, fmt="plain"):
    if fmt == "bin":
        return np.fromfile(path, dtype=np.float32)
    else:
        try:
            return np.loadtxt(path, delimiter=',')
        except ValueError:
             with open(path, 'r') as f:
                 data = []
                 for line in f:
                     parts = line.strip().rstrip(',').split(',')
                     data.extend([float(p) for p in parts])
                 return np.array(data)

def compute_sismos_stats(receivers_data):
    stats_per_receiver = []
    all_measures = []

    for i, recv in enumerate(receivers_data):
        m = np.array(recv["measures"])
        if m.size == 0: continue
        
        all_measures.extend(recv["measures"])
        stats_per_receiver.append(
            {
                "receiver_index": i,
                "min": np.min(m),
                "max": np.max(m),
                "mean": np.mean(m),
                "stdvar": np.std(m),
            }
        )

    all_measures = np.array(all_measures)
    if all_measures.size > 0:
        global_stats = {
            "min": np.min(all_measures),
            "max": np.max(all_measures),
            "mean": np.mean(all_measures),
            "stdvar": np.std(all_measures),
        }
    else:
        global_stats = {"min": 0, "max": 0, "mean": 0, "stdvar": 0}
        
    return stats_per_receiver, global_stats

def compute_snapshot_stats(snapshot_data_list):
    stats_per_snapshot = []
    
    global_min = float('inf')
    global_max = float('-inf')
    global_sum = 0.0
    global_sq_sum = 0.0
    global_count = 0

    for item in snapshot_data_list:
        # Use absolute values for pressure field stats
        m = np.abs(np.array(item["measures"]))
        if m.size == 0: continue
        
        # Local stats
        l_min = np.min(m)
        l_max = np.max(m)
        l_mean = np.mean(m)
        l_std = np.std(m)
        
        stats_per_snapshot.append({
            "snapshot_index": item["index"],
            "min": l_min,
            "max": l_max,
            "mean": l_mean,
            "stdvar": l_std
        })
        
        # Global accumulation
        if l_min < global_min: global_min = l_min
        if l_max > global_max: global_max = l_max
        global_sum += np.sum(m)
        global_sq_sum += np.sum(m**2)
        global_count += m.size

    if global_count > 0:
        global_mean = global_sum / global_count
        global_var = (global_sq_sum / global_count) - (global_mean**2)
        global_std = np.sqrt(global_var) if global_var > 0 else 0
    else:
        global_min, global_max, global_mean, global_std = 0, 0, 0, 0

    global_stats = {
        "min": global_min,
        "max": global_max,
        "mean": global_mean,
        "stdvar": global_std,
    }
    return stats_per_snapshot, global_stats


def plot_pressure_evolution(receivers_data, dt, output_dir):
    if not receivers_data: return
    n_measures = len(receivers_data[0]["measures"])
    time_axis = np.arange(n_measures) * dt
    for i, recv in enumerate(receivers_data):
        plt.figure()
        plt.plot(time_axis, recv["measures"])
        plt.title(f"Absolute Pressure Evolution - Receiver {i}")
        plt.xlabel("Time (s)")
        plt.ylabel("|Pressure|")
        plt.savefig(os.path.join(output_dir, f"pressure_recv_{i}.png"))
        plt.close()


def plot_fourier(receivers_data, dt, output_dir):
    if not receivers_data: return
    n_measures = len(receivers_data[0]["measures"])
    xf = fftfreq(n_measures, dt)[: n_measures // 2]
    for i, recv in enumerate(receivers_data):
        yf = fft(recv["measures"])
        plt.figure()
        plt.plot(xf, 2.0 / n_measures * np.abs(yf[0 : n_measures // 2]))
        plt.title(f"Fourier Transform (Abs) - Receiver {i}")
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Amplitude")
        plt.savefig(os.path.join(output_dir, f"fourier_recv_{i}.png"))
        plt.close()


def plot_propagation_speed(receivers_data, dt, src_coords, output_dir):
    if not receivers_data: return
    
    distances = []
    arrival_times = []
    
    src = np.array(src_coords)
    
    for i, recv in enumerate(receivers_data):
        rcv_coords = np.array(recv["coords"])
        dist = np.linalg.norm(rcv_coords - src)
        
        measures = np.array(recv["measures"])
        peak_idx = np.argmax(measures)
        t_peak = peak_idx * dt
        
        distances.append(dist)
        arrival_times.append(t_peak)
    
    distances = np.array(distances)
    arrival_times = np.array(arrival_times)
    
    slope, intercept, r_value, p_value, std_err = linregress(arrival_times, distances)
    speed = slope
    
    plt.figure()
    plt.scatter(arrival_times, distances, label='Receivers')
    plt.plot(arrival_times, intercept + slope * arrival_times, 'r', label=f'Fit: v={speed:.2f} m/s')
    plt.title(f"Propagation Speed Estimation")
    plt.xlabel("Arrival Time (s)")
    plt.ylabel("Distance (m)")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "propagation_speed.png"))
    plt.close()
    
    print(f"Estimated Propagation Speed: {speed:.4f} m/s (R^2={r_value**2:.4f})")
    with open(os.path.join(output_dir, "propagation_speed.txt"), "w") as f:
        f.write(f"Estimated Speed: {speed} m/s\n")
        f.write(f"Intercept: {intercept} m\n")
        f.write(f"R-squared: {r_value**2}\n")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        prog="Ad-Hoc Extended Statistics",
        description="Compute different statistics from ad-hoc workflow",
    )

    args_input = arg_parser.add_argument_group("inputs")
    args_input.add_argument("--sismos")
    args_input.add_argument("--snapshot_folder")
    args_input.add_argument("--params-file", help="Path to exported parameters file")
    
    args_input.add_argument("--dt", type=float, default=0.001)
    args_input.add_argument("--srcx", type=float, default=1010.0)
    args_input.add_argument("--srcy", type=float, default=1010.0)
    args_input.add_argument("--srcz", type=float, default=1010.0)
    args_input.add_argument("--snapshot-format", default="plain", choices=["plain", "bin"])

    args_actions = arg_parser.add_argument_group("actions")
    args_actions.add_argument("-s", "--stats", action="store_true")
    args_actions.add_argument("-ft", "--fourier-transform", action="store_true")
    args_actions.add_argument("-ps", "--propagation-speed", action="store_true")

    args = arg_parser.parse_args()

    params = {}
    if args.params_file and os.path.exists(args.params_file):
        params = load_params(args.params_file)
        print(f"Loaded parameters: {params}")

    dt = params.get("dt", args.dt)
    srcx = params.get("srcx", args.srcx)
    srcy = params.get("srcy", args.srcy)
    srcz = params.get("srcz", args.srcz)
    
    if not args.sismos and not args.snapshot_folder:
        raise Exception("you did not provide either sismos or snapshots")

    if args.sismos and not os.path.exists(args.sismos):
        raise Exception(f"couldn't find sismos at {args.sismos}")

    output_name = f"output_{datetime.now().strftime('%d_%m_%y_%H_%M_%S')}"
    os.makedirs(output_name, exist_ok=True)

    if args.sismos:
        n_recv, n_meas, data = load_sismos(args.sismos)

        plot_pressure_evolution(data, dt, output_name)

        if args.stats:
            r_stats, g_stats = compute_sismos_stats(data)
            with open(
                os.path.join(output_name, "receiver_stats.csv"), "w", newline=""
            ) as f:
                if r_stats:
                    writer = csv.DictWriter(f, fieldnames=r_stats[0].keys())
                    writer.writeheader()
                    writer.writerows(r_stats)
            with open(
                os.path.join(output_name, "global_sismos_stats.csv"), "w", newline=""
            ) as f:
                writer = csv.DictWriter(f, fieldnames=g_stats.keys())
                writer.writeheader()
                writer.writerow(g_stats)

        if args.fourier_transform:
            plot_fourier(data, dt, output_name)

        if args.propagation_speed:
            src_coords = [srcx, srcy, srcz]
            plot_propagation_speed(data, dt, src_coords, output_name)

    if args.snapshot_folder:
        if os.path.isdir(args.snapshot_folder):
            snaps_files = sorted(
                glob.glob(os.path.join(args.snapshot_folder, "*")), 
                key=lambda x: int(re.findall(r'\d+', os.path.basename(x))[-1]) if re.findall(r'\d+', os.path.basename(x)) else 0
            )
            
            if len(snaps_files) == 0:
                print(f"Warning: couldn't find any snapshot at {args.snapshot_folder}")
            else:
                snapshot_data_list = []
                for i, snap_path in enumerate(snaps_files):
                     try:
                         idx = int(re.findall(r'\d+', os.path.basename(snap_path))[-1])
                     except:
                         idx = i
                         
                     data = load_snapshot(snap_path, args.snapshot_format)
                     flat_data = data.flatten()
                     snapshot_data_list.append({"index": idx, "measures": flat_data})
                
                if args.stats:
                    s_stats, g_snaps_stats = compute_snapshot_stats(snapshot_data_list)
                    
                    with open(os.path.join(output_name, "snapshot_stats.csv"), "w", newline="") as f:
                        if s_stats:
                            writer = csv.DictWriter(f, fieldnames=s_stats[0].keys())
                            writer.writeheader()
                            writer.writerows(s_stats)
                            
                    with open(os.path.join(output_name, "global_snapshots_stats.csv"), "w", newline="") as f:
                        writer = csv.DictWriter(f, fieldnames=g_snaps_stats.keys())
                        writer.writeheader()
                        writer.writerow(g_snaps_stats)
                    
                    indices = [s["snapshot_index"] for s in s_stats]
                    means = [s["mean"] for s in s_stats]
                    maxs = [s["max"] for s in s_stats]
                    
                    plt.figure()
                    plt.plot(indices, means, label="Mean Pressure")
                    plt.plot(indices, maxs, label="Max Pressure")
                    plt.title("Snapshot Statistics Evolution")
                    plt.xlabel("Snapshot Index")
                    plt.ylabel("Pressure")
                    plt.legend()
                    plt.savefig(os.path.join(output_name, "snapshot_stats_evolution.png"))
                    plt.close()

        else:
            raise Exception(f"couldn't find directory {args.snapshot_folder}")
