import subprocess
import time
import csv
import os
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
PROXY_EXE = "./bin/semproxy"
VISUALIZER_SCRIPT = "visualizer.py"
SNAPSHOT_DIR = "snapshot"
ITERATION = 0
Z_SLICE = 10

DIMENSIONS = [10, 20, 40, 60, 80]
FORMATS = ["bin", "plain"] 
REPETITIONS = 2 

def run_benchmark():
    all_results = {f_type: [] for f_type in FORMATS}

    if not os.path.exists(SNAPSHOT_DIR):
        os.makedirs(SNAPSHOT_DIR)

    print(f"{'Format':>7} | {'Dim':>5} | {'Avg Time (s)':>15} | {'Status'}")
    print("-" * 45)

    for f_type in FORMATS:
        for dim in DIMENSIONS:
            durations = []
            

            gen_cmd = [
                PROXY_EXE, 
                "--ex", str(dim), "--ey", str(dim), "--ez", str(dim), 
                "--snapshot-folder-path", SNAPSHOT_DIR,
                "--snapshot-format", f_type
            ]
            
            subprocess.run(gen_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            for _ in range(REPETITIONS):
                viz_cmd = ["python3", VISUALIZER_SCRIPT, str(ITERATION), str(Z_SLICE), f_type]
                process = subprocess.run(viz_cmd, capture_output=True, text=True)
                
                try:
                    exec_time = float(process.stdout.strip().split('\n')[-1])
                    durations.append(exec_time)
                except:
                    continue

            if durations:
                avg_time = sum(durations) / len(durations)
                all_results[f_type].append((dim, avg_time))
                print(f"{f_type:>7} | {dim:5d} | {avg_time:15.6f}s | Success")

            for f in os.listdir(SNAPSHOT_DIR):
                os.remove(os.path.join(SNAPSHOT_DIR, f))

    return all_results

def save_and_plot(all_results):
    plt.figure(figsize=(10, 6))
    
    with open("benchmark_comparison.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Format", "Dimension", "Time_Seconds"])

        for f_type, data in all_results.items():
            if not data: continue
            dims, times = zip(*data)
            
            for d, t in data:
                writer.writerow([f_type, d, t])
            
            plt.plot(dims, times, marker='o', label=f"Format: {f_type}")

    plt.title("Data Processing Speed: Binary vs Plain Text")
    plt.xlabel("Grid Dimension (N)")
    plt.ylabel("Time (seconds)")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    
    
    
    plt.savefig("benchmark_comparison.png")
    print(f"\nBenchmark finished. Results saved in 'benchmark_comparison.png'")

if __name__ == "__main__":
    results = run_benchmark()
    save_and_plot(results)