import subprocess
import time
import csv
import os
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
PROXY_EXE = "./bin/semproxy"
VISUALIZER_SCRIPT = "../src/visualizer/visualizerBench.py"
SNAPSHOT_DIR = "snapshot"
ITERATION = 0
Z_SLICE = 10

DIMENSIONS = [10, 20, 40, 60, 80]
FORMATS = ["bin", "plain"] 
REPETITIONS = 3

def run_benchmark():
    all_results = {f_type: [] for f_type in FORMATS}

    if not os.path.exists(SNAPSHOT_DIR):
        os.makedirs(SNAPSHOT_DIR)

    print(f"{'Format':>7} | {'Dim':>5} | {'Avg Time (s)':>15} | {'Status'}")
    print("-" * 45)

    for f_type in FORMATS:
        for dim in DIMENSIONS:
            durations = []
            node_counts = []
            

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
                    output = process.stdout.strip().split('\n')[-1]
                    nodes, exec_time = output.split(',')
                    node_counts.append(int(nodes))
                    durations.append(float(exec_time))
                except:
                    continue

            if durations:
                avg_time = sum(durations) / len(durations)
                actual_nodes = node_counts[0] 
                all_results[f_type].append((actual_nodes, avg_time))
                print(f"{f_type:>7} | {actual_nodes:8d} nodes | {avg_time:15.6f}s | Success")
    
    return all_results

def save_and_plot(all_results):
    plt.figure(figsize=(10, 6))
    
    for f_type, data in all_results.items():

        nodes, times = zip(*data)
        plt.plot(nodes, times, marker='o', label=f"Format: {f_type}")

    plt.title("Performance vs Total Number of Nodes")
    plt.xlabel("Total Nodes ($N_{total}$)")
    plt.ylabel("Execution Time (s)")
    plt.xscale('log') 
    plt.yscale('log') 
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.savefig("benchmark_nodes.png")

if __name__ == "__main__":
    results = run_benchmark()
    save_and_plot(results)