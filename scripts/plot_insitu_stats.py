import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot in-situ stats")
    parser.add_argument("--folder", default="insitu_stats", help="Folder containing in-situ CSVs")
    args = parser.parse_args()

    folder = args.folder
    if not os.path.exists(folder):
        print(f"Folder {folder} does not exist.")
        exit(1)

    sismos_data_path = os.path.join(folder, "sismos_data.csv")
    if os.path.exists(sismos_data_path):
        df = pd.read_csv(sismos_data_path)
        time = df["time"]
        for col in df.columns:
            if col == "time": continue
            plt.figure()
            plt.plot(time, np.abs(df[col]))
            plt.title(f"Absolute Pressure Evolution - {col}")
            plt.xlabel("Time (s)")
            plt.ylabel("|Pressure|")
            plt.savefig(os.path.join(folder, f"pressure_{col}.png"))
            plt.close()
        print("Plotted pressure evolution.")

    fft_data_path = os.path.join(folder, "fft_data.csv")
    if os.path.exists(fft_data_path):
        df = pd.read_csv(fft_data_path)
        freq = df["freq"]
        for col in df.columns:
            if col == "freq": continue
            plt.figure()
            plt.plot(freq, df[col])
            plt.title(f"FFT - {col}")
            plt.xlabel("Frequency (Hz)")
            plt.ylabel("Amplitude")
            plt.savefig(os.path.join(folder, f"fourier_{col}.png"))
            plt.close()
        print("Plotted FFT.")

    snap_stats_path = os.path.join(folder, "snapshot_stats.csv")
    if os.path.exists(snap_stats_path):
        df = pd.read_csv(snap_stats_path)
        if not df.empty:
            plt.figure()
            plt.plot(df["index"], df["mean"], label="Mean Pressure")
            plt.plot(df["index"], df["max"], label="Max Pressure")
            # plt.plot(df["index"], df["min"], label="Min Pressure")
            plt.title("Snapshot Statistics Evolution")
            plt.xlabel("Snapshot Index")
            plt.ylabel("Pressure")
            plt.legend()
            plt.savefig(os.path.join(folder, "snapshot_stats_evolution.png"))
            plt.close()
            print("Plotted snapshot stats evolution.")
