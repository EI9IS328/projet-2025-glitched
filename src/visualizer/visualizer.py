import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import os
import nodeCalc
import binDecode
import numpy as np

def get_global_limits(file_type):
    folder = "../../snapshot"
    path = os.path.normpath(os.path.join(os.path.dirname(__file__), folder))
    ext = f".{file_type}"
    files = [f for f in os.listdir(path) if f.startswith("snapshot") and f.endswith(ext)]
    
    g_min, g_max = float('inf'), float('-inf')

    for f_name in files:
        try:
            it = int(f_name.replace("snapshot", "").replace(ext, ""))
            if file_type == "bin":
                data = binDecode.openBin(it)
            else:
                data = nodeCalc.getSnapshotData(it)
            
            if data:
                vals = [float(v) for v in data.values()]
                g_min, g_max = min(g_min, min(vals)), max(g_max, max(vals))
        except: continue
    return g_min, g_max

def get_slice_data(data, z_target, meta):
    nx = meta['ex'] * meta['order'] + 1
    ny = meta['ey'] * meta['order'] + 1
    nz = meta['ez'] * meta['order'] + 1
    num_nodes = nx * ny * nz 
    
    matrix = np.zeros((nx, ny))
    for idx in range(num_nodes):
        if idx not in data: break
        coords = nodeCalc.get_node_coords(idx, meta)
        if coords[2] == z_target:
            matrix[coords[0], coords[1]] = float(data[idx])
    return matrix

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 visualizer.py <iteration> <z_slice> <bin|txt>")
        sys.exit(1)

    it, z_slice, f_type = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3].lower()

    vmin, vmax = get_global_limits(f_type)
    
    meta = nodeCalc.getMetaDim(it, f_type)
    data = binDecode.openBin(it) if f_type == "bin" else nodeCalc.getSnapshotData(it)
    
    if data and meta:
        matrix = get_slice_data(data, z_slice, meta)
        
        plt.figure(figsize=(10, 8))
        norm = colors.SymLogNorm(linthresh=1e-10, linscale=1, vmin=vmin, vmax=vmax)
        im = plt.imshow(matrix.T, origin='lower', cmap='viridis', aspect='auto', norm=norm)
        plt.colorbar(im, label='Pressure (Pa)')
        plt.title(f"Iteration {it} | Z-Slice {z_slice} | Fixed Scale")
        plt.xlabel("X Index")
        plt.ylabel("Y Index")
        plt.show()