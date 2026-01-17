import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import os
import nodeCalc
import binDecode
import numpy as np
import time

def get_global_limits(file_type):
    cwd = os.getcwd() 
    path = os.path.join(cwd, "snapshot")
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
        sys.exit(1)

    it, z_slice, f_type = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3].lower()

    t1 = time.perf_counter()
    snapshot_path = os.path.join(os.getcwd(), "snapshot")
    meta = nodeCalc.getMetaDim(it, f_type) 
    data = binDecode.openBin(it) if f_type == "bin" else nodeCalc.getSnapshotData(it)
    
    total_nodes = 0
    if data and meta:
        nx = meta['ex'] * meta['order'] + 1
        ny = meta['ey'] * meta['order'] + 1
        nz = meta['ez'] * meta['order'] + 1
        total_nodes = nx * ny * nz
        matrix = get_slice_data(data, z_slice, meta)
    
    t2 = time.perf_counter()

    if total_nodes > 0:
        print(f"{total_nodes},{t2-t1}")
    else:
        print(f"0,0")
    