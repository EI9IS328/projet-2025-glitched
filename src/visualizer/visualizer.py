import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import os
import nodeCalc
import binDecode

def get_global_limits(file_type):
    folder = "../../snapshot"
    path = os.path.join(os.path.dirname(__file__), folder)
    ext = ".bin" if file_type == "bin" else ".txt"
    
    files = [f for f in os.listdir(path) if f.startswith("snapshot") and f.endswith(ext)]

    global_min = float('inf')
    global_max = float('-inf')
    
    for filename in files:
        try:
            it_str = filename.replace("snapshot", "").replace(ext, "")
            iteration = int(it_str)
            
            if file_type == "bin":
                data = binDecode.openBin(iteration)
            else:
                data = nodeCalc.getSnapshotData(iteration)
            
            if data:
                vals = [float(v) for v in data.values()]
                global_min = min(global_min, min(vals))
                global_max = max(global_max, max(vals))
        except ValueError:
            continue

    print(f"[+] Global Bounds Found: vmin={global_min}, vmax={global_max}")
    return global_min, global_max

def get_slice_unified(nodes, z_target, meta):
    nx = meta['ex'] * meta['order'] + 1
    ny = meta['ey'] * meta['order'] + 1
    prsMatrix = [[0.0 for _ in range(ny)] for _ in range(nx)]
    
    for key, value in nodes.items():
        coord = nodeCalc.get_node_coords(key, meta)
        if coord[2] == z_target:
            if 0 <= coord[0] < nx and 0 <= coord[1] < ny:
                prsMatrix[coord[0]][coord[1]] = float(value)
    return prsMatrix

def heatmap(data_slice, iteration, z_slice, file_type, vmin, vmax):
    plt.figure(figsize=(10, 8))
    norm = colors.SymLogNorm(linthresh=1e-10, linscale=1, vmin=vmin, vmax=vmax)
    
    im = plt.imshow(data_slice, origin='lower', cmap='viridis', aspect='auto', norm=norm)
    plt.colorbar(im, label='Pressure (Pa)')
    plt.title(f"Snapshot {iteration} ({file_type}) | Z={z_slice}\nGlobal Scale: [{vmin:.2e}, {vmax:.2e}]")
    plt.xlabel("X Node Index")
    plt.ylabel("Y Node Index")
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python visualizer.py <iteration> <z_slice> <bin|txt>")
        sys.exit(1)

    it_arg = int(sys.argv[1])
    z_arg = int(sys.argv[2])
    type_arg = sys.argv[3].lower()

    vmin, vmax = get_global_limits(type_arg)
    
    if vmin is None or vmax is None:
        sys.exit(1)

    meta = nodeCalc.getMetaDim(it_arg, type_arg)
    if meta:
        values = binDecode.openBin(it_arg) if type_arg == "bin" else nodeCalc.getSnapshotData(it_arg)
        
        if values:
            matrix = get_slice_unified(values, z_arg, meta)
            heatmap(matrix, it_arg, z_arg, type_arg, vmin, vmax)