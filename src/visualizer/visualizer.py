import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import nodeCalc
import binDecode

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

def heatmap(data_slice, iteration, z_slice, file_type):
    plt.figure(figsize=(10, 8))
    im = plt.imshow(data_slice, origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar(im, label='Pressure (Pa)')
    plt.title(f"Snapshot {iteration} ({file_type}) - Z-Slice {z_slice}")
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

    meta = nodeCalc.getMetaDim(it_arg, type_arg)
    if not meta:
        print(f"Error: Snapshot file snapshot{it_arg}.{type_arg} not found.")
        sys.exit(1)

    if type_arg == "bin":
        values = binDecode.openBin(it_arg)
    else:
        values = nodeCalc.getSnapshotData(it_arg)

    if values:
        matrix = get_slice_unified(values, z_arg, meta)
        heatmap(matrix, it_arg, z_arg, type_arg)