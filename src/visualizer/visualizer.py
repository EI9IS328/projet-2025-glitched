import matplotlib.pyplot as plt
import matplotlib.colors as colors
import nodeCalc
import binDecode

def get_slice_data(nodes, z_target, iteration):
    meta = nodeCalc.getMetaDim(iteration)
    if not meta:
        raise FileNotFoundError(f"Snapshot {iteration} not found.")
    
    nx = meta['ex'] * meta['order'] + 1
    ny = meta['ey'] * meta['order'] + 1
    
    prsMatrix = [[0.0 for _ in range(ny)] for _ in range(nx)]
    
    for key, value in nodes.items():
        coord = nodeCalc.get_node_coords(key, meta)
        
        if coord[2] == z_target:
            if 0 <= coord[0] < nx and 0 <= coord[1] < ny:
                prsMatrix[coord[0]][coord[1]] = float(value)
            
    return prsMatrix

def heatmap(data_slice):
    plt.figure(figsize=(10, 8))
    im = plt.imshow(data_slice, origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar(im, label='Pressure (Pa)')
    plt.xlabel("X-axis (Nodes)")
    plt.ylabel("Y-axis (Nodes)")
    plt.title("Pressure Slice")
    plt.show()

if __name__ == "__main__":
    ITERATION = 1200
    Z_PLANE = 25

    meta = nodeCalc.getMetaDim(ITERATION)
    if meta and meta['is_bin']:
        data = binDecode.openBin(ITERATION)
    else:
        data = nodeCalc.getSnapshotData(ITERATION) 

    if data:
        matrix = get_slice_data(data, Z_PLANE, ITERATION)
        heatmap(matrix)