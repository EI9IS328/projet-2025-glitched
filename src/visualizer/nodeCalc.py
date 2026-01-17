import os
import struct

def getMetaDim(iteration):
    folder = "../../snapshot"
    base_path = os.path.join(os.path.dirname(__file__), folder, f"snapshot{iteration}")
    
    if os.path.exists(base_path + ".bin"):
        with open(base_path + ".bin", "rb") as f:
            data = f.read(16)
            ex, ey, ez, order = struct.unpack("<iiii", data)
            return {"ex": ex, "ey": ey, "ez": ez, "order": order, "is_bin": True}
            
    elif os.path.exists(base_path + ".txt"):
        with open(base_path + ".txt", "r") as f:
            line = f.readline().strip().split(",")
            return {"ex": int(line[0]), "ey": int(line[1]), "ez": int(line[2]), "order": int(line[3]), "is_bin": False}
            
    return None

def get_node_coords(globalIdx, meta):
    nx = meta['ex'] * meta['order'] + 1
    ny = meta['ey'] * meta['order'] + 1
    nz = meta['ez'] * meta['order'] + 1

    if meta['is_bin']:
        i = globalIdx // (ny * nz)
        remainder = globalIdx % (ny * nz)
        j = remainder // nz
        k = remainder % nz
    else:
        k = globalIdx // (nx * ny)
        remainder = globalIdx % (nx * ny)
        j = remainder // nx
        i = remainder % nx

    return [i, j, k]