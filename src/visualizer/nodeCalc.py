import os
import struct

def getMetaDim(iteration, file_type):
    folder = "../../snapshot"
    ext = ".bin" if file_type == "bin" else ".txt"
    path = os.path.join(os.path.dirname(__file__), folder, f"snapshot{iteration}{ext}")
    
    if not os.path.exists(path):
        return None

    if file_type == "bin":
        with open(path, "rb") as f:
            data = f.read(16)
            ex, ey, ez, order = struct.unpack("<iiii", data)
            return {"ex": ex, "ey": ey, "ez": ez, "order": order, "is_bin": True}
    else:
        with open(path, "r") as f:
            line = f.readline().strip().split(",")
            return {"ex": int(line[0]), "ey": int(line[1]), "ez": int(line[2]), "order": int(line[3]), "is_bin": False}

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