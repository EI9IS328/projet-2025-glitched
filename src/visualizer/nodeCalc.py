import os
import struct

def getMetaDim(iteration, file_type):
    folder = "../../snapshot"
    ext = ".bin" if file_type == "bin" else ".txt"
    path = os.path.normpath(os.path.join(os.path.dirname(__file__), folder, f"snapshot{iteration}{ext}"))
    
    if not os.path.exists(path):
        return None

    if file_type == "bin":
        with open(path, "rb") as f:
            data = f.read(16)
            if len(data) < 16: return None
            ex, ey, ez, order = struct.unpack("<iiii", data)
            return {"ex": ex, "ey": ey, "ez": ez, "order": order, "is_bin": True}
    else:
        with open(path, "r") as f:
            line = f.readline().strip().split(",")
            if not line or len(line) < 4: return None
            return {"ex": int(line[0]), "ey": int(line[1]), "ez": int(line[2]), "order": int(line[3]), "is_bin": False}

def get_node_coords(idx, meta):
    nx = meta['ex'] * meta['order'] + 1
    ny = meta['ey'] * meta['order'] + 1
    
    z = idx // (nx * ny)
    remainder = idx % (nx * ny)
    y = remainder // nx
    x = remainder % nx
    return [x, y, z]

def getSnapshotData(iteration):
    folder = "../../snapshot"
    path = os.path.normpath(os.path.join(os.path.dirname(__file__), folder, f"snapshot{iteration}.txt"))
    meta = getMetaDim(iteration, "txt")
    if not meta: return {}
    
    nx, ny = meta['ex'] * meta['order'] + 1, meta['ey'] * meta['order'] + 1
    order = meta['order']
    
    iterData = {}
    with open(path, 'r') as f:
        next(f) # Skip header
        for el, row in enumerate(f):
            nodes = row.strip().split(",")
            for nodeIdx, val in enumerate(nodes):
                if not val: continue
                lx, ly = nodeIdx % (order + 1), (nodeIdx // (order + 1)) % (order + 1)
                lz = nodeIdx // ((order + 1)**2)
                elemZ = el // (meta['ex'] * meta['ey'])
                tmp = el % (meta['ex'] * meta['ey'])
                elemY, elemX = tmp // meta['ex'], tmp % meta['ex']
                gx, gy, gz = elemX * order + lx, elemY * order + ly, elemZ * order + lz
                gIdx = gx + gy * nx + gz * nx * ny
                iterData[gIdx] = val
    return iterData