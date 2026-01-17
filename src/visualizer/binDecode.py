import struct
import os

# Header size for: ex, ey, ez, order
HEADER_FORMAT = "<iiii"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)

def openBin(iteration):
    file = f"snapshot{iteration}.bin"
    dir = os.path.dirname(__file__)
    path = os.path.normpath(os.path.join(dir, '..', '..', "snapshot", file))

    iterData = {}

    if not os.path.exists(path):
        print(f"File not found: {path}")
        return {}

    with open(path, "rb") as f:
        header_bytes = f.read(HEADER_SIZE)
        ex_val, ey_val, ez_val, order_val = struct.unpack(HEADER_FORMAT, header_bytes)
        

        data_bytes = f.read()
        num_floats = len(data_bytes) // 4
        values = struct.unpack("<" + "f" * num_floats, data_bytes)

        for global_idx, val in enumerate(values):
            iterData[global_idx] = val

    return iterData