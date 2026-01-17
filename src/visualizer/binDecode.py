import struct
import os

def openBin(iteration):
    folder = "../../snapshot"
    path = os.path.normpath(os.path.join(os.path.dirname(__file__), folder, f"snapshot{iteration}.bin"))
    if not os.path.exists(path): return {}

    with open(path, "rb") as f:
        f.seek(16)
        data = f.read()
        num_floats = len(data) // 4
        values = struct.unpack(f"<{num_floats}f", data)
        return {i: v for i, v in enumerate(values)}