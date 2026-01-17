import struct
import os

def openBin(iteration):
    folder = "../../snapshot"
    path = os.path.normpath(os.path.join(os.path.dirname(__file__), folder, f"snapshot{iteration}.bin"))
    
    iterData = {}
    if not os.path.exists(path):
        return {}

    with open(path, "rb") as f:
        f.seek(16) 
        data_bytes = f.read()
        num_floats = len(data_bytes) // 4
        values = struct.unpack("<" + "f" * num_floats, data_bytes)
        
        for idx, val in enumerate(values):
            iterData[idx] = val
    return iterData