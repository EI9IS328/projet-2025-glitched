import glob
import re
import struct
import os

FILEMAGIC = 0xCAFEBABE
FILEMAGIC_FORMAT = "<I"
FILEMAGIC_SIZE = 4
HEADER_FORMAT = "<iiii"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
ROW_FORMAT = "<" + "f" * 27
ROW_SIZE = struct.calcsize(ROW_FORMAT)

def openBin():
    file = "snapshot600.bin"
    dir = os.path.dirname(__file__)
    path = os.path.normpath(os.path.join(dir, '..', '..', "snapshot", file))

    iterData = {}

    with open(path, "rb") as f:
        ex, ey, ez, order = struct.unpack(HEADER_FORMAT, f.read(HEADER_SIZE))
        
        nx = order * ex + 1
        ny = order * ey + 1
        nz = order * ez + 1

        def globalNodeIndex(el, i, j, k):
            elemZ = el // (ex * ey)
            tmp = el % (ex * ey)
            elemY = tmp // ex
            elemX = tmp % ex

            ix = elemX * order + i
            iy = elemY * order + j
            iz = elemZ * order + k

            return ix + iy * nx + iz * nx * ny

        row_bytes = f.read(ROW_SIZE)
        el = 0
        while len(row_bytes) != 0:
            nodeIdx = 0
            nodes = struct.unpack(ROW_FORMAT, row_bytes)

            for nodePrs in nodes:
                nodeIdx +=1
                x = nodeIdx % (order + 1)
                y = (nodeIdx // (order + 1) % (order + 1))
                z = nodeIdx // (order + 1) * (order + 1)
                iterData[globalNodeIndex(el,x,y,z)] = nodePrs

            el += 1
            row_bytes = f.read(ROW_SIZE)
        
        return iterData



if __name__ == "__main__":
    print(openBin())