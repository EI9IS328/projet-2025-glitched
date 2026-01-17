import os
import struct

def getMetaDim(iteration):
    folder = "../../snapshot"
    base_path = os.path.join(os.path.dirname(__file__), folder, f"snapshot{iteration}")
    
    if os.path.exists(base_path + ".bin"):
        with open(base_path + ".bin", "rb") as f:
            data = f.read(16)
            ex, ey, ez, order = struct.unpack("<iiii", data)
            return ex, ey, ez, order, True # True means it IS binary
            
    elif os.path.exists(base_path + ".txt"):
        with open(base_path + ".txt", "r") as f:
            line = f.readline().strip().split(",")
            return int(line[0]), int(line[1]), int(line[2]), int(line[3]), False # False means text
            
    return None

ex, ey, ez, order, is_bin = getMetaDim(0)


def getSnapshotData(iter):
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
    
    iterData = {}
    iter = str(iter)
    file = "snapshot" + iter + ".bin"
    dir = os.path.dirname(__file__)
    path = os.path.normpath(os.path.join(dir, '..', '..', "snapshot", file))
    with open(file=path, newline='') as snapshot:
        next(snapshot)
        el = 0
        for row in snapshot:
            nodeIdx = 0
            nodes = row.split(",")
            for nodePrs in nodes:
                nodePrs = nodePrs.replace("\n", "")
                nodeIdx += 1
                x = nodeIdx % (order + 1)
                y = (nodeIdx // (order + 1)) % (order + 1)
                z = nodeIdx // ( (order + 1) * (order + 1)) 
                iterData[globalNodeIndex(el,x,y,z)] = nodePrs
            el+=1
    return iterData



def nodeCoord(globalIdx):
    nodesPerDim = [ex * order + 1 , ey * order + 1, ez * order + 1]
    
    k = globalIdx // (nodesPerDim[0] * nodesPerDim[1])
    remainder = globalIdx % (nodesPerDim[0] * nodesPerDim[1])
    j = remainder // nodesPerDim[0]
    i = remainder % nodesPerDim[0]

    nodeIdx = [i, j, k]

    return nodeIdx

def nodeCoordBin(globalIdx, ex, ey, ez, order):
    nx = ex * order + 1
    ny = ey * order + 1
    nz = ez * order + 1


    k = globalIdx // (nx * ny)
    remainder = globalIdx % (nx * ny)
    j = remainder // nx
    i = remainder % nx

    return [i, j, k]