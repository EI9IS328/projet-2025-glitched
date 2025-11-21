import os

ex = 2
ey = 1
ez = 1
order = 2


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


def getSnapshotData(iter):
    iterData = {}
    iter = str(iter)
    file = "snapshot" + iter + ".bin"
    dir = os.path.dirname(__file__)
    path = os.path.normpath(os.path.join(dir, '..', '..', "snapshot", file))
    with open(file=path, newline='') as snapshot:
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

