import matplotlib.pyplot as plt
import nodeCalc

def getSlice(nodes, z):
    prsMatrix = [(nodeCalc.ey * nodeCalc.order + 1) * [0] for i in range(nodeCalc.ex * nodeCalc.order + 1)]
    for key, value in nodes.items():
        nodeCoord = nodeCalc.nodeCoord(key)
        if nodeCoord[2] == z:
            prsMatrix[nodeCoord[0]][nodeCoord[1]] = float(value)
    return prsMatrix


def heatmap(slice):

    pass

if __name__ == "__main__":
    values = nodeCalc.getSnapshotData(400)
    matrix = getSlice(values, 0)