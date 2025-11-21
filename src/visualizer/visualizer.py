import matplotlib.pyplot as plt
import matplotlib.colors as colors
import nodeCalc

def getSlice(nodes, z):
    prsMatrix = [(nodeCalc.ey * nodeCalc.order + 1) * [0] \
                  for i in range(nodeCalc.ex * nodeCalc.order + 1)]
    for key, value in nodes.items():
        nodeCoord = nodeCalc.nodeCoord(key)
        if nodeCoord[2] == z:
            prsMatrix[nodeCoord[0]][nodeCoord[1]] = float(value)
    return prsMatrix


def heatmap(slice):
    colors_list = ['#0099ff', '#33cc33']
    cmap = colors.ListedColormap(colors_list)
    plt.imshow(slice, \
               norm=colors.SymLogNorm(linthresh=1e-10, \
               linscale=1, vmin=-1e20, vmax=1e20))
    plt.colorbar()

    plt.title("Customized heatmap with annotations")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.show()

if __name__ == "__main__":
    values = nodeCalc.getSnapshotData(500)
    matrix = getSlice(values, 1)
    heatmap(matrix)

