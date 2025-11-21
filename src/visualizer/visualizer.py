import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
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

    min_val = min(min(row) for row in slice)
    max_val = max(max(row) for row in slice)

    plt.imshow(slice, \
               norm=colors.SymLogNorm(linthresh=1e-10, \
               linscale=1, vmin=min_val, vmax=max_val))
    plt.colorbar()

    plt.title("Customized heatmap with annotations")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.show()

def surface3d(nodes):
    Zlevels = [i for i in range(nodeCalc.ez * nodeCalc.order + 1)]

    nx = np.arange(nodeCalc.order * nodeCalc.ex + 1)
    ny = np.arange(nodeCalc.order * nodeCalc.ey + 1)
    Xg, Yg = np.meshgrid(nx, ny)



    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for z in Zlevels:
        Pi = getSlice(nodes, z)
        Pi = np.array(Pi)

        ax.contourf(
            Xg, Yg, np.full_like(Pi, z), 
            Pi,
            zdir='z',
            offset=z,
            levels=20,
            cmap='viridis'
        )

    plt.show()

if __name__ == "__main__":
    values = nodeCalc.getSnapshotData(1200)
    surface3d(values)
    # matrix = getSlice(values, 1)
    # heatmap(matrix)

