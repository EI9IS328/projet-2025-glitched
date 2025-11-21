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

    min_val = min(np.min(np.array(getSlice(nodes, Z))) for Z in Zlevels)
    max_val = max(np.max(np.array(getSlice(nodes, Z))) for Z in Zlevels)

    norm = colors.SymLogNorm(linthresh=1e-10, linscale=1, vmin=min_val, vmax=max_val)
    cmap = plt.cm.viridis

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for z in Zlevels:
        Pi = getSlice(nodes, z)
        Pi = np.array(Pi)
        Zg = np.full_like(Pi, z)
        facecolors = cmap(norm(Pi))

        ax.plot_surface(
            Xg, Yg, Zg,
            facecolors=facecolors,
            rstride=1,
            cstride=1,
            shade=False  
        )

    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])  # required
    plt.colorbar(mappable, ax=ax, shrink=0.5, pad=0.1)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

if __name__ == "__main__":
    values = nodeCalc.getSnapshotData(600)
    surface3d(values)
    # matrix = getSlice(values, 1)
    # heatmap(matrix)

