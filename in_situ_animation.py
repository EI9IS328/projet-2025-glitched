import glob
import re
import struct

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

FILEMAGIC = 0xCAFEBABE
FILEMAGIC_FORMAT = "<I"
FILEMAGIC_SIZE = 4
HEADER_FORMAT = "<ifffiii"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
ROW_FORMAT = "<ffff"
ROW_SIZE = struct.calcsize(ROW_FORMAT)

SNAPSHOT_FOLDER_PATH = "build/snapshots"
OUT_FOLDER_PATH = "out"
files = glob.glob(SNAPSHOT_FOLDER_PATH + "/snapshot*.bin")

def file_name_to_timestep(file_name: str) -> int:
  match = re.search(r"snapshot(\d+)\.bin", file_name)
  return int(match.group(1))
sorted_files = sorted(files, key=file_name_to_timestep)

fig, ax = plt.subplots()
images = []
for i, file_name in enumerate(sorted_files):
  timestep = file_name_to_timestep(file_name)
  print(f"generating image for timestep {timestep}")
  with open(file_name, "rb") as f:
    # parse header
    sig = struct.unpack(FILEMAGIC_FORMAT, f.read(FILEMAGIC_SIZE))[0]
    assert(sig == FILEMAGIC), "snapshot file magic number doesn't match expected"

    sliced_dim, domain_size_x, domain_size_y, domain_size_z, nb_nodes_x, nb_nodes_y, nb_nodes_z = struct.unpack(HEADER_FORMAT, f.read(HEADER_SIZE))
    domain_size_xyz = domain_size_x, domain_size_y, domain_size_z
    nb_nodes_xyz = nb_nodes_x, nb_nodes_y, nb_nodes_z
    step_xyz = step_x, step_y, step_z = domain_size_x / (nb_nodes_x-1), domain_size_y / (nb_nodes_y-1), domain_size_z / (nb_nodes_z-1)

    slice_data = np.empty(nb_nodes_xyz, dtype=np.float32)


    row_bytes = f.read(ROW_SIZE)
    while len(row_bytes) != 0:
      x, y, z, value = struct.unpack(ROW_FORMAT, row_bytes)
      nx, ny, nz = int(x/step_x), int(y/step_y), int(z/step_z)
      slice_data[nx, ny, nz] = value
      row_bytes = f.read(ROW_SIZE)


  idxs = [slice(None, None), slice(None, None), slice(None, None)]
  idxs[sliced_dim] = int(domain_size_xyz[sliced_dim] / 2 / step_xyz[sliced_dim])

  im = ax.imshow(slice_data[*idxs], vmin=0, vmax=0.1, cmap="magma", animated=True)
  if i == 0:
      ax.imshow(slice_data[*idxs], vmin=0, vmax=0.1, cmap="magma")  # show an initial one first
  images.append([im])

ani = animation.ArtistAnimation(fig, images, interval=16, blit=True, repeat_delay=1000, repeat=False)
ani.save("movie.mp4")


# To save the animation using Pillow as a gif
# writer = animation.PillowWriter(fps=15,
#                                 metadata=dict(artist='Me'),
#                                 bitrate=1800)
# ani.save('in_situ.gif', writer=writer)
