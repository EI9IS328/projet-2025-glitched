import glob
import re

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

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
  with open(file_name) as f:
    # parse header
    sliced_dim = int(f.readline())
    domain_size_xyz = domain_size_x, domain_size_y, domain_size_z = [int(n) for n in f.readline().split(",")]
    nb_nodes_xyz = nb_nodes_x, nb_nodes_y, nb_nodes_z = [int(n) for n in f.readline().split(",")]
    step_xyz = step_x, step_y, step_z = domain_size_x / (nb_nodes_x-1), domain_size_y / (nb_nodes_y-1), domain_size_z / (nb_nodes_z-1)

    slice_data = np.empty(nb_nodes_xyz, dtype=np.float32)

    line = f.readline()
    while len(line) != 0:
      x, y, z, value = line.split(",")
      x, y, z, value = int(x), int(y), int(z), float(value)
      nx, ny, nz = int(x/step_x), int(y/step_y), int(z/step_z)
      slice_data[nx, ny, nz] = value
      line = f.readline()


  idxs = [slice(None, None), slice(None, None), slice(None, None)]
  idxs[sliced_dim] = int(domain_size_xyz[sliced_dim] / 2 / step_xyz[sliced_dim])

  im = ax.imshow(slice_data[*idxs], vmin=0, vmax=0.1, animated=True)
  if i == 0:
      ax.imshow(slice_data[*idxs])  # show an initial one first
  images.append([im])

ani = animation.ArtistAnimation(fig, images, interval=50, blit=True, repeat_delay=1000, repeat=False)
ani.save("movie.mp4")


# To save the animation using Pillow as a gif
# writer = animation.PillowWriter(fps=15,
#                                 metadata=dict(artist='Me'),
#                                 bitrate=1800)
# ani.save('in_situ.gif', writer=writer)
