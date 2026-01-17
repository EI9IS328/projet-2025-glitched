# FUnTiDES: Fast Unstructured Time Dynamic Equation Solver

**FUnTiDES** is a collection of simplified codes that represent real scientific applications. It serves as a standard tool for evaluating and comparing the performance of various high-performance computing (HPC) systems, particularly those used for scientific simulations.

---

## Included Applications

The current implementation includes two proxy applications for solving the 2nd-order acoustic wave equation in 2D and 3D:

- **SEM (Spectral Element Method)**
  A benchmark designed to simulate wave propagation using SEM, a Galerkin-based finite element method for solving partial differential equations (PDEs).

- **FD (Finite Difference Method)**
  A benchmark that uses finite-difference stencil operators to simulate wave propagation and solve PDEs.

A key feature of these proxy applications is their adaptability to different programming models and HPC architectures. They are also easy to build and run, making them accessible to both researchers and developers.

---

## Supported Programming Models

The SEM proxy currently supports:

- [Kokkos](https://kokkos.github.io/kokkos-core-wiki/) — for performance portability

> **Note**: Kokkos is included as a Git submodule and will be compiled automatically when enabled.

---

## Supported Data Containers

The current SEM proxy supports the following data container:

- `std::vector` (default for serial )

---

## Quick Start: Build and Run

### Step 1: Compile and Install

```sh
mkdir build
cd build
cmake ..
make install
```

By default, this builds the applications in sequential mode using `std::vector`.
Both SEM and FD applications are compiled.

### Step 2: Run Examples

```sh
# Run SEM simulation with 100 x 100 x 100 elements
./src/main/semproxy -ex 100

# Run FD simulation
./src/main/fdproxy
```

---

## CMake Options

The following options can be used to configure your build:

| Option                 | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `COMPILE_FD`           | Enable compilation of the FD proxy (default: ON)                            |
| `COMPILE_SEM`          | Enable compilation of the SEM proxy (default: ON)                           |
| `ENABLE_CUDA`          | Enable CUDA backend (used by Kokkos)                                        |
| `ENABLE_PYWRAP`        | Enable Python bindings via pybind11 (experimental)                          |
| `USE_KOKKOS`           | Enable Kokkos support (serial by default, CUDA/OpenMP with flags)           |
| `USE_VECTOR`           | Use `std::vector` for data arrays (enabled by default unless Kokkos is used)|

---
## Quick Start, Nix:
Use the given flake (eg. `nix develop`). It is made for compiling with Ninja instead of Make mostly for compilation speed.

Once in the shell you can than do (for development version)

```sh
mkdr -p build/debug
cd build/debug
cmake ../.. -DUSE_KOKKOS=ON -DUSE_VECTOR=OFF -DCOMPILE_FD=OFF -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug -G Ninja
cd ../..
cmake --build build/debug
```

---
## Quick Start, Cremi:
```
# clone the repo
git submodule update --init --recursive
export PATH="/usr/local/cuda-12/bin/:${PATH}"
mkdir build
cd build
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DENABLE_CUDA=ON -DUSE_KOKKOS=ON
make
```

---
## Benchmarking Export Modes

The `scripts/benchmark_insitu_vs_adhoc.py` script compares different snapshot export strategies:

| Mode | Description |
|------|-------------|
| `base` | No snapshot export (baseline) |
| `adhoc-plain` | Ad-hoc export in plain text format |
| `adhoc-bin` | Ad-hoc export in binary format |
| `insitu` | In-situ grayscale image export |
| `rgb` | In-situ RGB image export with colormap |

### Running Benchmarks

```sh
# Run all modes with default grid sizes (10, 20, 30, 40, 50)
python scripts/benchmark_insitu_vs_adhoc.py

# Run specific modes
python scripts/benchmark_insitu_vs_adhoc.py --modes base insitu rgb

# Custom grid sizes and multiple runs for averaging
python scripts/benchmark_insitu_vs_adhoc.py --sizes 10 20 30 --runs 3

# Specify binary path and output file
python scripts/benchmark_insitu_vs_adhoc.py --bin ./build/bin/semproxy --output results.csv
```

### Visualizing Results

```sh
# Display plots interactively
python scripts/plot_benchmark.py results.csv

# Save to image file
python scripts/plot_benchmark.py results.csv --output benchmark_results.png
```
