# MadVoro - Massively distributed Construction of Voronoi Diagrams
![C++ project](https://img.shields.io/badge/C++-2874a6)
![Linux](https://img.shields.io/badge/Linux-0e6655)
![macOS](https://img.shields.io/badge/macOS-27ae60)

MadVoro provides the core 3D Voronoi tessellation engine used by the [RICH](https://github.com/maormizrachi/RICH) astrophysical simulation code. It can also be used as a standalone library in other projects.

<img src="examples/fox/fox.png?raw=true" alt="An example for a fox mesh, in a 16 processors construction." width="500"/><img src="examples/pyramid/pyramid.png?raw=true" alt="An example for a mesh construction of a pyramidal space." width="300"/>

## Directory Structure

```
MadVoro/
├── Voronoi3D.hpp/.cpp        Core Voronoi tessellation engine
├── VoronoiPayload.hpp        Per-point payload for ghost exchange
├── delaunay/                  Delaunay triangulation
├── geometry/                  Face-sphere intersection tests
├── utils/                     Exact geometric predicates (Shewchuk)
├── exception/                 Error types
├── io/                        I/O utilities (HDF5, VTK, binary)
├── examples/                  Example programs
├── install_deps.sh            Script to clone external dependencies
└── CMakeLists.txt             Standalone build configuration
```

## External Dependencies

MadVoro depends on several external libraries. When used inside RICH, these are already present as submodules. For standalone usage, run:

```bash
./install_deps.sh
```

This clones into `deps/`:
- **[mpi_utils](https://github.com/maormizrachi/mpi_utils)** — MPI serialization, exchange, collectives
- **[spatial_ds](https://github.com/maormizrachi/spatial_ds)** — Spatial data structures (OctTree, KDTree, RangeTree)
- **[MeshDecomposer3D](https://github.com/maormizrachi/MeshDecomposer3D)** — Domain decomposition, Hilbert ordering, load balancing

Additionally, you need:
- **Boost** ≥ 1.74 (container, multiprecision)
- **OpenMP**
- **MPI** (optional, for parallel builds)
- **HDF5** (optional, for HDF5 I/O)
- **VTK** ≥ 9.3 (optional, for VTK output)

## Building Standalone

Run CMake from the MadVoro repository root. The standalone dependency directory
defaults to `deps/`, so `MADVORO_DEPS_DIR` can usually be omitted after running
`install_deps.sh`. Standalone single-config builds default to `Release`; pass
`-DCMAKE_BUILD_TYPE=Debug` or `-DCMAKE_BUILD_TYPE=RelWithDebInfo` when you want
debug symbols or debug behavior.

### Serial Library and Examples

```bash
git clone git@github.com:maormizrachi/MadVoro.git
cd MadVoro
./install_deps.sh
cmake -S . -B build \
      -DMADVORO_BUILD_EXAMPLES=ON
cmake --build build -j"$(nproc)"
```

### MPI Library and Examples

```bash
./install_deps.sh
cmake -S . -B build \
      -DMADVORO_WITH_MPI=ON \
      -DMADVORO_BUILD_EXAMPLES=ON
cmake --build build -j"$(nproc)"
```

Example executables are written under per-example build directories, such as:

```text
build/examples/uniform_points/example_uniform_parallel
build/examples/fox/example_fox
```

For convenience, CMake also copies each built example executable back beside its
source files, such as `examples/fox/example_fox`. These generated files are
ignored by git.

Example `input/` directories are copied beside the matching executable. For
example, the fox input files are copied to `build/examples/fox/input/`.

### Building with VTK

Enable VTK output with `MADVORO_WITH_VTK=ON`:

```bash
cmake -S . -B build-vtk \
      -DMADVORO_WITH_MPI=ON \
      -DMADVORO_WITH_VTK=ON \
      -DMADVORO_BUILD_EXAMPLES=ON
cmake --build build-vtk -j"$(nproc)"
```

MadVoro requires VTK 9.3 or newer. With MPI enabled, CMake also requires VTK's
parallel MPI components. If CMake cannot find VTK automatically, point it at
your VTK package configuration:

```bash
cmake -S . -B build-vtk \
      -DMADVORO_WITH_MPI=ON \
      -DMADVORO_WITH_VTK=ON \
      -DVTK_DIR=/path/to/vtk/lib/cmake/vtk-9.3 \
      -DMADVORO_BUILD_EXAMPLES=ON
```

You can also use `-DCMAKE_PREFIX_PATH=/path/to/vtk` if that is how VTK is
installed on your system.

### Building with HDF5

Enable HDF5 I/O with `MADVORO_WITH_HDF5=ON`:

```bash
cmake -S . -B build-hdf5 \
      -DMADVORO_WITH_MPI=ON \
      -DMADVORO_WITH_HDF5=ON \
      -DMADVORO_BUILD_EXAMPLES=ON
cmake --build build-hdf5 -j"$(nproc)"
```

When `MADVORO_WITH_MPI=ON`, MadVoro asks CMake to prefer a parallel HDF5 build.
If CMake cannot find the right HDF5 installation, provide one of:

```bash
-DHDF5_ROOT=/path/to/hdf5
```

or:

```bash
-DCMAKE_PREFIX_PATH=/path/to/hdf5
```

### Building with Both VTK and HDF5

```bash
cmake -S . -B build-full \
      -DMADVORO_WITH_MPI=ON \
      -DMADVORO_WITH_VTK=ON \
      -DMADVORO_WITH_HDF5=ON \
      -DVTK_DIR=/path/to/vtk/lib/cmake/vtk-9.3 \
      -DHDF5_ROOT=/path/to/hdf5 \
      -DMADVORO_BUILD_EXAMPLES=ON
cmake --build build-full -j"$(nproc)"
```

## Using Inside RICH

When used as a submodule inside RICH (at `source/3D/tessellation/voronoi/`), the CMakeLists.txt is not used. RICH's own build system compiles the MadVoro source files directly via `GLOB_RECURSE`. External dependency include paths are provided by RICH's CMake configuration.

## License

BSD 3-Clause. See [LICENSE](LICENSE) for details.
