# MadVoro

**Massively distributed construction of 3D Voronoi diagrams.**

MadVoro provides the core 3D Voronoi tessellation engine used by the [RICH](https://github.com/maormizrachi/RICH) astrophysical simulation code. It can also be used as a standalone library in other projects.

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

```bash
./install_deps.sh
mkdir build && cd build
cmake .. -DMADVORO_DEPS_DIR=../deps \
         -DMADVORO_WITH_MPI=ON \
         -DMADVORO_BUILD_EXAMPLES=ON
make -j$(nproc)
```

## Using Inside RICH

When used as a submodule inside RICH (at `source/3D/tessellation/voronoi/`), the CMakeLists.txt is not used. RICH's own build system compiles the MadVoro source files directly via `GLOB_RECURSE`. External dependency include paths are provided by RICH's CMake configuration.

## License

See the [RICH repository](https://github.com/maormizrachi/RICH) for license information.
