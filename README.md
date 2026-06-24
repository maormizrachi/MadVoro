# MadVoro - Massively Distributed Construction of Voronoi Diagrams
![C++ project](https://img.shields.io/badge/C++-2874a6)
![Linux](https://img.shields.io/badge/Linux-0e6655)
![macOS](https://img.shields.io/badge/macOS-27ae60)

MadVoro is a C++ framework for construction of Voronoi diagrams of 3D points, in distributed memory (using MPI).
MadVoro is a free and open-source project, released under the **Creative Commons Attribution 4.0 International (CC-BY 4.0)** license.  
You are free to use, modify, and distribute the code for any purpose — **provided that you cite the following paper** (an available BibTeX citation appears by the end of this file).

<img src="examples/fox/fox.png?raw=true" alt="An example for a fox mesh, in a 16 processors construction." width="500"/><img src="examples/pyramid/pyramid.png?raw=true" alt="An example for a mesh construction of a pyramidal space." width="300"/>


## Requirements
- CMake >= 3.14
- Any C++ compiler with C++17 support (tested with `g++` and `icpx`)
- Boost >= 1.74.0

### Optional
- For parallel support: Any MPI implementation (tested with OpenMPI, IntelMPI and MPICH).
- Recommended: [VCL](https://github.com/vectorclass/version2) for vectorization acceleration.
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) (C++ binding) >= 1.8.0, for output.
- [VTK](https://vtk.org/) >= 9.3.0, for visualization.

## Quick Start

The easiest way to use MadVoro is to add it directly to your CMake project as a subdirectory. No separate build or install step is needed — MadVoro's sources are compiled as part of your project.

**1. Add MadVoro to your source tree** (e.g. as a git submodule or by copying):
```bash
cd your_project
git submodule add https://github.com/maormizrachi/MadVoro.git external/MadVoro
```

**2. In your `CMakeLists.txt`:**
```cmake
# Enable the features you need before adding MadVoro
set(MADVORO_WITH_MPI ON CACHE BOOL "")

add_subdirectory(external/MadVoro)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE MadVoro::madvoro)
```

That's it. All include paths, compile definitions (e.g. `MADVORO_WITH_MPI`), C++17 requirement, and link dependencies (Boost, MPI, etc.) are automatically propagated to your target through `MadVoro::madvoro`.

**3. Use MadVoro with your own vector type:**

MadVoro is templated on a user-provided 3D vector type. Your type must have public `double x, y, z` members and be constructible from three doubles via brace initialization (`Vec3{x, y, z}`).

```cpp
#include <madvoro/Voronoi3D.hpp>

// Your own 3D vector type
struct Vec3 {
    double x, y, z;
};

int main() {
    Vec3 ll{0, 0, 0}, ur{1, 1, 1};
    MadVoro::Voronoi3D<Vec3> voronoi(ll, ur);

    std::vector<Vec3> points = { {0.1, 0.2, 0.3}, {0.8, 0.5, 0.7} };
    voronoi.Build(points);

    Vec3 cm = voronoi.GetCellCM(0);
    // ... use cm.x, cm.y, cm.z
}
```

`MadVoro::Face<Vec3>` is also templated, so face vertices are returned in your type.

### CMake Options

Set these **before** `add_subdirectory()` (or pass as `-D` flags):

| Option | Default | Description |
|---|---|---|
| `MADVORO_WITH_MPI`  | `OFF` | Build with MPI support |
| `MADVORO_WITH_HDF5` | `OFF` | Build with HDF5 output support |
| `MADVORO_WITH_VTK`  | `OFF` | Build with VTK visualization support |
| `MADVORO_WITH_VCL`  | `OFF` | Build with VCL vectorization acceleration |
| `MADVORO_BUILD_EXAMPLES` | `OFF` | Build example programs |
| `BUILD_SHARED_LIBS`  | `OFF` | Build shared library (`.so`) instead of static (`.a`) |

### Dependency Hints

If CMake cannot find a dependency automatically, pass these hints:

| Dependency | Hint |
|---|---|
| Boost | `-DBOOST_ROOT=/path/to/boost` |
| MPI | `-DMPI_CXX_COMPILER=/path/to/mpicxx` |
| HDF5 | `-DHDF5_ROOT=/path/to/hdf5` |
| VTK | `-DVTK_DIR=/path/to/vtk/lib/cmake/vtk` |
| VCL | `-DVCL_DIR=/path/to/vcl` |

## Standalone Build & Install

If you prefer to build and install MadVoro as a standalone library (e.g. system-wide), you can do so and then use `find_package` in your project.

### Build
```bash
git clone https://github.com/maormizrachi/MadVoro.git
cd MadVoro
cmake -B build -DMADVORO_WITH_MPI=ON -DCMAKE_INSTALL_PREFIX=$(pwd)/install
cmake --build build -j$(nproc)
cmake --install build
```

### Use via `find_package`
```cmake
find_package(MadVoro REQUIRED)
target_link_libraries(your_target PRIVATE MadVoro::madvoro)
```

You may need to pass `-DMadVoro_DIR=/path/to/install/lib/cmake/MadVoro` if the install location is non-standard.

## Usage
### API
MadVoro offers a wide API by merely giving the list of points to build and the construction zone (usually a box used to clip the Voronoi cells), including a cell's vertices, faces, the list of a cell's neighbors, cell's center of mass and faces center of mass.
The user API is based on two class templates:
- `Voronoi3D<Vec3>`, representing a distributed three-dimensional Voronoi diagram.
- `Face<Vec3>`, representing a face defined by vertices of your vector type.

These templates accept any 3D vector type with public `double x, y, z` members. MadVoro does **not** define its own vector type — you use your project's existing one.

> [!IMPORTANT]
> You can find the full API of each class in the include files (`include/madvoro/`). Pay attention to methods documentations.

### Vec3 Requirements

Your vector type must satisfy:
1. Public `double x, y, z` members (for reading coordinates)
2. Brace-constructible from three doubles: `Vec3{x, y, z}` (for returning results)

No operators, no math functions, no inheritance required.

### Examples
To build and run examples:
```bash
cmake -B build -DMADVORO_BUILD_EXAMPLES=ON -DMADVORO_WITH_MPI=ON
cmake --build build -j$(nproc)
```

The examples include a simple `Vector3D` type (in `examples/Vector3D.hpp`) as a reference implementation.

#### Serial Examples
```bash
./build/examples/example_faces_information
./build/examples/example_uniform_serial
```
> [!WARNING]
> Serial examples should not be run when your project is compiled with MPI support.

#### Parallel Examples
```bash
mpirun -n 16 ./build/examples/example_pentagon
mpirun -n 8  ./build/examples/example_uniform_parallel
```

## Cleaning
To start fresh, simply remove the build directory:
```bash
rm -rf build
```

## Support and Contact
If you run into problems or difficulties in compiling or running, or have any questions or suggestions, feel free to contact me by email: maor.mizrachi@mail.huji.ac.il.

## Reference
If you wish to cite our work, we would appreciate it if you used the following BibTeX citation:
```
@article{10.1093/rasti/rzaf039,
    author  = {Mizrachi, Maor and Raveh, Barak and Steinberg, Elad},
    title   = {madvoro: parallel construction of Voronoi diagrams in distributed memory systems},
    journal = {RAS Techniques and Instruments},
    volume  = {4},
    pages   = {rzaf039},
    year    = {2025},
    month   = {09},
    issn    = {2752-8200},
    doi     = {10.1093/rasti/rzaf039},
    url     = {https://doi.org/10.1093/rasti/rzaf039},
    eprint  = {https://academic.oup.com/rasti/article-pdf/doi/10.1093/rasti/rzaf039/64231338/rzaf039.pdf},
    abstract = {Voronoi diagrams are essential geometrical structures with numerous applications, particularly astrophysics-driven finite volume methods. While serial algorithms for constructing these entities are well-established, parallel construction remains challenging. This is especially true in distributed memory systems, where each host manages only a subset of the input points. This process requires redistributing points across hosts and accurately computing the corresponding Voronoi cells. In this paper, we introduce a new distributed construction algorithm, which is implemented in our open-source C++ 3D Voronoi construction framework. Our approach leverages Delaunay triangulation as an intermediate step, which is then transformed into a Voronoi diagram. We introduce the algorithms we implemented for the precise construction and our load-balancing approach and compare the running time with other state-of-the-art frameworks. madvoro is a versatile tool that can be applied in various scientific domains, such as mesh decomposition, computational physics, chemistry, and machine learning.}
}
```
