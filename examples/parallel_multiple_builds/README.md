# Parallel Multiple Builds Example

Demonstrates building a parallel Voronoi diagram and optionally exporting the result to VTK. Each MPI rank reads its own input file from `input/<rank>`.

## Prerequisites

Create an `input/` directory with one file per MPI rank (named `0`, `1`, ..., `N-1`). Each file contains points in the format `(x,y,z)`, one per line.

## Building

```bash
MADVORO_PREFIX=... BOOST_DIR=... bash build.sh
```

To enable VTK output:

```bash
MADVORO_PREFIX=... BOOST_DIR=... VTK_DIR=... HDF5_DIR=... bash build.sh
```

## Running

```bash
mpirun -np <N> ./test
```

If built with VTK support, output is written to `output.vtk`.
