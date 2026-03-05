# Parallel Output VTK Example

Benchmarks repeated parallel Voronoi builds. Runs 15 iterations, reporting the time for the first build and the average time for subsequent (incremental) builds.

Each MPI rank reads its own input file from `input/<rank>`.

Requires MadVoro built with **MPI** support. Does not require VTK or HDF5 (despite the directory name, this example does not produce VTK output).

## Prerequisites

Create an `input/` directory with one file per MPI rank (named `0`, `1`, ..., `N-1`). Each file contains points in the format `(x,y,z)`, one per line.

## Building

```bash
MADVORO_PREFIX=... BOOST_DIR=... bash build.sh
```

## Running

```bash
mpirun -np <N> ./test
```
