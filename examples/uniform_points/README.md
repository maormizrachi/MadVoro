# Uniform Points Example

Benchmarks Voronoi construction on randomly generated uniform points. Includes both a serial and a parallel (MPI) version.

Does not require VTK or HDF5.

## Building

```bash
MADVORO_PREFIX=... BOOST_DIR=... bash build.sh
```

This produces two binaries: `test_serial` and `test_parallel`.

## Running

Serial:

```bash
./test_serial [N]
```

Parallel:

```bash
mpirun -np <P> ./test_parallel [N]
```

`N` is the number of points (default: 10000). In the parallel version, each rank generates `N` points.
