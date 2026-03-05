# Pentagon (Prism) Example

Constructs a Voronoi diagram inside an 80-sided prism (approximating a cylinder) using custom boundary faces. Optionally exports the result to VTK.

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

If built with VTK support, output is written to `pentagon_example.vtu`. Open it in [ParaView](https://www.paraview.org/) to visualize.
