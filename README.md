# MadVoro - Massively distributed Construction of Voronoi Diagrams
![C++ project](https://img.shields.io/badge/C++-3393ff)

MadVoro is a C++ framework for construction of Voronoi diagrams of 3D points, in distributed memory (using MPI).

<img src="examples/fox/fox.png?raw=true" alt="An example for a fox mesh, in a 16 processors construction." width="500"/><img src="examples/pyramid/pyramid.png?raw=true" alt="An example for a mesh construction of a pyramidal space." width="300"/>


## Requirements
- Any C++ Compiler (tested with g++ and icpx) 
- Boost >= 1.74.0

### Optional
- For parallel support: Any MPI implementation (tested with openmpi, IntelMPI and MPICH). 
- Recommended: [VCL](https://github.com/vectorclass/version2) for vectorization accelartion.
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) (C++ binding) >= 1.8.0, for output.
- [VTK](https://vtk.org/) >= 1.9.3, for visualization.

## Build & Install
First, clone this git repository and change the directory to the cloned
```
git clone https://github.com/maormizrachi/MadVoro.git
cd MadVoro
```
Then, run the configuration file, and compile:
```
./configure --with-boost(=DIR)
make
make install
```
Follow the instructions printed in `make install` to install the library successfuly (in particular, it is possible that the path to library should be added to LD_LIBRARY_PATH).

By default, running `make install` will install the binary files in /usr/local/bin, the libraries in /usr/local/lib and so on. That's probably not what you want if you don't have permissions to write in those directories (and the install phase will fail).
You can configure another prefix, where the project will be installed in, by using the `--prefix(=DIR)` flag:
```
./configure <...> --prefix(=DIR)
```
For additional info, and for more flags and settings, you can always use `./configure --help`.

One may not supply `--with-boost(=DIR)`, or supply `--with-boost` only, if the boost directory is defined in an environment variable called `BOOST_DIR` or `BOOST_ROOT`, assuming one of those contains boost's `include` (and `lib`, if exists) directories. If not, you must supply the lib directory and include directory apart:
```
./configure --with-boost-libdir(=LIBDIR) --with-boost-include(=INCDIR)
```
Assuming `LIBDIR` contains boost's libraries (shared object files), and `INCDIR` contains boost's header files.

### With MPI
To compile with MPI, two options are possible when running `./configure`:
```
./configure <...> MPICXX=/path/to/mpicxx
```
Where `mpicxx` is the MPI C++ compiler.
Or, alternatively:
```
./configure <...> --with-mpi(=DIR)
```
Where `DIR` is a directory containing `bin`, `include` and `lib` directory of the desired MPI package.
Supplying `--with-mpi` with no argument, searches for the MPI implementation of the MPI that implements the command `mpicxx`.

### Build with output
MadVoro supports an output system, allowing you to extract the distributed Voronoi diagram into a 3-dimensional .pvtu file (parallel VTK). These files can be read by programs such as paraview.
Use `--with-vtk` to build with vtk support. You may supply your VTK path after the `=`.
You can also use `--with-hdf5` to use HDF5.
When adding VTK a new method will be added to the `Voronoi3D` class, called `ToVTK`. The latter method gets a filename to print the VTK file to. It may also get additional data to supply with each point (You may add double/float64 fields by specifying the field name and field values for each point, by giving each in a dedicated vector parameter to this method).  
A similar method called `ToHDF5` will be added in case you compiled the project with HDF5.

## Running
### API
MadVoro offers a wide API by merely giving the list of points to build and the construction zone (usually a box used to clip the Voronoi cells), including a cell's vertices, faces, the list of a cell's neighbors, cell's center of mass and faces center of mass.
The API is available to the user by using three classes:
- `Voronoi3D`, representing a distributed three-dimensional Voronoi diagram.
- `Face`, representing an area limited by points (composed by indices of local points as points).
- `Vector3D`, representing a three-dimensional point containing x, y, z coordinates and basic operations.
> [!IMPORTANT]
> You can find the full API of each class in the include files generated after you install the library. Pay attantion to methods documentations.

You should write your own C++ project, and include the three headers files (found in your include) allowing you to use the code. You can see examples in the examples directory.
After that, you should add an include path link path to your compilation. Specifically, you should add the following flags in compilation:
```
-I /path/to/madvoro/include -L /path/to/madvoro/lib -lmadvoro
```
Where /path/to/madvoro is the prefix you have set in configuration (or alternatively, the default of the configure script).
This is enough for compiling, yet not enough for running. As MadVoro is linked by a dynamic (shared object) library, you should add /path/to/madvoro/lib to your LD_LIBRARY_PATH environment variable. The way of doing that differs by different shells. In bash, you can run:
```
export LD_LIBRARY_PATH=/path/to/madvoro/lib:$LD_LIBRARY_PATH
``` 
It is recommended to add this line to your `~/.bashrc` file, in this case. 
> [!TIP]  
> You may use other techniques to run, avoiding the use of `LD_LIBRARY_PATH`. For example, using RUNPATH or different link-flags. Instructions are shown when running `make install`. 
### Examples
You can find examples in the `examples` directory.

## Cleaning
You can clean the project by writing the following make command:
```
make distclean
```
> [!NOTE]
> To recompile with different configuration (for example, add a support of HDF5 if you didn't have one currently), it is required to clean and reconfigure the project again.

## Reference
If you wish to cite us, we'll be thankful if you'd used the following BibTeX format:
```
TODO
```