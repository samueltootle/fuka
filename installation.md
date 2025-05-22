\page install Installation

To utilize the FUKA initial data codes and utilities, it is first
necessary to compile the FUKA library which is composed both of the components
related to the KADATH spectral library as well as FUKA specific solvers and supporting
utilities along with the initial data exporters which are utilized when importing
FUKA initial data into an evolution framework.

# 1. Dependencies

FUKA initial data solvers require the following dependencies:
1. Compiler for C/C++ code (tested with Intel and GCC)
2. CMake > 2.8
3. Boost libraries including `ptree`
4. MPI (e.g. OpenMPI, Intel-MPI, MPT)
5. Git
6. ScaLAPACK - Parallel LAPACK package (MKL also works)
7. FFTW3
8. GSL

For compiling the FUKA PythonTools, the following additional dependencies are needed
1. Python3
2. Boost libraries including `Python`

#### For compiling with `GRHayLEOS` support:
FUKA currently uses an extended version of the `GRHayL` library (pull request pending). To use GRHayL
with FUKA, build using `additional_eos_tools` branch from `S. Tootle's` fork of [GRHayL](https://github.com/samueltootle/GRHayL/)

After installing `GRHayL`, set the environment variable `GRHAYL_ROOT` to
the installation directory in order for FUKA's `cmake` build scripts to
find it.  Finally, `-DGRHAYL_EOS=ON` needs to be added to all cmake
commands to ensure `GRHayL` support is active at build time.

# 2. Repository sync

One can contain the latest version of FUKA using git using the following:

`git clone https://bitbucket.org/fukaws/fuka`

<b>
Note:

1. The latest build is available under the `fuka` branch as opposed to
   master or main.

2. Stable releases are available as specific versions, e.g. `git checkout
   fukav2,2`.

</b>

# 3. Environment Setup

In order to build the FUKA library, it is mandatory that one sets the
`HOME_KADATH` environment variable.  If one is running linux, is at a
bash command line interface, and is currently in the root directory of
FUKA (e.g. `$HOME/fuka`), one can run

```bash
export HOME_KADATH=`pwd`
```

Furthermore, if one wants to utilize the `compile` build script (see
below), it is necessary to set these additional environment variables:

1. `KAD_CC` - e.g. `export KAD_CC=gcc`
2. `KAD_CXX` - e.g. `export KAD_CC=g++`
3. `KAD_NUMC` - e.g. `export KAD_NUMC=7` (number of parallel build tasks)
4. `GRHAYL_ROOT` - (optional) the installation path for the `GRHayL` library

Finally, one can ensure these are loaded whenever starting a new terminal
session by setting the same commands in the user's RC file (e.g.
`$HOME/.bashrc`, for bash).

# 4. Cmake Setup

To ensure `cmake` is able to find all of the necessary libraries, it is
important to set these dependencies manually especially on high
performance computing clusters where packages are installed in
non-default locations.  To do so, one needs to modify the file
`Cmake/CMakeLocal.cmake`.  The default looks like the follow:

```
set (PGPLOT_LIBRARIES "/usr/lib/libpgplot.so.5")
set (GSL_LIBRARIES "/usr/lib/x86_64-linux-gnu/libgsl.so")
set (SCALAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so")
set (FFTW_LIBRARIES "/usr/lib/x86_64-linux-gnu/libfftw3.so")
set (BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so")
set (LAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/liblapack.so")
```

# 5. Compilation - Script

To compile the FUKA library using the automated script, one can go to
`$HOME_KADATH/build_release` and run the shell script using

`. compile`

<b>
This script may need to be modified to include `-DGRHAYL_EOS=ON` if
GRHayL support is required.
</b>

So long as the previous steps have been done properly, compilation should
be successful and the static library will be located in
`$HOME_KADATH/lib/libkadath.a`

# 6. Compilation - Manually

- From `$HOME_KADATH/build_release` or `$HOME_KADATH/build_debug`, create
a build directory where all of the temporary `cmake` and `make` files can
be stored.
- Invoke `cmake (options) ..` where options consist of:
    - `-DPAR_VERSION = On/Off (On)` Set to On to build the MPI parallel
        version of the library. The initial data codes within this branch
        are only designed for use with MPI.
    - `-DCMAKE_BUILD_TYPE = Release/Debug`
        Specifies the build type
    - `-DMPI_CXX_COMPILER=mpicxx`
        Path to the MPI C++ wrapper (when not automatically detected by cmake)
    - `-DMPI_C_COMPILER=mpicc`
        Path to the MPI C wrapper (when not automatically detected by cmake)
    - `-DGRHAYL_EOS=OFF`
        Set to `ON` to enable `GRHayLEOS` support
-  Once cmake has been successfully invoked, use `make -j $KAD_NUMC` or specify
the number of parallel tasks manually to start the compilation.

# 7. Compilation - Initial Data

The latest solvers can be found in the
`$HOME_KADATH/codes/FUKAv2_Solvers` base directory. As discussed in the
[Organization](index.html) section, all of the initial data solvers have
a similar structure.

To compile an initial data solver once the FUKA library has been
successfully compiled, enter the base directory for the solver, e.g.
`$HOME_KADATH/codes/FUKAv2_Solvers/BH` and execute the `compile` script.

Alternatively, one can follow the instructions for manual compilation as
discussed previously.

After successful compilation, the binary files can be found in

`<basedir>/bin/Release`

# 8. Compilation - PythonTools

Within `$HOME_KADATH/codes/PythonTools` one can build the Python
libraries that allow for analyzing FUKA initial data solutions within
Python.  Furthermore, a set of Python utilities are available to quickly
enable plotting 1D and 2D plots.

So long as the dependencies are met, compilation is can be done by using
the `compile` script or by manually compilation as discussed previously.

<b>

Note: The handling of Python and Boost libraries utilizes cmake's
built-in find capabilities to find Python and Boost libraries along with
the compiled Boost modules. If you have Boost libraries installed
correctly, but cmake does not find them, you can set the environment
variable `BOOST_ROOT` to the base directory where Boost has been
installed to.

</b>