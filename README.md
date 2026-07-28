\mainpage FUKA Reference
# Frankfurt University/Kadath Initial Data branch
#### Author(s)    : L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément

## Overview - Updated
Included are the Frankfurt initial data solvers and utilities based on the Kadath
  spectral solver library.  The original solvers written by the aforementioned authors
  (hereafter denoted as FUKAv1) are deprecated, but can be obtained by checking out
  the previous version of FUKA `git checkout fukav2.1` which are then located in
  `./codes/FUKAv1/[BH, NS, BHNS, BNS, BBH]` respectively.

  The FUKA solvers from now are can be found in ./codes/FUKA/[BH, NS, BHNS, BNS, BBH] respectively.
	The latest FUKA solvers includes support for polytropic equations of state as well as tabulated EOS
  in the standard LORENE format.  Examples and additional details can be found in the [eos](https://github.com/samueltootle/fuka/tree/fuka/eos/) directory.

  As of FUKAv2.3, FUKA now supports *stellar collapse* format equations of state using the `GRHayLEOS` library.
  The `GRHayL` can be found [here](https://github.com/GRHayL/GRHayL) and needs to be installed separately before
  compiling the Kadath library and, after, the FUKA codes.

## FUKAv2.4 Notes:
1. The solution INFO file for new solutions now contains a `metadata`
   section which includes
    * build date
    * FUKA version
    * GIT information
2. `metadata` will take the place of the temporary fix, `new_initial_data` control.
    * `fuka_version` is now used to distinguish initial data generations.
      The other parameters are just for information and debugging.
    * `new_initial_data` can still be used to distinguish solutions from \< v2.3.  This currently only effects BNS solutions before and after v2.3, as noted below in FUKAv2.3 notes.
    * Alternatively, one can manually add metadata to their previous solutions, at a minimum:

```
metadata
{
    fuka_version v2.2.0
}
```

## FUKAv2.3 Notes:
### Refactored features
1. Many components of the Kadath library have been refactored to support copy constructor operations
    * Note: This only works with memory pools turned off, e.g. `-DDEFAULT_KAD_MEM` and should only be used
    for exporting ID solutions to an evolution framework
2. The EOS module has been overhauled to enable easier incorporation of EOS backends, see e.g. `include/EOS/FUKA_EOS_Wrapper.hh`
3. The algorithm for computing shells around objects has been revised again enabling more accurate and efficient calculations of
highly asymmetric binaries at very close and very large separations.
4. FUKA Python readers have been refactored to use the new EOS backends.
    * New readers have been added to support new ID solutions
5. The BNS grid has been refactored to allow for external spherical shells. **Support for this new feature breaks support for old ID solutions to a degree**
    * **Old ID solutions are assumed by default**
    * **New solutions are distinguished based on the control `new_initial_data on` in the `sequence_controls` section of the INFO file.**

### New Features
* A new solver that computes solutions of isolated neutron stars in quasi-isotropic coordinates is now available, see [NS_isotropic](./codes/FUKA/NS_isotropic/)
  * This solver supports irrotational, uniform, and differential rotation models.
  * Currently differential rotation models are limited to the
  [KEH law](https://ui.adsabs.harvard.edu/link_gateway/1989MNRAS.239..153K/doi:10.1093/mnras/239.1.153)
* A new solver that computes solutions of differentially rotating neutron stars in XCTS coordinates is now available,
see [NS_DIFFROT](./codes/FUKA/NS_DIFFROT/)
* Support for temperature dependent EOS' in *stellar collapse* format is now available
  * Users must first build and install the [GRHayL](https://github.com/GRHayL/GRHayL) library before building the Kadath library and, after, the FUKA solvers.
  * At compile time for the Kadath library and FUKA solvers, `cmake` searches the environment variables `GRHAYL_ROOT` for `GRHayL` resources, so this must be set manually by users.
  * Also, `cmake` must be ran with `-DGRHAYL_EOS=ON`
* A new suite of exporters are now available for all ID solutions.  These exporters leverage the new copy constructors in Kadath
to allow for multi-threaded import by evolution frameworks.
    * Note: This only works with memory pools turned off, e.g. compiling Kadath and executables with `-DDEFAULT_KAD_MEM`
* A new code for computing head-on initial data for two Neutron stars at rest is now available [BNS_headon](./codes/FUKA/BNS_headon/)

## FUKA Maintainer(s):

Samuel D. Tootle - tootle@itp.uni-frankfurt.de

## KADATH Maintainer:
Philippe Grandclément - philippe.grandclement@obspm.fr

License      : GPLv3+ for all other code

## REQUIRED CITATIONS:

1) L. Jens Papenfort, Samuel D. Tootle, Philippe Grandclément, Elias R. Most, Luciano Rezzolla: https://arxiv.org/abs/2103.09911

2) Philippe Grandclément, http://dx.doi.org/10.1016/j.jcp.2010.01.005

# 1. Purpose

This collection of ID solvers aims at delivering consistent initial data (ID)
solutions to the eXtended Conformal Thin-Sandwich (XCTS) formulation of
Einstein's field equations for a variety of compact object configurations.

As each solver has their own specific nuances and considerations, we have included
a README in each solver directory to provide a basis for getting started with the
respective solver.

Additionally, each initial data has a respective exporter which can be seen
in `src/Utilities/Exporters`.  These exporters allow one to compile an interface code
for an evolution framework along with the Kadath static library located in `$HOME_KADATH/lib`, in
order to export data based on input grid points.

# 2. Modifications from base Kadath

In addition to the solving routines included within `$HOME_KADATH/codes`, we also note the major overall modifications
and additions that differ from base Kadath.
1.  This branch includes memory optimizations that inspired portions of the optimization (now main) branch
2.  Modification/addition of numerical spaces for the BH, BBH, BNS, and BHNS
3.  Addition of an equation of state infrastructure utilizing Margherita standalone to handle
tabulated and polytropic EOS - see [include/EOS](https://github.com/samueltootle/fuka/tree/fuka/include/EOS)
4.  Addition of the Configurator framework to enable extensibility of solvers by managing controls,
stages, and key variables - see [include/Configurator](https://github.com/samueltootle/fuka/tree/fuka/include/Configurator)
5.  Addition of exporters for all the previously mentioned ID types - see [src/Utilities/Exporters](https://github.com/samueltootle/fuka/tree/fuka/src/Utilities/Exporters)

**Note: as of summer 2021, the FUKA solvers are based on the deprecated branch of Kadath.  Given the optimizations and changes made
within the FUKA branch conflict with those implimented in the `master` branch (previously the `optimized` branch), a considerable
level of effort is required to merge FUKA with the new `master` branch as well as test to see which optimizations provide better results.
Currently, there is no timeline for when this will be done.**

# 3. Public Thorns for use with the Einstein Toolkit

The following workspace includes the FUKA ID respository (including versioned branches)
as well as available thorns for use with the Einstein Toolkit in order to import FUKA ID:
https://bitbucket.org/fukaws/



# 4. Compiling KADATH

## Get the Sources

`git clone git@bitbucket.org:fukaws/fuka.git`

## Set the following environment variables based on your setup in your ~/.bashrc file

- HOME_KADATH - e.g., export HOME_KADATH=/home/user/lib/fuka/
- KAD_CC - e.g. gcc
- KAD_CXX - e.g. g++
- KAD_NUMC - number of parallel compiling jobs CMake can run, e.g. 7

## Build Process

1. Go to build_release.
2. Create a build directory and enter it
3. Invoke `cmake (options) ..`
  - where `..` denotes the location where the CMakeList.txt file is
  - The important cmake options are the following (the value in parentheses corresponds to the default settings) :
    - `-DPAR_VERSION = On/Off (On)`
      - Set to On to build the MPI parallel version of the library. The initial data codes within this branch are only designed for use with MPI.
    - `-DCMAKE_BUILD_TYPE = Release/Debug`
      - Specifies the build type which results in different compiler options being used
    - `-DMPI_CXX_COMPILER = <compiler pile>`
      - Path to the MPI C++ wrapper (when not automatically detected by cmake)
    - `-DMPI_C_COMPILER  = <compiler pile>`
      - Path to the MPI C wrapper (when not automatically detected by cmake)
    - (optional) `-DDEFAULT_KAD_MEM`
      - Deactivate memory pools to allow for thread-safe exporters.  ***Do not use with ID solvers!***
    - (optional) `-DGRHAYL_EOS=ON`
      - Enable support for stellar collapse tables using the GRHayLEOS library.
      - GRHayL must be installed separately by the user
      - The environment variables `GRHAYL_ROOT` must be set by the user to the directory GRHayL was installed to.

Example using GNU+mpi compilers:
    `cmake -DCMAKE_BUILD_TYPE=Release -DPAR_VERSION=On -DMPI_CXX_COMPILER=mpic++ -DMPI_C_COMPILER=mpicc -DGRHAYL_EOS=OFF ..`

In most HPC systems, `cmake` will likely not find the dependency libraries that the user may intend.  Therefore,
one must specify them manually through the [CMakeLocal.cmake](https://github.com/samueltootle/fuka/tree/fuka/Cmake/CMakeLocal.cmake) file
(the `fftw` and `scalapack` libraries must usually be provided in this way).
Some working [CMakeLocal.cmake](https://github.com/samueltootle/fuka/tree/fuka/Cmake/CMakeLocal.cmake) files
are provided for HPC systems in Germany as well as examples for personal computers.

Once cmake has been successfully invoked, use make -j $KAD_NUMC to start the compilation.

## Compiling the library with the compile script

A script called [compile](https://github.com/samueltootle/fuka/tree/fuka/build_release/compile) is also provided that can be used to facilitate the installation process. So long as the above environment variables are set and the libraries are found, no additional input is necessary.
Run the compile script within the `build_release` directory using

`. compile`

in order to build the library.

## Compiling FUKA solvers
The above mentioned [compile script](https://github.com/samueltootle/fuka/tree/fuka/build_release/compile) has been added as a symbolic link to the FUKAv1 and FUKAv2 solver directories for convenience to compile the individual solvers.

# 5. Dependencies

1. C++ compiler - gcc is recommended
    - Must support c++17 standards
    - Must support `<filesystem>`
2. Cmake
3. git
4. FFTW3
5. GSL
6. scaLAPACK
7. MPI
8. Boost
9. Boost::python (to compile [PythonTools](https://github.com/samueltootle/fuka/tree/fuka/codes/PythonTools/))
10. (optional) [GRHayL](https://github.com/GRHayL/GRHayL)