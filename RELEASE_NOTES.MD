# Frankfurt University/Kadath Initial Data (FUKA) - Release Notes

## v2.3
### Removed Features
* FUKAv1 Codes have been removed.
    * FUKA solvers are now stored in `codes/FUKA`
### New Features
* More robust algorithm for determining grid structure for BHNS and BBH systems
* Added support for GRHaYLEOS
    * Needs to be installed separately and the base root set as GRHAYL_ROOT environment variable
    * FUKA needs to be compiled with GRHAYL_EOS macro defined, e.g.
    `cmake -DPAR_VERSION=ON -DCMAKE_BUILD_TYPE=Release -DGRHAYL_EOS=ON  ..`
    * GRHaYLEOS support enables native support for stellar collapse 3D tables
* (beta) Multi-threaded support for **export only**
    * `include/memory.hpp`: Need to define `#define DEFAULT_KAD_MEM`
### Solvers in development
* (beta) Support for differentially rotating neutron star solutions in XCTS coordinates
* (beta) Support for differentially rotating neutron star solutions in Quasi-isotropic coordinates
* (alpha) Head-on initial data for BNS and BHNS
    * Black holes can be boosted

## v2.2
* With this release, user demanded shells around compact objects have been removed.
In its place is a significantly more robust solution
which determines an optimal configuration of spherical grids around each compact object in binary spaces.
* A minimal Config file has been introduced and is now the default for all solvers. In this way the user is only
concerned with the parameters necessary to generate initial data.
    * The "full" Config file must still be used for the ID solutions
    * The "full" Config file can be obtained by running `solve full` if desired for a new sequence
* (Beta) A 1D sequences are now supported such that a series of initial data solutions can be automated

## v2.1
This release includes a rewrite and reorganization of the Python readers for FUKA initial data among many bug fixes and quality of life changes.
The old python readers have thus been removed from the solvers and are now centrally located in the [codes/PythonTools](https://bitbucket.org/fukaws/fuka/src/fuka/codes/PythonTools) directory.

Along side of the FUKAv2.1 release is the official inclusion of the thorns to import FUKA initial data into the Einstein Toolkit - [ticket](https://bitbucket.org/einsteintoolkit/tickets/issues/2692/inclusion-of-fuka-importer-thorns). Please note the import thorns are not activated by default due to the need for Boost libraries.  For an example thorn to provide Boost libraries within the Einstein Toolkit, see [https://github.com/dradice/Boost](https://github.com/dradice/Boost).

## v2.0
The release of FUKAv2 is considerable step forward in reliable generation of extremal spinning, asymmetric
binary initial data using the KADATH spectral software.  The v2 solvers aim to maximize convergence by using super-imposed
isolated solutions to setup the initial guess for binary ID.  Additionally, v2 aims to automate the generation of ID by minimizing
the workflow for the user.  Finally, v2 allows quite for considerable flexibility in setting up the config file to make ID generation
as efficient as possible within the KADATH framework including reusing previously solved implicit isolated solutions.  For more details,
please see the documentation in the [FUKAv2](https://bitbucket.org/fukaws/fuka/src/fuka/codes/FUKAv2_Solvers/).

  - FUKAv1 Specific: There have been some core changes to various utilities used when constructing ID that have required some refactoring of the v1 solvers.  Spot testing has been done to ensure these codes function as originally intended - however - these codes, overall, remain unchanged.

## v1.0
The first release of the FUKA codes were released in 2021 alongside the code paper which can be found here https://arxiv.org/abs/2103.09911

The original codes can be found in [FUKAv1](https://bitbucket.org/fukaws/fuka/src/fuka/codes/FUKAv1_Solvers/) and exist primarily
for historical purposes and should not be used by new users.  However, those that wish to develop more complicated codes within
KADATH may find them useful as their structure is less abstract.

**NOTE: After the v2.0 release, these codes are no longer modified or maintained.**