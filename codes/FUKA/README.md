\mainpage FUKA Initial Data
# FUKAv2

Author(s)    : Samuel D. Tootle,  L. Jens Papenfort  
Maintainer(s): Samuel D. Tootle - tootle@itp.uni-frankfurt.de  
License      : GPLv3+ for all other code  

---

# 1. Overview

The FUKAv2 encompasses a rewrite of the original initial data (ID) solvers that were made publicly
available in 2021.  Although portions of the original codes have been reused, a significant
overhaul of the solvers has been done not only to make the interface more accessible to the user, 
but to also provide a coherent interface between the isolated and binary solvers.

Aside from the abstraction and interface, the key component that makes the v2 binary solvers
considerably more reliable and functional is the use of super-imposed boosted isolated objects
to construct the "initial guess" for the binary solver.  This allows *nearly* arbitrary configurations
to be constructed from scratch without iterative changes to the binary objects as was required in
the FUKAv1 binary solvers.

The *nearly* caveat is a function not only of resolution, but also the domain decomposition 
near the excision or stellar surfaces due to Gibbs phenomena that occur.  In the vast majority of test
cases this does not prohibit the construction of ID even at extreme mass ratios `q^{-1} > 20`.
It remains an on-going effort to optimize the domain decomposition to minimize these effects.

# 2. Nomenclature and Units

- In all FUKA solvers we work in geometric units such that `G = c = 1` unless otherwise specifically stated [^1]
- There will be mentions of *implicit solutions* throughout the various solver READMEs.  Each subsequent stage for
a given solver builds on the solution of previous stages.  These solutions usually have some useful characteristics 
(e.g. a non-rotating TOV solution prior to generating a rotating one; all stages of the isolated solutions when generating binary ID) 
which can be reused or analyzed later. However, these solutions do not match what has been requested by the user, 
they are simply a consequence of the ID solving process.
- In the case of binary ID, we use the convention throughout the documentation and within the FUKAv2 solvers that `q <= 1`
- When fixing the spin, it is common to use the dimensionless spin parameter `chi = S / M^2` where 
    - S is the spin angular momentum as measured (quasi)-locally on the compact object
    - `M` is the compact object's gravitational mass.  
        - In the case of NS ID, `M` is the ADM mass as computed from the isolated solution for a given rest mass and spin
        - In the case of BH ID, `M` is the Christodoulou mass
- For details on the formulations used within these solvers, please review the original [publication](https://arxiv.org/abs/2103.09911) 

[^1]: There exists a variable in the config file call `qpig` which is equivalent to `4*PI*G`.  This was originally designed to allow an arbitrary rescaling of the system of equations, matter, etc; however, this has not been implemented consistently as there was no clear benefit other than confusion when working in different unit systems

# 3. Organization

The solvers are structured as follows:

1. Each directory contains the files needed for the respective solver
    1. `CMakeLists.txt` is the file needed by CMake to compile the codes
    2. `compile` is a symbolic link to the script stored in `$HOME_KADATH/build_release` to ease compiling
    3. `src` directory contains the relevant source files:
        - `solve.cpp`: the one and done solve code for each solver
        - `reader.cpp`: the reader can provide diagnostics from ID solutions that are computed from the ID
2. The related `solve.cpp` in most cases is simply a frontend with minimal functionality.  
Most of the business is stored in `$HOME_KADATH/include/Solvers`
3. The `reader.cpp` code is currently shared between v1 and v2 solvers, therefore symbolic links are used

# 4. Getting Started

First and foremost, see \ref install for instructions on how to obtain and build
the FUKA library, initial data codes, and related PythonTools.

Once the intended features have been successfully compiled, the \ref tutorials are
a great place to get acquainted with FUKA!

# 5. Acknowledgements

I have been developing FUKAv2 since before the initial public release in an 
effort to automate ID construction across the space of parameters.
To that end, I am very grateful to my original collaborators 
L. Jens Papenfort and Elias R. Most for their helpful discussions.  

Furthermore, I would like to thank Konrad Topolski for being an early adopter of the v2 codes and for his continued feedback and support throughout the many changes since
then.

# 6. Outstanding Tasks

The following are on the list of things to do based on expected level of effort:

1. Rewrite exporters to make them more readable and consistent
2. Add documentation for `PythonTools` to make them more accessible
3. Implement solvers for BH, BBH, and BHNS that allow for maximal spinning black holes

# 7. Contributing

1. Feedback and bugs are always a welcomed contribution.  One can always e-mail the maintainers for issues experienced or configurations that provide challenges as this can help to diagnose limitations in the solvers
2. Help with the current outstanding tasks is of course possible whether through code or experience with certain issues, e.g. optimizing a fix grid setup to avoid spectral artifacts.  For code contributions, forking and submitting pull requests is possible.