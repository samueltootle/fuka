\page tutorials Tutorials

# 1. Introduction

Welcome to the tutorials for the FUKA initial data solvers. The following
offers a brief introduction to all of the current initial data codes and is organized
to provide an iterative and intuitive introduction to their uses.

The user interface for all the solvers is in the form of a Configurator (config) file.  
See the \subpage configreadme for details.
To this end, using all the ID solvers is the same:

1. Generate the initial config file by running `solve` for the first time
2. Modify the initial config file based on ID characteristics you are interested in
3. Rerun (using parallelization) using this config file, e.g. `mpirun solve initial_bh.info`

For this reason, it is recommended to learn about using FUKA solvers by generating ID in the following order

-  \ref bhxcts
-  \ref bbhxcts
-  \ref nsxcts
-  \ref bnsxcts
-  \ref bhnsxcts

Furthermore, details related to the use of polytropic and tabulated equations of state can be found in \ref eos

# 2. Maintainers

Author/Maintainer: Samuel. D. Tootle - tootle@itp.uni-frankfurt.de

# 3. Feedback

The tutorials are living documents that serve both as test cases to ensure all the
solvers still perform well prior to publishing a new stable release as well as to
minimize the learning curve for new users.  Therefore, user feedback is greatly
appreciated for the benefit of future users.
