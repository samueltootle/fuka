\page bhxcts Black Hole
# Black Hole (BH) Initial Data

# Overview

The isolated BH code was created to test many aspects of the initial data creation such as
formalism implementation, initial guess generation, and import into an evolution framework.
It has now become a critical piece in generating binary ID that include a BH companion.
The codes included here are simply the means for a user to generate an isolated BH and analyze it.  
This is an excellent way for a new user to get acquainted with how the FUKAv2 solvers work 
since the interface is similar amongst the v2 solvers by design.

# Organization

1. `CMakeLists.txt` is the file needed by CMake to compile the codes
2. `compile` is a symbolic link to the script stored in `$HOME_KADATH/build_release` to ease compiling
3. `src` directory contains the relevant source files:
    - `solve.cpp`: the one and done solve code
    - `reader.cpp`: the reader can provide diagonstics from ID solutions that are computed from the ID

# Basic Usage

1. Generate the initial config file by running `solve` for the first time
2. Modify the initial config file based on ID characteristics you are interested in
3. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bh.info`

# Your first run!

1. Generate the initial config file by running `solve` for the first time
2. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bh.info`
3. This will result in the generation of a pair of files containing the solution: `BH_TOTAL_BC.0.5.0.09.<info/dat>`

We can deconstruct the name to make it understandable:

- `BH_TOTAL_BC.` denotes a converged BH solution after the `TOTAL_BC` stage.
        This is meant to distinguish the solution from other stages suched as the boosted BH stage.
        This also distinquishes it from checkpoints that can be turned on which are saved to file during each iteration of the solver
- `0.5.0.`: In this case, we have a 0.5M BH with a dimensionless spin of zero
- `0.09.`: No additional spherical shells have been added and the global resolution is nine collocation points in the radial and theta direction with 8 points in the phi direction
- `info`: The info file contains all the steering parameters, values used in creating the domain decomposition, stored fields, stages, settings, and controls
- `dat`: The dat file is a binary file that contains the numerical space and the variable fields

The default configuration for an isolated BH has a mass of 0.5M and non-spinning.
Aside from the diagnostics observed during the solver stage, we can use the reader
to verify the ID.  This can be done by running:

`./bin/Release/reader BH_TOTAL_BC.0.5.0.0.09.info`

Which results in the following:

```
                   RES = [+9,+9,+8]
            Coord R_IN = +0.20812
               Coord R = +0.42977
           Coord R_OUT = +1.66493

               Areal R = +1.00000
                 LAPSE = [+0.43494, +0.43494]
                   PSI = [+1.52539, +1.52539]

                  Mirr = +0.50000[+0.50000]
                   Mch = +0.50000[+0.50000]
                   Chi = -0.00000[+0.00000]

                 Omega = -3.02175e-17
            Komar mass = +4.97314e-01
              Adm mass = +5.00087e-01, Diff: +5.56053e-03
           Adm moment. = +3.14400e-14
                    Px = +1.22876e-15
                    Py = -2.36096e-16
                    Pz = +0.00000e+00
```

The first block contains information related to the numerical space specifically the 

- resolution `[r,theta,phi]`
- the coordinate radii related to the domain decomposition.

The second block includes the 

- proper radius of the horizon 
- `[min, max]` of the lapse and conformal factor on the horizon.

The third block contains the 

- measured values of the fixing parameters: 
    - Christodoulou mass (MCH)
    - irreducible mass (MIRR)
    - the dimensionless spin parameter (CHI)
    
- brackets are the values stored in the config file that the solution was fixed by

Finally, the fourth block contains the 

- variable spin frequency (Omega) 
- Komar mass
- ADM quantities

<b>
Note: 
1. In all blocks, brackets denote the values stored in the config file that the solution was fixed by
2. The `Diff` noted by the ADM mass is the symmetric difference between the ADM and Komar mass.

</b>

# Understanding the BH INFO file

Using your favorite text editor, you can open up the `initial_bh.info`.  We will go through the file,
but we'll discuss only the details relevant to the BH case.  For details on all the parameters you can
see more in [Configurator README](https://bitbucket.org/fukaws/fuka/src/fuka/include/Configurator/).

<b>
Notes: 

1. It is always best practice to generate new ID using the `initial_bh.info`.  Using old initial
data unless for very small changes in `chi` is inefficient.

2. In FUKAv2.2 a minimal Config file was introduced such that only the basic fixing parameters most
relevant to users are shown.  This minimal Config file can be bypassed by running: 
    
    solve full

to obtain the full Config file. Although useful for development, there is little advantage to using
the full Config.

</b>

## BH Fixing parameters

```
bh
{
    chi 0
    mch 0.5
    res 9
}
```

The above includes parameters that can be fixed by the user as well as parameters that are automated in the background
and should not be changed.  Those most relevant to generating ID of interest are

- `chi`: the dimensionless spin parameter `[-0.84, 0.84]`
- `mch`: the Christodoulou mass
- `res`: the final global resolution of interest (this will be discussed more below)

## Fields

```
fields
{
    conf on
    lapse on
    shift on
}
```

Within the full Config or the output solutions, `Fields` document the fields that are used in the solver and stored in the `dat` file.  Changing this has no impact.

## Stages

```
stages
{
    total_bc on
}
```

For the isolated BH, only the `TOTAL_BC` is relevant.  Old stages available in the original v1 solver are since deprecated.

## Relevant Controls

```
sequence_controls
{
    checkpoint off
    fixed_lapse off
}
```

- `checkpoint`: this will result in checkpoints being saved to file during each solving iteration - mainly helpful for high resolution binary ID
- `fixed_lapse`: toggling this control enables a fixed lapse on the horizon - not recommended

## Sequence Settings

```
sequence_settings
{
    solver_max_iterations 15
    solver_precision 1e-08
    initial_resolution 9
}
```

- `solver_max_iterations`: set the number of iterations not to exceed
- `solver_precision`: Determines what is the maximum precision allowed by the solver that determines whether or not the solution has converged
- `initial_resolution`: by default all solutions are ran at at a default resolution of `[9,9,8]` prior to regridding to a higher resolution.  In the event one wants to increase the default `initial_resolution`, it can be done here

# Your Second run!

Now that you've generated the simplest case and we have a better understanding of the config file, we can try something more interesting

1. Open the initial config file in your favorite editor
2. Set
    - `mch 1` 
    - `chi 0.84`
    - `res 11`
3. Run (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bh.info`

Many things will run differently this time!  Most noteably are the following:

1. Even though the resolution was set to `res 11`, the initial solution is always constructed using `initial_resolution` value first.
Currently the default initial resolution is `[9,9,8]`.
2. Given the `chi` used is close to the maximum possible with the maximally sliced, flat background XCTS, an interative procedure
is done automatically to avoid the solution diverging.  This is done in three steps to ensure convergence
    
    1. `chi = 0.5`
    2. `chi = 0.8`
    3. `chi = 0.84`

3. Once the desired `chi` is achieved, the solution is regridded to a new grid with the desired final resolution
4. With the new initial guess based on the low resolution solution, a final solving stage is conducted to obtain the desired ID
5. This results in the converged dataset of `BH_TOTAL_BC.1.0.84.0.11.info/dat`, however, the other implicit solutions have been saved as well.

We can of course verify that the ID matches our expectation using

`./bin/Release/reader BH_TOTAL_BC.1.0.84.0.11.info`

```
                   RES = [+11,+11,+10]
            Coord R_IN = +0.32502
               Coord R = +0.54163
           Coord R_OUT = +2.60012

               Areal R = +1.75647
                 LAPSE = [+0.28547, +0.32448]
                   PSI = [+1.75317, +1.82613]

                  Mirr = +0.87823[+0.87823]
                   Mch = +1.00000[+1.00000]
                   Chi = +0.84000[+0.84000]

                 Omega = -3.03075e-01
            Komar mass = +1.00517e+00
              Adm mass = +1.00629e+00, Diff: +1.11392e-03
           Adm moment. = +8.39740e-01
                    Px = +4.19768e-15
                    Py = -3.99361e-15
                    Pz = +0.00000e+00
```

# How BH ID is Generated

## Initial Setup

In initial data generation, the initial setup is the most crucial.  When starting the BH solver from scratch we need to determine
an initial guess for the input mass and aspin - `mch`, `chi`.  For this a crude estimate is generated based on a fixed value of the conformal
factor of `conf = 1.55` to determine the domain decomposition and initial guess of the excision radius of 

    `rmid = 2 * mch / conf^2`

In the domain near the excision region, `conf = 1.55` and `conf = 1` everywhere else. `lapse = 1` and `shift = 0`.

## TOTAL_BC Stage

With the initial guess generated, a rotating solution is now found based on the input mass `mch` and dimensionless spin `chi` fixing parameters.

## Automated iterative spin

For highly spinning BH configurations using a conformally flat background, it has been found to be necessary to obtain lower spinning solutions prior to a highly spinning one.  Therefore, an initial spin of `chi = 0.5` is used before attempting spins up to `chi = 0.8`.  For spins near the XCTS spin limit of approximately `chi = 0.84`, an intermediate solution
is obtained at `chi = 0.8` before attempting `chi > 0.8`.

## Automated resolution increase

Once the above procedures are completed for the last stage and the input `res` is higher than the `initial_resolution`, 
the low resolution solution is interpolated onto the higher resolution numerical grid and is used as the initial guess before 
running the last stage again.  No additional iterative procedures are required at this point.
