\page nsxcts Neutron Star
# Neutron Star (NS) Initial Data

# 1. Overview

The isolated NS code was created to test many aspects of the initial data creation such as
formalism implementation, initial guess generation, and import into an evolution framework.
It has now become a critical piece in generating binary ID that include a NS companion.
The codes included here are simply the means for a user to generate an isolated NS and analyze it.  
This is an excellent way for a new user to get acquainted with how the FUKAv2 solvers work 
with ID that requires an EOS.

# 2. Organization

1. `CMakeLists.txt` is the file needed by CMake to compile the codes
2. `compile` is a symbolic link to the script stored in `$HOME_KADATH/build_release` to ease compiling
3. `src` directory contains the relevant source files:
    - `solve.cpp`: the one and done solve code
    - `reader.cpp`: the reader can provide diagnostics from ID solutions that are computed from the ID

# 3. Base Usage

1. Generate the initial config file by running `solve` for the first time
2. Modify the initial config file based on ID characteristics you are interested in
3. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_ns.info`

# 4. Your first run!

1. Generate the initial config file by running `solve` for the first time
2. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_ns.info`
3. This will result in the generation of a pair of files containing the solution: 
`NS_TOTAL_BC.togashi.1.4.0.0.09.<info/dat>`

We can deconstruct the name to make it understandable:

- `NS_TOTAL_BC.` denotes a converged NS solution after the TOTAL_BC stage.
        This is meant to distinguish the solution from other stages such as the NOROT_BC and boosted NS stage.
        This also distinguishes it from checkpoints that can be turned on which our saved to file during each iteration of the solver
- `togashi` denotes the eosfile used in the ID construction
- `1.4.0.`: In this case, we have a 1.4M NS with a dimensionless spin of zero
- `0.09.`: the number of *additional* spherical `nshells` is `0` and the resolution in each domain is `9` collocation points in the radial and theta direction with `8` points in the phi direction
- `info`: The info file contains all the steering parameters, values used in creating the domain decomposition, 
EOS parameters, stored fields, stages, settings, and controls
- `dat`: The dat file is a binary file that contains the numerical space and the variable fields

The default configuration for an isolated NS has a mass of 1.4M and non-spinning.
Aside from the diagnostics observed during the solver stage, we can use the reader
to verify the ID.  This can be done by running:

`./bin/Release/reader NS_TOTAL_BC.togashi.1.4.0.0.09.info`

Which results in the following:

```
                   RES = [+9,+9,+8]
            Coord R_IN = +3.03469
               Coord R = [+6.33309, +6.33309]
           Coord R_OUT = +12.13876

               Areal R = +7.81046 [+11.53605km]
         Baryonic Mass = +1.55255
              ADM Mass = +1.40000 [+1.40000]
          ADM Momentum = +0.00000
                   Chi = +0.00000 [+0.00000]
                 Omega = +0.00000
       Central Density = +1.37405e-03
        Central log(h) = +2.31726e-01
      Central Pressure = +2.34280e-04
    Central dlog(h)/dx = +3.40893e-14
Central Euler Constant = -2.22137e-01
     Integrated log(h) = +186.08372

                    Mk = +1.40096, Diff: +6.84955e-04
                    Px = +0.00000
                    Py = +0.00000
                    Pz = +0.00000
```

The first block contains information related to the numerical space specifically the 

- resolution `[r,theta,phi]`
- the fixed coordinate radii related to the nucleus
- the extremal coordinate radii along the adaptive surface
- the fixed coordinate radii related to the spherical boundaries outside the stellar surface

The second block includes the 

- proper radius of the NS as measured on the stellar surface
- total baryonic mass
- ADM mass
- ADM spin angular momentum
- `chi` the dimensionless spin parameter
- variable spin frequency `omega`
- Central density
- Central log specific enthalpy `h`
- Central pressure
- Central dlog(h)/dx
- Euler constant at the steller center
- Integral of log specific enthalpy over the whole NS

Finally, the third block includes

- Komar mass
- ADM linear momenta `P<x,y,z>`

<b>
Note: 

1. In all blocks, brackets denote the values stored in the config file that the solution was fixed by
2. The `Diff` noted by the Komar mass is the symmetric difference between the ADM and Komar mass.

</b>

# 5. Understanding the NS INFO file

Using your favorite text editor, you can open up the `initial_ns.info`.  We will go through the file,
but we'll discuss only the details relevant to the NS case.  For details on all the parameters you can
see more in [Configurator README](https://bitbucket.org/fukaws/fuka/src/fuka/include/Configurator/).

<b>
Notes: 

1. It is always best practice to generate new ID using the `initial_ns.info`.  Using old initial
data unless for very small changes in `chi` is inefficient.
2. In FUKAv2.2 a minimal Config file was introduced such that only the basic fixing parameters most
relevant to users are shown.  This minimal Config file can be bypassed by running: 
    
    solve full

to obtain the full Config file. Although useful for development, there is little advantage to using
the full Config.

</b>

## NS Fixing parameters

```
ns
{
    chi 0
    madm 1.3999999999999999
    res 9
    eosfile togashi.lorene
    eostype Cold_Table
}
```

The above includes parameters that can be fixed by the user as well as parameters that are automated in the background
and should not be changed.  Those most relevant to generating ID of interest are

- `chi`: the dimensionless spin parameter `~[-0.6, 0.6]`
- `madm`: the total gravitational mass is a fixing parameter iff  the `mb_fixing` control is set to off
- `res`: the final resolution used in each numerical domain (this will be discussed more below)

## EOS Parameters

- `eosfile`: The EOS file containing the cold table or piecewise polytrope to be used
- `eostype`: specify the use of a cold table `Cold_Table` or a piecewise polytrope `Cold_PWPoly`
- `h_cut`: (optional) in the event one wants to cut the cold table at a certain specific enthalpy value (default: 0)
- `interpolation_pts`: (optional) number of points to use in order to interpolate a cold table (default: 2000).  Not relevant for `Cold_PWPoly`

For example, if one wanted to change from the default EOS to a `Gamma = 2` polytrope, they would set in the config file

- `eosfile gam2.polytrope`
- `eostype Cold_PWPoly`

The other parameters can be ignored.

## Fields

```
fields
{
    conf on
    lapse on
    shift on
    logh on
}
```

Within the full Config or the output solutions, `Fields` documents the fields that are used in the solver and stored in the `dat` file.  Changing this has no impact.

## Stages

```
stages
{
    norot_bc on
    total_bc on
}
```

For the isolated NS, only the `NOROT_BC` and `TOTAL_BC` stages are relevant.  These will be discussed in the
sections below related to these stages.
Old stages available in the original v1 solver are since deprecated.

## Relevant Controls

```
sequence_controls
{
    checkpoint off
    sequences on
}
```

- `checkpoint`: this will result in checkpoints being saved to file during each solving iteration - mainly helpful for high resolution binary ID
- `sequences`: this toggle is enabled by default and essentially tells the driver routine to start from scratch.  If this is enabled when attempting to use a previous solution, the previous solution (i.e. the `dat` file) is ignored and the solver starts from the most basic initial guess
- `fixed_mb`: **(not fully tested)** toggling this control enables using fixed baryonic mass `mb` instead of fixing the gravitational mass `madm`

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
- `initial_resolution`: (optional) by default all solutions are ran at at a default resolution of `[9,9,8]` prior to regridding to a higher resolution.  In the event one wants to increase the default `initial_resolution`, it can be done here


# 6. Your Second run!

Now that you've generated the simplest case and we have a better understanding of the config file, we can try something more interesting

1. Open the initial config file in your favorite editor
2. Set
    - `madm 2.3` 
    - `chi 0.6`
    - `res 11`
3. Run (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_ns.info`

Many things will run differently this time!  Most notably are the following:

1. Even though the resolution was set to `res 11`, the initial solution is always constructed using `initial_resolution` first.
Currently the default initial resolution is `res 9`.
2. For rotating solutions, two automated procedures are in place to ensure convergence
    1. Iterative mass increase: in the event the desired mass is above Mtov 
        
        - the initial spinning solution is obtained using Mtov
        - once the desired spin is obtained at Mtov, the mass is increased to the desired value and the `TOTAL_BC` stage is reran

    2. Iterative spin increase: the system of equations is sensitive to large changes in the spin, therefore, an automated iterative
procedure is done to maximize the chance for convergence by solving for

        1. `chi = 0.1`
        2. `chi = 0.6`

3. Once the desired `chi` and `madm` is achieved, the solution is interpolated onto a new numerical grid with the desired final resolution
4. With the high resolution initial guess, the `TOTAL_BC` stage is reran to obtain the desired ID
5. This results in the converged dataset of `NS_TOTAL_BC.togashi.2.3.0.6.0.0.11.info/dat`, 
however, the other implicit solutions have been saved as well.


We can of course verify that the ID matches our expectation using

`./bin/Release/reader NS_TOTAL_BC.togashi.2.3.0.6.0.11.info`

```
                   RES = [+11,+11,+10]
            Coord R_IN = +1.47139
               Coord R = [+2.94504, +3.97120]
           Coord R_OUT = +5.95079

               Areal R = +6.25539 [+9.23921km]
         Baryonic Mass = +2.65411
              ADM Mass = +2.30000 [+2.30000]
          ADM Momentum = +3.17400
                   Chi = +0.60000 [+0.60000]
                 Omega = +0.06396
       Central Density = +4.56602e-03
        Central log(h) = +1.46938e+00
      Central Pressure = +1.02698e-02
    Central dlog(h)/dx = +6.70501e-14
Central Euler Constant = -7.39008e-01
     Integrated log(h) = +583.85557

                    Mk = +2.29825, Diff: +7.63310e-04
                    Px = +0.00000
                    Py = +0.00000
                    Pz = +0.00000
```

# 7. How NS ID is Generated

## Initial Setup

In initial data generation, the initial setup is the most crucial.  When starting the NS solver from scratch, we need to determine
an initial guess and what better way than a 1D TOV solution!  In this way, a root finder is attempted for the input ADM mass for the 
specific EOS.  Either a solution is found for your input mass (i.e. the input mass is less than the maximum non-rotating TOV 
mass, Mtov) or the maximum mass is used as the initial guess for the non-rotating solution.  
In the event an `madm` above Mtov is used with no spin, an error is given.

Once the 1D TOV solution is obtained, this solution is then interpolated onto the scalar fields associated with the low resolution 3D numerical grid.

## NOROT_BC Stage

Now that the initial setup is complete, the first stage to run is the non-rotating stage.  In this stage, all equations related to the shift are not taken into account while the fluid and space-time are able to find an equilibrium solution.

## TOTAL_BC Stage

With a stable solution to start from, a rotating solution is now found based on the `chi` fixing parameter.

### Automated iterative spin

For highly spinning NS configurations it has been found to be most reliable to obtain a slowly spinning solution prior to a highly spinning one.  Therefore, an initial spin of `chi = +/- 0.1` is used before attempting higher spins.

### Automated iterative mass

In the event the desired mass is above Mtov, the initial solution is constructed up to the desired `chi` using Mtov.  Once
this solution has been obtained, the ADM mass is updated to the desired ADM mass and only the `TOTAL_BC` stage is reran

### Automated resolution increase

Once the above procedures are completed for the last stage and the input `res` is higher than the `initial_resolution`, 
the low resolution solution is interpolated onto the higher resolution numerical grid and is used as the initial guess before 
running the last stage again.  No additional iterative procedures are required at this point.
