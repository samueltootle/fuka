\page bnsxcts Binary Neutron Star
# Binary Neutron Star (BNS)

# Overview

Here lies the binary neutron star initial data solver and diagnostic codes.  The initial data is constructed using maximal
slicing with a flat spacial metric with which it is possible to achieve maximal spins up to 
approximately `[-0.6, 0.6]`.  The v2 code is a considerable improvement over the v1 code as it uses superimposed NS solutions
to initialize the binary instead of building the binary from equal mass, non-boosted TOV solutions followed by iterative
changes to the binary components.  When comparing to v1, generating equal-mass non-rotating is automated (i.e. the user no
longer needs to perform the setup and solve separately), however, generating the ID is slightly slower due to needing to perform
a stage of fixed orbital frequency prior to solving the hydro fields consistently.  Most importantly, 
when generating arbitrary spin and unequal mass, the cost savings is roughly (N)x faster 
where N is, in the case of v1, the number of iterative solutions needed to achieve a given mass ratio and spin configurations.

Note:  When referring to `Mtot` below, we will be referring to the sum of the ADM masses of the TOV solution as measured
at infinite separation (i.e. isolated TOV solutions) - `Mtot := (MADM_MINUS + MADM_PLUS)`
where plus and minus simply refer to their location on the x-axis.

# Organization

1. `CMakeLists.txt` is the file needed by CMake to compile the codes
2. `compile` is a symbolic link to the script stored in `$HOME_KADATH/build_release` to ease compiling
3. `src` directory contains the relevant source files:
    - `solve.cpp`: the one and done solve code
    - `reader.cpp`: the reader can provide diagnostics from ID solutions that are computed from the ID

# Basic Usage

1. Generate the initial config file by running `solve` for the first time
2. Modify the initial config file based on ID characteristics you are interested in
3. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bns.info`

# Your first run!

1. Generate the initial config file by running `solve` for the first time
2. Rerun (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bns.info`
3. This will result in the generation of a pair of files containing the solution: `BNS_ECC_RED.togashi.30.2.0.0.2.8.q1.0.0.09.<info/dat>`

We can deconstruct the name to make it understandable:

- `BNS_ECC_RED.` denotes a converged BNS solution after the eccentricity reduction stage is completed which uses 
3.5th order PN estimates for the orbital frequency and radial infall velocity. This is meant to distinguish the solution 
from earlier stages which will be discussed later.This also distinguishes it from checkpoints that can be turned on which 
are saved to file during each iteration of the solver
- `togashi`: the leading name of the eosfile
- `30.2`: separation distance in geometric units!
- `0.0.`: In this case, both NSs are not spinning
- `2.8.q1`: the total mass `2.8` and mass ratio `1`
- `0.0.09`: `nshells = 0` for `ns1`, `nshells = 0` for `ns2` where the resolution in each domain is `9` collocation points in the radial and theta direction with `8` points in the phi direction
- `info`: The info file contains all the steering parameters, values used in creating the domain decomposition, stored fields, stages, settings, and controls
- `dat`: The dat file is a binary file that contains the numerical space and the variable fields

The default configuration for BNS has a total mass of 2.8M and non-spinning.
Aside from the diagnostics observed during the solver stage, we can use the reader
to verify the ID.  This can be done by running:

`./bin/Release/reader BNS_ECC_RED.togashi.28.0.0.2.8.q1.0.0.09.info`

Which results in the following:

```
###################### NS_MINUS ######################
            Center_COM = (-15.10000, 0, 0)
            Coord R_IN = +3.16655
               Coord R = [+5.95732,+6.12761] ([+8.79896,+9.05047] km)
           Coord R_OUT = +9.25540
               Areal R = +7.79075 [+11.50693km]
         Baryonic Mass = +1.55255 (+0.52518,+1.02737,)
     Isolated ADM Mass = +1.40000
      Quasi-local Madm = +1.37754 Diff:+0.01605
         Quasi-local S = +0.00000
                   Chi = +0.00000 [+0.00000]
                 Omega = +0.00004
       x(max(Density)) = -15.10000 (+0.00000)
       Central Density = +1.37425e-03
        Central log(h) = +2.31793e-01
      Central Pressure = +2.34396e-04
    Central dlog(h)/dx = +6.32016e-15
Central Euler Constant = -0.28323
     Integrated log(h) = +186.22988

###################### NS_PLUS ######################
            Center_COM = (+15.10000, 0, 0)
...

###################### Binary ######################
                   RES = [+9,+9,+8]
                     Q = +1.00000
            Separation = +30.20 [+10.79] (+44.61km)
         Orbital Omega = +0.00901
            Komar mass = +2.77644
              Adm mass = +2.77347, Diff: +0.00053
     Total mass (Minf) = +2.80000
           Adm moment. = +7.75769
        Binding energy = -2.65311e-02
            Minf * Ome = +2.52151e-02
            E_b / Minf = -9.47539e-03
                    Px = +5.50470e-15
                    Py = -6.51024e-15
                    Pz = +0.00000e+00
                  COMx = -0.00000, A-COMx = -0.00000
                  COMy = +0.00000, A-COMy = +0.00000
                A-COMz = +0.00000
```

The first two blocks contain information related to the component NSs - the second has been abbreviated since it contains identical information.  
These details are covered in the [NS README](https://bitbucket.org/fukaws/fuka/src/fuka/codes/FUKAv2_Solvers/NS/).
The only additional parameter is the `Center_COM`.  This is the coordinate center of each object when shifted by the
"center-of-mass" of the binary or, more specifically, the location of the axis of rotation for the binary that can approximate a quasi-stationary solution.

The third block contains information specifically related to the binary
- `res`: resolution in each numerical domain `[r,theta,phi]`
- `q`: mass ratio `q <= 1`
- `Separation`: coordinate separation `d [d/Mtot] (d in km)`
- `Orbital Omega`: orbital frequency of the binary
- Komar mass
- ADM quantities
- `E_b / mu`: Binding energy normalized by the reduced mass
- `COM<x/y>`: shift in the coordinates to find a helical killing vector
- `A-COM<x/y/z>`: A numerical calculation of the COM based on the analytical prescription from Osokine+ (REF)

**Note: The `Diff` noted by the ADM mass is the symmetric difference between the ADM and Komar mass.**

# Understanding the BNS INFO file

Using your favorite text editor, you can open up the `initial_bns.info`.  We will go through the file,
but we'll discuss only the details relevant to the BNS case.  For details on all the parameters you can
read more in the [Configurator README](https://bitbucket.org/fukaws/fuka/src/fuka/include/Configurator/).

<b>
Notes: 

1. It is always best practice to generate new ID using the `initial_bns.info`.  Using old initial
data unless for very small changes in `chi` is inefficient.

2. In FUKAv2.2 a minimal Config file was introduced such that only the basic fixing parameters most
relevant to users are shown.  This minimal Config file can be bypassed by running: 
    > `solve full`

    to obtain the full Config file. Although useful for development, there is little advantage to using
the full Config.
</b>

## BNS Fixing parameters

```
binary
{
    distance 30.2
    outer_shells 0
    res 9
    ns1
    {
        chi 0
        madm 1.4
        res 9
        eosfile togashi.lorene
        eostype Cold_Table
    }
    ns2
    {
        chi 0
        madm 1.4
        res 9
        eosfile togashi.lorene
        eostype Cold_Table
    }
}
```

The above includes parameters that can be fixed by the user as well as parameters that are automated in the background
and should not be changed.  The parameters for each NS are simply copied from the isolated solution which is discussed
in detail in the [NS README](https://bitbucket.org/fukaws/fuka/src/fukav2//codes/FUKAv2_Solvers/NS/) - 
the same fixing applies also in the BNS.  The relevant parameters to discuss are

- `res` The resolution shown for the individual compact objects is the highest resolution the *isolated* dataset will be ran at.  This can be important for TOV solutions as the total baryonic mass is sensitive to the resolution.  `res 11` is the minimum recommended for production runs

The fixing parameters most relevant to the binary are


- `distance`: this is in geometric units! It is important to pick something reasonable.  A general rule for a binary with
a few orbits is `distance = 10 * Mtot`, however this strongly depends on `q` and the spins of component objects
- `outer_shells`: This allows for additional shells to be placed near the compactified domain.  This can be helpful
for more accurate quasi-equilibrium ID at lower resolution, but otherwise can be ignored and left to `0`
- `q`: this parameter is computed.  Changing it by hand does nothing
- `res`: global resolution of the binary!
- `adot`: (optional) This is the radial infall velocity parameter when performing eccentricity reduction.  This will be discussed more in the eccentricity reduction section below
- `ecc_omega`: (optional) This is the fixed orbital velocity parameter used when performing eccentricity reduction.  This will be discussed more in the eccentricity reduction section below

## Fields

```
fields
{
    conf on
    lapse on
    shift on
    logh on  ; log specific enthalpy
    phi on   ; fluid velocity potential
}
```

Within the full Config or the output solutions, `Fields` document the fields that are used in the solver and stored in the `dat` file.  Changing this has no impact.

## Stages

```
stages
{
    total on        ; hydrostatic-equilibrium stage
    ecc_red on      ; eccentricity reduction stage - hydro-rescaling with fixed `global_omega`
    total_bc on     ; hydro-rescaling stage with fixed `global_omega`
}
```

## Relevant Controls

```
sequence_controls
{
    centralized_cos off
    checkpoint off
    corot_binary off
    fixed_lapse off
    resolve off
    sequences on
    use_pn off
}
```

- `checkpoint`: this will result in checkpoints being saved to file during each solving iteration - mainly helpful for high 
resolution binary ID where walltimes or server failures are a concern prior to a converged solution being obtained
- `corot_binary`: the objects are no longer fixed based on `chi` and instead provide a corotating ID solution
- `fixed_lapse`: toggling this control enables a fixed lapse on the horizon - not recommended
- `sequences`: this toggle is enabled by default and essentially tells the driver routine to start from scratch.  
If this is enabled when attempting to use a previous solution, the previous fields and numerical space (i.e. the `dat` file) is ignored
- `use_pn`: toggle whether to always use 3.5PN estimates.  It is important to ensure this is off if the user wants to specify their own `adot` and `global_omega`
parameters by hand (e.g. for iterative eccentricity reduction)
- `resolve`: force all implicit compact object solutions to be resolved regardless of an existing previous solution
- `centralized_cos`: stores all implicit COs into `$HOME_KADATH/COs`

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
- `initial_resolution`: the resolution to use when solving the binary initial before regridding to a higher
resolution

# Your Second run!

Now that you've generated the simplest case and we have a better understanding of the config file, we can try something more interesting

1. Open the initial config file in your favorite editor
1. Set `distance 30.2`
1. (optional) set `res 11`
2. For `ns1` set:
    - `madm 1.18` 
    - `chi 0`
    - `res 11`
3. For `ns2` set:
    - `madm 2.42` 
    - `chi 0.52`
    - `res 11`
3. Run (using parallelization) using this config file, e.g. `mpirun ./bin/Release/solve initial_bns.info`

This time around we see the iterative `chi` increase being done for the primary NS as well as a regrid of the solution
to the higher resolution before being imported into the initial binary setup. Overall, the main changes 
observed are related to the isolated NS solvers 
(see the [NS README](https://bitbucket.org/fukaws/fuka/src/fukav2//codes/FUKAv2_Solvers/NS/) for details), 
but the binary solver itself is consistent when compared to the equal mass case.

In the event you changed the resolution to 11pts, the solver will solve the binary at the `initial_resolution` until a
solution in hydrostatic equilibrium has been obtained.  At this point the solution will be regridded to a higher resolution
at which point the solver will go through all the activated stages (e.g. `TOTAL`, `TOTAL_BC`, and `ECC_RED`).

This results in the converged dataset of `BNS_ECC_RED.togashi.28.0.52.0.3.6.q0.487603.0.0.11.info/dat`, however, the other implicit solutions have been saved as well:

1. `BNS_TOTAL_FIXED_OMEGA.`: is the initial solution after the import of the two isolated solutions have been solved in the binary space for a fixed COM and orbital frequency.  The hydro fields are simply rescaled to enforce the
specified baryonic mass, but the fluid is not in hydrostatic equilibrium
2. `BNS_TOTAL.`: this is the solution in complete hydrostatic equilibrium with the ADM linear momenta and the orbital frequency being fixed by the force-balance equations for each NS
3. `BNS_TOTAL_BC.`: is the iterative solution such that the orbital frequency is set to a constant and the remaining linear momenta are driven to zero by varying the COM. In this stage, the fluid deviates from hydrostatic equilibrium since the orbital frequency is a constant
4. `BNS_ECC_RED.`:  The final solution is one where the orbital frequency and radial infall velocity is fixed to either 3.5th order PN estimates based on the COM obtained in the `TOTAL_BC` stage or the values for `adot` and `ecc_omega` are used in the case of iterative eccentricity reduction

We can of course verify that the ID matches our expectation using

`./bin/Release/reader BNS_ECC_RED.togashi.30.2.0.0.52.3.6.q0.487603.0.0.09.info`

```
###################### NS_MINUS ######################
            Center_COM = (-20.30551, 0, 0)
            Coord R_IN = +3.27559
               Coord R = [+5.85351,+6.22200] ([+8.64563,+9.18989] km)
           Coord R_OUT = +9.40079
               Areal R = +7.74194 [+11.43484km]
         Baryonic Mass = +1.28557 (+0.47860,+0.80697,)
     Isolated ADM Mass = +1.18000
      Quasi-local Madm = +1.15688 Diff:+0.01960
         Quasi-local S = +0.00000
                   Chi = +0.00000 [+0.00000]
                 Omega = +0.00010
       x(max(Density)) = -15.10000 (+0.00000)
       Central Density = +1.22102e-03
        Central log(h) = +1.83864e-01
      Central Pressure = +1.57833e-04
    Central dlog(h)/dx = -4.68108e-16
Central Euler Constant = -0.29379
     Integrated log(h) = +142.05251

###################### NS_PLUS ######################
            Center_COM = (+9.89449, 0, 0)
            Coord R_IN = +1.88069
               Coord R = [+3.68691,+4.64383] ([+5.44556,+6.85893] km)
           Coord R_OUT = +9.40079
               Areal R = +7.20294 [+10.63873km]
         Baryonic Mass = +2.90870 (+0.84142,+2.06728,)
     Isolated ADM Mass = +2.42000
      Quasi-local Madm = +2.35065 Diff:+0.02866
         Quasi-local S = +3.04533
                   Chi = +0.52000 [+0.52000]
                 Omega = +0.05413
       x(max(Density)) = +15.10000 (+0.00000)
       Central Density = +2.70298e-03
        Central log(h) = +7.61540e-01
      Central Pressure = +2.10043e-03
    Central dlog(h)/dx = +2.41079e-16
Central Euler Constant = -0.65064
     Integrated log(h) = +538.45907

###################### Binary ######################
                   RES = [+9,+9,+8]
                     Q = +0.48760
            Separation = +30.20 [+8.39] (+44.61km)
         Orbital Omega = +0.00987
            Komar mass = +3.56746
              Adm mass = +3.56919, Diff: +0.00024
     Total mass (Minf) = +3.60000
           Adm moment. = +13.28549
        Binding energy = -3.08070e-02
            Minf * Ome = +3.55439e-02
            E_b / Minf = -8.55751e-03
                    Px = -3.84335e-14
                    Py = +6.42085e-14
                    Pz = +0.00000e+00
                  COMx = -5.20551, A-COMx = +5.21435
                  COMy = +0.05651, A-COMy = +0.00037
                A-COMz = +0.00000
```

# How BNS ID is Generated

## Initial Setup

To generate the initial setup for a BNS ID, we first need to make some guesses based on the input 
ADM masses and the coordinate separation

  1. `COM` - An estimate of the center-of-mass: this is purely Newtonian
  2. `global_omega` - An estimate of the orbital frequency: 3.5th PN estimate

Once these estimates are computed an interface code is ran
which 

- solves each NS configuration in isolation (see the [NS README](https://bitbucket.org/fukaws/fuka/src/fuka/codes/FUKAv2_Solvers/NS/) for more details).  
- obtains boosted isolated solutions using the estimated `global_omega`

At this point, the binary numerical space and fields are constructed and the isolated solutions are interpolated onto 
the new grid using the idea of superimposed solutions.  Specifically:

- a decay parameter `decay_limit := w` is chosen such that `w = distance / 2`
- the fields are interpolated such that the solutions decay exponetially away from each object as 
`decay_rate = exp(-(r_NS / w)^4)` where `r_NS` is the coordinate distance to the respective NS
- The resulting value at a given point is then simply the sum of the background with the deviations computed from the isolated solutions

For example, if we wanted to compute the initial guess for the lapse at a given point `x`, it would be

`lapse(x) = 1. + decay_rate_NS1 * (lapse_NS1(x) - 1.) + decay_rate_NS2 * (lapse_NS2(x) - 1.)`

This is then repeated for all fields in all numerical domains, with the compactified domain set the asymptotic 
values of the fields, i.e. `lapse = psi = 1`, `log(h) = shift = 0`.

## Initial Solution

Once the initial guess has been setup for the BNS, the first solver stage solves the full XCTS system of equations consistently 
and all using a fixed orbital velocity.  As a consequence the matter is simply rescaled to achieve the desired baryonic mass
for each NS.  The reason for this is two fold:

  1. The spacetime fields need to settle to a more accurate estimate before resolving the matter consistently
  2. The fluid velocity potential field `phi` needs to initialize more accurately to the binary configuration

The output file from this stage, `BNS_TOTAL_FIXED_OMEGA.`, can readily be discarded once a solution from a later
stage is obtained.

## Hydrostatic Equilibrium

Now with a more stable background, the matter becomes a variable field and is allowed to solve consistently in order to obtain
hydrostatic equilibrium in the binary configuration.  In this stage, the force balance equations are introduced which provides
a quasi-equilibrium estimate of `global_omega` and one component of the ADM linear momenta is minimized.

**Note: in the event one wants to later increase the resolution or make iterative changes (e.g. make small changes to the MADM, MB, CHI of one or both stars),
it can only be done using the solution from this stage, `BNS_TOTAL.*<info/dat>`, as the initial starting point.  Reusing the solutions from the hydro-rescaling stage
will more often than not cause the solution to diverge or lead to unphysical results in the numerical
evolution.  Therefore, these solutions can be useful to retain.**
to retain.**

### Automated resolution increase

In the event the binary resolution was set to something higher than the sequence_setting `initial_resolution`, the automated
increase resolution will take place here and resolve the binary in hydro-static equilibrium.  This is very important as
increasing the resolution and repeating a matter-rescaling stage will result in a very inconsistent description of the fluid
as it will include numerical errors from the interpolated solution at lower resolution!

## Minimize Linear Momentum

In this stage we return to simply rescaling the matter now that it has obtained a more consistent solution
and we fix it to the current `global_omega`, however, we now minimize the x and y linear momentum components using the force-balance
equations to determine the most accurante COM for the given binary.  This is very important prior to the eccentricity reduction
stage as a precise estimate of the COM is required for the post-Newtonian estimates.

The output file from this stage, `BNS_TOTAL_BC.`, can be discarded when no longer needed.

## Eccentricity Reduction

Using either 3.5PN estimates or, in the case of iterative eccentricity reduction, `ecc_omega` and `adot` from the config file,
a final stage of matter rescaling is performed based on the changes introduced by `ecc_omega` and `adot`.  This is the recommended
solution to use for evolutions and it is stored with a filename title of `BNS_ECC_RED.*<info/dat>`
