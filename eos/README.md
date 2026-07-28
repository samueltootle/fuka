\page eos Equation of State support
# Equations of State (EOS)

The FUKA solvers support both piecewise polytropic and tabulated EOS in
generating ID for Neutron Stars.  For both cases, a file is required
which includes either the tabulated EOS or describes the piecewise
polytrope, as described below.

## Piecewise Polytropic Equations of State

The file format for storing a piecewise polytrope is quite basic.  All
lines starting with a `#` or are blank are ignored. Most importantly, the
order of parameters is strictly defined

```
#h4 pwpolytrope
num_pieces: 4
rhomin: 1.e-19
rhomax: 1.0
ktab0: 3.99873692e-8
Ptab0: 0

gamma_tab
1.35692395 2.909 2.246 2.144

rho_tab
0.0 8.87824e+14 5.01187234e+14 1.0e+15

#units can be:
#geometrised
#cgs
#cgs_cgs_over_c2
units: cgs_cgs_over_c2
```

1. In the first block is an ordered list (the order matters) of:
    - the number of pieces used to construct the piecewise polytrope
    - the extrema of the allowed densities in geometrized units. Note,
    values above or below the min/max will be set to the min/max value.
    - `K_0` constant factor at the initial boundary and `P0` (pressure at the initial boundary)
2. secondly the values of the polytropic exponent for each piece is
   defined in a single row
3. thirdly the values of the density are defined for each piece within a single row
4. finally the units must be specified for the input values of `K_0` and `rho_tab`

## Tabulated Equations of State

The format for tabulated EOS' that are currently supported is the well
known format used by the [LORENE](https://lorene.obspm.fr/) spectral
solver. Details on the format of the table can be found
[here](https://lorene.obspm.fr/Refguide/classLorene_1_1Eos__tabul.html)

## EOS Framework

The backend that handles the EOS is constructed in two parts

### 1. EOS Management

As for FUKAv2.3, an EOS wrapper was written by Samuel Tootle to interface
between multiple potential EOS backends which currently includes
Margherita and now the [GRHayL](https://github.com/GRHayL/GRHayL).

#### MargheritaEOS
The construction and evaluation of the EOS is handled by default by a
standalone version of the Margherita EOS framework developed by [Elias R.
Most](emost@th.physik.uni-frankfurt.de) For code details see
`$HOME_KADATH/include/EOS/standalone/`

#### GRHayLEOS
In FUKAv2.3, support has been added to link against the GRHayL library to
leverage the GRHayLEOS module. While this includes support for hybid
EOS', it provides additional support to ingest 3D temperature dependent
EOSs in *stellar collapse* format. More information can be found
[here](https://github.com/GRHayL/GRHayL).

### 2. EOS Interface

To interface between KADATH and the EOS backend, a user-defined module
was written by Samuel Tootle and L. Jens Papenfort allowing the use of an
EOS manager within KADATH's System_of_equations framework. For details
see `$HOME_KADATH/include/EOS/EOS.hh`.

To initialize an EOS manager, the following parameters need to be set
within the
[config](https://github.com/samueltootle/fuka/tree/fuka/include/Configurator/)
file

```
  {"eostype",EOSTYPE}, // Options: Cold_PWPoly, Cold_Table, (for GRHayL: ghl_eos_hybrid, ghl_eos_tabulated)
  {"eosfile",EOSFILE}, // Polytrope or Tabulated EOS file
  {"h_cut",HCUT}, // value to cut the specific enthalpy (0 is default)
  {"interpolation_pts", INTERP_PTS} // number of points to use for interpolating the table
```

If no path is specified for the `eosfile`, e.g. `eosfile togashi.lorene`,
the EOS framework will only search the `$HOME_KADATH/eos/` directory for
the EOS file specified.  In the event you want to specify an EOS not
located in the default directory, a path must be included.  For example,
to read an EOS from the current directory, one would need to write
`eos_file ./togashi.lorene`.

### 3. EOS File Formats

#### Polytropes
Both Margherita and GRHayL (`eostype: Cold_PWPoly, ghl_eos_simple, ghl_eos_hybrid`)
support poltryopic equations of state which can
specified in a `<eos>.polytrope` file.  See, e.g. `eos/apr4.polytrope`
and `eos/gam2.polytrope`.

#### Cold Tables
Margherita supports ingesting 1D cold slices in the common `LORENE`
format.  Examples of these tables include `eos/dd2.lorene` and
`eos/togashi.lorene`, for instance.  These tables can be used for
computing initial data with `eostype Cold_Table`.

#### Stellar Collapse
Support for stellar collapse tables is now possible when compiling with
GRHayL. Use of stellar collapse tables requires specifying a
configuration file, e.g. `dd2.stellarcollapse`, which stores the location
of the 3D table along with initialization parameters such as table bounds
and the temperature slice to use (`T_beta`), in units of MeV.

For cold beta equilibrium slices, excellent agreement has been found
between `dd2.lorene` and `dd2.stellarcollapse` computed initial data. Use
of stellar collapse tables to compute initial data using a hot, beta
equilibrium slice has not yet been fully explored.