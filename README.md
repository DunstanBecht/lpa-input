<div align="center">
  <img width="250" src="https://dunstan.becht.network/views/signatures/mines.svg" alt="Mines Saint-Etienne">
</div>

# Line Profile Analysis - Input

This project is related to the analysis of crystals containing dislocations by X-ray diffraction. It was developed and used for a study conducted during a research internship at the laboratory of material and structural sciences of the *École Nationale Supérieure des Mines de Saint-Étienne*. This repository contains the distribution of one of the three published python packages that have been proposed to conduct line profile analyses based on simulation results:
* [`lpa-input`](https://github.com/DunstanBecht/lpa-input) (line profile analysis input generator)
* [`lpa-xrd`](https://github.com/DunstanBecht/lpa-xrd) (line profile analysis x-ray diffraction simulation program)
* [`lpa-output`](https://github.com/DunstanBecht/lpa-output) (line profile analysis output analyzer)

The repository [`lpa-workspace`](https://github.com/DunstanBecht/lpa-workspace) contains the parameters and the scripts for the generation of the data used in the study. You can then easily replicate the results obtained or use it as inspiration to take the code in hand and conduct your own calculations. The software is placed in the public domain and you can use it as you wish. However, users are encouraged to contribute to the development and report issues.

# Features

The package `lpa.input` can be used to:
* generate dislocation distributions according to different models
* export the distributions in standardized files for input to an X-ray diffraction simulation program
* export the distributions in dislocation maps
* export statistical spatial analyses of the distributions

# Installation

The package is indexed on [PyPI](https://pypi.org/project/lpa-input/) and installable directly via pip:
```bash
pip install -U lpa-input
```

# Examples

### Distribution maps
<div align="center">
<img width="49%" src="https://raw.githubusercontent.com/DunstanBecht/lpa-input/e8440e96ed40d712cf4617009e1117e86617ab48/tests/maps/rho5e13m-2_circle_1000nm_RDD_d5e-5nm-2_S0.svg" alt="RDD">
<img width="49%" src="https://raw.githubusercontent.com/DunstanBecht/lpa-input/e8440e96ed40d712cf4617009e1117e86617ab48/tests/maps/rho5e13m-2_square_2000nm_RRDD-E_s0200nm_f2_S0.svg" alt="RRDD-E">
</div>
<div align="center">
<img width="49%" src="https://raw.githubusercontent.com/DunstanBecht/lpa-input/e8440e96ed40d712cf4617009e1117e86617ab48/tests/maps/rho5e13m-2_circle_2000nm_RCDD-R_d5e-5nm-2_s0500nm_t020nm_IDBC_S0.svg" alt="RCDD-R IDBC">
<img width="49%" src="https://raw.githubusercontent.com/DunstanBecht/lpa-input/e8440e96ed40d712cf4617009e1117e86617ab48/tests/maps/rho5e13m-2_square_2000nm_RCDD-R_d5e-5nm-2_s0500nm_t020nm_PBCG1_S0.svg" alt="RCDD-R PBCG1">
</div>

### Input data files
```
# please keep the structure of this file unchanged
 1  1  0 # z: direction of 'l' (line vector) [uvw]
-1  1  0 # x: direction of 'L' (Fourier variable) [uvw]
 1  1  0 # b: Burgers vector direction [uvw]
 2  0  0 # g: diffraction vector direction (hkl)
0.250000 # C: contrast coefficient [1]
0.404940 # a: cell parameter [nm]
     400 # s: Cylinder radius [nm]
    11.8 # a3: step size of 'L' along x [nm]
   0.345 # nu: Poisson's number [1]
      25 # number of dislocations
# Burgers vector and dislocation (x,y) coordinates
 1 -1.615480114652098E+02 -1.878580545878664E+02
 1 -4.953218476034917E+01  3.963586159012320E+02
 1  3.830929977892965E+02  1.008633196410215E+02
 1  3.294056743470483E+02  3.433102423286451E+01
 1  1.248962568459768E+02 -2.974464812848093E+02
 1  2.832615001129226E+02 -1.729577967167920E+02
 1 -1.955257421410394E+02 -1.549100125777942E+02
 1 -1.888803055949505E+01 -1.458035771343398E+02
 1 -3.270779371095376E+02 -9.196824173171110E+01
 1  2.661326728113684E+02 -1.150221372536317E+02
 1  8.957885498951875E+01 -2.039959039495639E+02
 1  2.787662920148084E+02  4.797066749674389E+00
-1  2.356960793049869E+02 -2.945596912275971E+02
-1  3.780080652647775E+02  8.097447219821595E+01
-1 -3.050154881474039E+01 -2.373117927818869E+02
-1  1.361749016452587E+02  2.700021649518127E+02
-1  1.481167036513976E+02 -1.719318023856983E+02
-1 -2.979589763217694E+02 -7.942577202164772E+01
-1 -7.145235936725004E+01  2.212698723775893E+02
-1 -2.213599709014841E+02  1.168708833157339E+02
-1  3.714582865817929E+02  6.680297635948119E+01
-1  1.354115169218568E+02  1.341973772108176E+02
-1 -1.510359965039211E+02 -2.773050142168986E+02
-1 -6.979442695040841E+01 -9.258073741251835E+01
 1 -2.732036831486906E+02 -2.420388628682046E+02
```

### Spatial analysis
![Ripley’s K function](https://raw.githubusercontent.com/DunstanBecht/lpa-input/4a69a406d1e89073dbd3768c10ab3109aedce3ac/tests/analyses/40000_rho5e13m-2_RRDD-E_s0200nm_f2_circle_1000nm_S0_KKKK.svg)
![Pair correlation function](https://raw.githubusercontent.com/DunstanBecht/lpa-input/4a69a406d1e89073dbd3768c10ab3109aedce3ac/tests/analyses/40000_rho5e13m-2_RRDD-E_s0200nm_f2_circle_1000nm_S0_gggg.svg)
![Symmetric and antisymmetric functions](https://raw.githubusercontent.com/DunstanBecht/lpa-input/4a69a406d1e89073dbd3768c10ab3109aedce3ac/tests/analyses/40000_rho5e13m-2_RRDD-E_s0200nm_f2_circle_1000nm_S0_GaGs.svg)

# Physical aspects

Two geometries are proposed:
* circle (intersection of a plane with a cylinder) centered in `(0, 0)`
* square (intersection of a plane with a cuboid) bottom left corner at `(0, 0)`

A dislocation associates:
* a Burgers vector sense `b`
* a position `p`

A distribution is mainly characterized by the following elements:
* the shape of the region of interest
* the model used for the random generation of dislocations
* the generated dislocations

A sample is a set of distribution and is mainly characterized by:
* the number of generated distribution stored
* the shape of the region of interest
* the model used for the random generation of dislocations
* the stored distributions

# Abbreviations

Some abbreviations are used in the program:

### Distribution models
* `RDD` random dislocation distribution
* `RRDD` restrictedly random dislocation distribution
* `RCDD` random cell dislocation distribution

### Distribution model variants
* `R` randomly distributed Burgers vectors
* `E` evenly distributed Burgers vectors
* `D` dipolar Burgers vectors

### Boundary conditions and considerations
* `ISD` image screw dislocations
* `PBC` periodic boundary conditions
* `GBB` generation beyond boundaries
* `NEC` no edge correction
* `WOA` weighting by overlapping area

# User guide

The directory `tests/` contains several examples of package module usage. The docstrings are carefully written and it is recommended to refer to the documentation with the `help()` command.
