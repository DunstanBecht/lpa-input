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
 1 -1.61548011465209754078E+02 -1.87858054587866405427E+02
 1 -4.95321847603491747236E+01  3.96358615901232042233E+02
 1  3.83092997789296475730E+02  1.00863319641021547568E+02
 1  3.29405674347048261552E+02  3.43310242328645074394E+01
 1  1.24896256845976836303E+02 -2.97446481284809294721E+02
 1  2.83261500112922590233E+02 -1.72957796716792046254E+02
 1 -1.95525742141039444277E+02 -1.54910012577794219624E+02
 1 -1.88880305594950499426E+01 -1.45803577134339775512E+02
 1 -3.27077937109537572269E+02 -9.19682417317111031707E+01
 1  2.66132672811368365728E+02 -1.15022137253631683507E+02
 1  8.95788549895187458105E+01 -2.03995903949563853530E+02
 1  2.78766292014808357180E+02  4.79706674967438928547E+00
-1  2.35696079304986938041E+02 -2.94559691227597113539E+02
-1  3.78008065264777542325E+02  8.09744721982159489926E+01
-1 -3.05015488147403850405E+01 -2.37311792781886936154E+02
-1  1.36174901645258700000E+02  2.70002164951812687832E+02
-1  1.48116703651397642716E+02 -1.71931802385698262015E+02
-1 -2.97958976321769398510E+02 -7.94257720216477167696E+01
-1 -7.14523593672500396679E+01  2.21269872377589337020E+02
-1 -2.21359970901484132355E+02  1.16870883315733863128E+02
-1  3.71458286581792890502E+02  6.68029763594811925032E+01
-1  1.35411516921856758700E+02  1.34197377210817563764E+02
-1 -1.51035996503921097656E+02 -2.77305014216898598534E+02
-1 -6.97944269504084076061E+01 -9.25807374125183457636E+01
 1 -2.73203683148690629423E+02 -2.42038862868204574852E+02
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

### Models
* `RDD` random dislocation distribution
* `RRDD` restrictedly random dislocation distribution
* `RCDD` random cell dislocation distribution

### Model variants
* `R` randomly distributed Burgers vectors
* `E` evenly distributed Burgers vectors
* `D` dipolar Burgers vectors

### Boundary conditions
* `PBCG` periodic boundary conditions applied when generating the distribution
* `PBCR` periodic boundary conditions applied when running the simulation
* `IDBC` image dislocations boundary conditions

# User guide

The directory `tests/` contains several examples of package module usage. The docstrings are carefully written and it is recommended to refer to the documentation with the `help()` command.
