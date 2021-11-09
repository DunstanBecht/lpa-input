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
<img width="49%" src="https://raw.githubusercontent.com/DunstanBecht/lpa-input/fc3d0ca5296e33fdcc5f51d7f4c19a581c2d6714/tests/maps/rho5e13m-2_circle_1000nm_RDD_d5e-5nm-2_S0.svg" alt="RDD">
<img width="49%" src="https://raw.githubusercontent.com/DunstanBecht/lpa-input/fc3d0ca5296e33fdcc5f51d7f4c19a581c2d6714/tests/maps/rho5e13m-2_square_2000nm_RRDD-E_s0200nm_f2_S0.svg" alt="RRDD-E">
</div>
<div align="center">
<img width="49%" src="https://raw.githubusercontent.com/DunstanBecht/lpa-input/fc3d0ca5296e33fdcc5f51d7f4c19a581c2d6714/tests/maps/rho5e13m-2_circle_2000nm_RCDD-R_d5e-5nm-2_s0200nm_t020nm_ISD_S0.svg" alt="RCDD-R ISD">
<img width="49%" src="https://raw.githubusercontent.com/DunstanBecht/lpa-input/fc3d0ca5296e33fdcc5f51d7f4c19a581c2d6714/tests/maps/rho5e13m-2_square_2000nm_RCDD-R_d5e-5nm-2_s0200nm_t020nm_PBC1_S0.svg" alt="RCDD-R PBC1">
</div>

### Input data files
```
   1.2.1 # v: lpa-input version
5.17E+13 # d: dislocation density [m^-2]
 1  1  0 # z: direction of 'l' (line vector) [uvw]
-1  1  0 # x: direction of 'L' (Fourier variable) [uvw]
 1  1  0 # b: Burgers vector direction [uvw]
 2  0  0 # g: diffraction vector direction (hkl)
0.250000 # C: contrast coefficient [1]
0.404940 # a: cell parameter [nm]
     400 # s: Cylinder radius [nm]
    11.6 # a3: step size of 'L' along x [nm]
   0.345 # nu: Poisson's number [1]
      26 # number of dislocations in this file
# Burgers vector senses and dislocation (x,y) coordinates [1], [nm], [nm]
 1 -2.604424178510505E+02 -3.028586084465366E+02
 1 -4.912383207593395E+01  3.930909606266370E+02
 1  3.202751109125161E+02  8.432420083752132E+01
 1  3.208662976752835E+02  3.344104093784984E+01
 1  1.284915410549065E+02 -3.060088246581881E+02
 1  2.129037512154038E+02 -1.299977713465352E+02
 1 -1.152378285855275E+02 -9.129996531477097E+01
 1 -4.364950050196162E+01 -3.369463689327758E+02
 1 -2.791018700891352E+02 -7.847826265191925E+01
 1  2.045136262814916E+02 -8.839047886863894E+01
 1  1.120985146258099E+02 -2.552794163887354E+02
 1  3.771948070067862E+02  6.490844548543851E+00
 1  2.415271223288822E+02 -3.018470005359157E+02
-1  2.339563106390758E+02  5.011662584016895E+01
-1 -3.854997415041839E+01 -2.999311127738565E+02
-1  1.021921315581637E+02  2.026224835001035E+02
-1  2.012644520339315E+02 -2.336249669436741E+02
-1 -2.246750810670323E+02 -5.989076747432608E+01
-1 -7.692139318221516E+01  2.382060858907951E+02
-1 -3.337560467370262E+02  1.762123650236804E+02
-1  1.876340919287669E+02  3.374407372277079E+01
-1  2.242855098885634E+02  2.222744996705513E+02
-1 -5.545625349147973E+01 -1.018187552559464E+02
-1 -2.197207935554015E+02 -2.914546903103951E+02
-1 -2.656264525089730E+02 -2.353259801332774E+02
-1 -1.457174998193881E+02  1.306350685124747E+02
```

### Spatial analysis
![Ripley’s K function](https://raw.githubusercontent.com/DunstanBecht/lpa-input/fc3d0ca5296e33fdcc5f51d7f4c19a581c2d6714/tests/analyses/40000_rho5e13m-2_RRDD-E_s0200nm_f2_square_2000nm_S0_KKKK_GBB.svg)
![Pair correlation function](https://raw.githubusercontent.com/DunstanBecht/lpa-input/fc3d0ca5296e33fdcc5f51d7f4c19a581c2d6714/tests/analyses/40000_rho5e13m-2_RRDD-E_s0200nm_f2_square_2000nm_S0_gggg_GBB.svg)
![Symmetric and antisymmetric functions](https://raw.githubusercontent.com/DunstanBecht/lpa-input/fc3d0ca5296e33fdcc5f51d7f4c19a581c2d6714/tests/analyses/40000_rho5e13m-2_RRDD-E_s0200nm_f2_square_2000nm_S0_GaGs_GBB.svg)

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

The directory `tests/` contains several examples of package module usage. The docstrings are carefully written and it is recommended to refer to the documentation with the `help()` python command.

The installation from PyPI does not allow the modification of the code. To edit the package and contribute to the development use the following commands in your working directory.
```bash
pip uninstall lpa-input
git clone https://github.com/DunstanBecht/lpa-input.git
pip install -e lpa-input
cd lpa-input
git branch <name_of_your_new_branch>
```
