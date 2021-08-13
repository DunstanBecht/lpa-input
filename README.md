<div align="center">
  <img width="250" src="https://dunstan.becht.network/views/signatures/mines.svg" alt="Mines Saint-Etienne">
</div>

# Line Profile Analysis - Input

This repository is related to the analysis of crystals containing dislocations by X-ray diffraction. It is part of a project conducted during a research internship at the laboratory of material and structural sciences of the *École Nationale Supérieure des Mines de Saint-Étienne*. Three python packages have been developed to conduct line profile analyses based on simulation results:
* [`lpa.input`](https://github.com/DunstanBecht/lpa-input) (line profile analysis input generator)
* [`lpa.xrd`](https://github.com/DunstanBecht/lpa-xrd) (line profile analysis x-ray diffraction simulation program)
* [`lpa.output`](https://github.com/DunstanBecht/lpa-output) (line profile analysis output analyzer)

# Features

The package `lpa.input` can be used to:
* generate dislocation distributions according to different models
* export the distributions in standardized files for input to an X-ray diffraction simulation program
* export the distributions in dislocation maps
* export a spatial analysis of the distributions

# Installation

The package is indexed on [PyPI](https://pypi.org/project/lpa-input/) and installable directly via pip:
```bash
pip install -U lpa-input
```

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
