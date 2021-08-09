<div align="center">
  <img width="250" src="https://dunstan.becht.network/views/signatures/mines.svg" alt="Mines Saint-Etienne">
</div>

# Project

This repository is related to the analysis of crystals containing dislocations by X-ray diffraction. It is part of a larger project conducted during a research internship at the laboratory of material and structural sciences of the *École Nationale Supérieure des Mines de Saint-Étienne*.

# Features

The tools developed can be used to:
* generate dislocation distributions according to different models
* export the distributions in standardized files for input to an X-ray diffraction simulation program
* export the distributions in dislocation maps
* export a spatial analysis of the distributions

# Physical aspects

A dislocation associates:
* a Burgers vector
* a position

Two geometries are proposed:
* circle (intersection of a plane with a cylinder) centered in (0, 0)
* square (intersection of a plane with a cuboid) bottom left corner in (0, 0)

A distribution is characterized by the following elements:
* the geometry of the region of interest
* the model used for the random generation of dislocations
* the generated dislocations

# Abbreviations

Some abbreviations are used in the program:

### Models
* **urdd**: uniformly random dislocation distribution
* **rrdd**: restrictedly random dislocation distribution
* **rcdd**: random cell dislocation distribution

### Model variants
* **r**: randomly distributed Burgers vectors
* **e**: evenly distributed Burgers vectors
* **d**: dipolar Burgers vectors

### Boundary conditions
* **pbcg**: periodic boundary conditions applied when generating the distribution
* **pbcr**: periodic boundary conditions applied when runnning the simulation
* **idbc**: image dislocations boundary conditions

# Content of the repository

The **disldist** directory is the python package containing the distribution generation and analysis tools. The directory **exported** contains the exported files. The directory **maths** provides the LaTeX code of the functions implemented in **disldist**. The **remote.py** script simplifies the remote management of calculations on a supercomputer. The script **run.py** is used as hello world and can be executed as is. The script **settings.py** gathers the standard parameters of the generated distributions. The files **slurm.job** and **slurm.py** are to be used only on a supercomputer.

# Examples of use

### Generation
To create a uniformly random dislocation distribution with evenly distributed Burgers vectors in a cylindrical geometry with a radius of 1000 nm:
```python
from disldist import sets
r = {'density': 0.03, 'variant': 'e'}
d = sets.Distribution('circle', 1000, 'urdd', r)
```

To create a sample of 100 uniformly random dislocation distribution with evenly distributed Burgers vectors in a cylindrical geometry with a radius of 1000 nm:
```python
from disldist import sets
r = {'density': 0.03, 'variant': 'e'}
s = sets.Sample(500, 'circle', 1000, 'urdd', r)
```

### Export

To export a dislocation map of a distribution `d`.
```python
from disldist import maps
maps.export(d)
```

To make standardized files for input to an X-ray diffraction simulation program from a sample `s`:
```python
from disldist import data
data.export(s)
```

To make a spatial analysis of a sample `s`:
```python
from disldist import analyze
analyze.export(s)
```

# License

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
