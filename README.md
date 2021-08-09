<div align="center">
  <img width="250" src="https://dunstan.becht.network/views/signatures/mines.svg" alt="Mines Saint-Etienne">
</div>

# Project

This repository is related to the analysis of crystals containing dislocations by X-ray diffraction. It is part of a project conducted during a research internship at the laboratory of material and structural sciences of the *École Nationale Supérieure des Mines de Saint-Étienne*.

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

# User guide

### Installation

The project is indexed on [PyPI](https://pypi.org/project/lpa-input/) and installable directly via pip.
```bash
pip install lpa-input
```

### Generation
To create a uniformly random dislocation distribution with evenly distributed Burgers vectors in a cylindrical geometry with a radius of 1000 nm:
```python
from lpa.input import sets
r = {'density': 0.03, 'variant': 'e'}
d = sets.Distribution('circle', 1000, 'urdd', r)
```

To create a sample of 100 uniformly random dislocation distribution with evenly distributed Burgers vectors in a cylindrical geometry with a radius of 1000 nm:
```python
from lpa.input import sets
r = {'density': 0.03, 'variant': 'e'}
s = sets.Sample(500, 'circle', 1000, 'urdd', r)
```

### Export

To export a dislocation map of a distribution `d`.
```python
from lpa.input import maps
maps.export(d)
```

To make standardized files for input to an X-ray diffraction simulation program from a sample `s`:
```python
from lpa.input import data
data.export(s)
```

To make a spatial analysis of a sample `s`:
```python
from lpa.input import analyze
analyze.export(s)
```
