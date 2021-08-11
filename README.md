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
* **rdd**: random dislocation distribution
* **rrdd**: restrictedly random dislocation distribution
* **rcdd**: random cell dislocation distribution

### Model variants
* **r**: randomly distributed Burgers vectors
* **e**: evenly distributed Burgers vectors
* **d**: dipolar Burgers vectors

### Boundary conditions
* **pbcg**: periodic boundary conditions applied when generating the distribution
* **pbcr**: periodic boundary conditions applied when running the simulation
* **idbc**: image dislocations boundary conditions

# User guide

### Installation

The project is indexed on [PyPI](https://pypi.org/project/lpa-input/) and installable directly via pip.
```bash
pip install lpa-input
```

### Generation
To create a random dislocation distribution with evenly distributed Burgers vectors in a cylindrical geometry with a radius of 1000 nm:
```python
from lpa.input import sets
r = {'density': 0.03, 'variant': 'e'}
d = sets.Distribution('circle', 1000, 'rdd', r)
```

To create a sample of 100 random dislocation distributions with evenly distributed Burgers vectors in a cylindrical geometry with a radius of 1000 nm:
```python
from lpa.input import sets
r = {'density': 0.03, 'variant': 'e'}
s = sets.Sample(500, 'circle', 1000, 'rdd', r)
```

### Exportation

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

### Parallelization

To parallelize the spatial analysis of distributions on a supercomputer equipped with Slurm Workload Manager, it is necessary to create two files. The first one is the python script that will be executed on each core.

`slurm.py`
```python
#!/usr/bin/env python
# coding: utf-8

"""
This script is executed on each core during a parallel analysis.
"""

import time
from lpa.input import sets
from lpa.input import parallel
import settings

n = 1 # number of distribution per core

p = [
    [n, *settings.circle, *settings.rrdde13],
    [n, *settings.circle, *settings.rrdde14],
]

if parallel.rank == parallel.root:
    t1 = time.time()
for args in p:
    s = sets.Sample(*args)
    if parallel.rank == parallel.root:
        print("- analysis of "+s.fileName()+" ", end="")
        t2 = time.time()
    parallel.export(s)
    if parallel.rank == parallel.root:
        print("("+str(round((time.time()-t2)/60))+" mn)")
if parallel.rank == parallel.root:
    print("total time: "+str(round((time.time()-t1)/60))+" mn")
```

The second file is used to submit the task to Slurm.

`slurm.job`
```bash
#!/bin/bash
#SBATCH --job-name=disldist
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=10:00:00
#SBATCH --partition=intensive.q
ulimit -l unlimited
###unset SLURM_GTIDS

SCRIPT=slurm.py
echo ------------------------------------------------------
echo number of nodes in the job resource allocation: $SLURM_NNODES
echo nodes allocated to the job: $SLURM_JOB_NODELIST
echo directory from which sbatch was invoked: $SLURM_SUBMIT_DIR
echo hostname of the computer from which sbatch was invoked: $SLURM_SUBMIT_HOST
echo id of the job allocation: $SLURM_JOB_ID
echo name of the job: $SLURM_JOB_NAME
echo name of the partition in which the job is running: $SLURM_JOB_PARTITION
echo number of nodes requested: $SLURM_JOB_NUM_NODES
echo number of tasks requested per node: $SLURM_NTASKS_PER_NODE
echo ------------------------------------------------------
echo generating hostname list
COMPUTEHOSTLIST=$( scontrol show hostnames $SLURM_JOB_NODELIST |paste -d, -s )
echo ------------------------------------------------------
echo creating scratch directories on nodes $SLURM_JOB_NODELIST
SCRATCH=/scratch/$USER-$SLURM_JOB_ID
srun -n$SLURM_NNODES mkdir -m 770 -p $SCRATCH || exit $?
echo ------------------------------------------------------
echo transferring files from frontend to compute nodes $SLURM_JOB_NODELIST
srun -n$SLURM_NNODES cp -rvf $SLURM_SUBMIT_DIR/$SCRIPT $SCRATCH || exit $?
echo ------------------------------------------------------
echo load packages
module load anaconda/python3
python3 -m pip install -U lpa-input
echo ------------------------------------------------------
echo run -mpi program
cd $SCRATCH
mpirun --version
mpirun -np $SLURM_NTASKS -npernode $SLURM_NTASKS_PER_NODE -host $COMPUTEHOSTLIST python3 $SLURM_SUBMIT_DIR/$SCRIPT
echo ------------------------------------------------------
echo transferring result files from compute nodes to frontend
srun -n$SLURM_NNODES cp -rvf $SCRATCH $SLURM_SUBMIT_DIR || exit $?
echo ------------------------------------------------------
echo deleting scratch from nodes $SLURM_JOB_NODELIST
srun -n$SLURM_NNODES rm -rvf $SCRATCH || exit 0
echo ------------------------------------------------------
```

Finally, to start the simulation enter the following command.

```bash
sbatch slurm.job
```
