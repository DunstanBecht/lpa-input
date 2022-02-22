#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module parallel.

This script should be submitted on a supercomputer with the command:
$ sbatch test_parallel.job
"""

from lpa.input import parallel
from test_sets import *
import warnings
import time

"""
The following lines instantiate a distribution and a sample of 1250
distributions on each core. The total number of distributions generated
is obtained by multiplying by the variable parallel.size giving the
number of cores.
"""
d = sets.Distribution('circle', 1000, *rdd)
s = sets.Sample(1250, 'square', 2000, *rrdd, S=0+parallel.rank)

"""
The following line defines the radius interval to study with the
statistical analysis functions.
"""
r = np.linspace(0, 3000, 200)

"""
The following line exports the pooled statistical analysis of
the distributions instantiated on each core.
"""
parallel.export(d)

"""
The following line disables the warning about random seeds.
(Read the warning in the module parallel.)
"""
warnings.filterwarnings("ignore")

"""
The following lines perform a benchmark while analysing the sample.
"""
t1 = time.time()
parallel.export(s, edgcon='NEC', intrad=r, expfmt='svg')
t2 = time.time()
parallel.export(s, edgcon='WOA', intrad=r, expfmt='svg')
t3 = time.time()
parallel.export(s, edgcon='PBC', intrad=r, expfmt='svg')
t4 = time.time()
parallel.export(s, edgcon='GBB', intrad=r, expfmt='svg')
t5 = time.time()

"""
The following lines display the running times on the main core only.
"""
if parallel.rank == parallel.root:
    print('NEC', t2-t1)
    print('WOA', t3-t2)
    print('PBC', t4-t3)
    print('GBB', t5-t4)
