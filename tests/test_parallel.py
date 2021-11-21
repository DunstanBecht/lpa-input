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

d = sets.Distribution('circle', 1000, *rdd)
s = sets.Sample(1250, 'square', 2000, *rrdd, S=0+parallel.rank)
r = np.linspace(0, 3000, 200)

# export the pooled statistical analysis of the samples s of each core
parallel.export(d)
warnings.filterwarnings("ignore") # (read the warning in the module parallel)
t1 = time.time()
parallel.export(s, edgcon='NEC', intrad=r, expfmt='svg')
t2 = time.time()
parallel.export(s, edgcon='WOA', intrad=r, expfmt='svg')
t3 = time.time()
parallel.export(s, edgcon='PBC', intrad=r, expfmt='svg')
t4 = time.time()
parallel.export(s, edgcon='GBB', intrad=r, expfmt='svg')
t5 = time.time()
if parallel.rank == parallel.root:
    print('NEC', t2-t1)
    print('WOA', t3-t2)
    print('PBC', t4-t3)
    print('GBB', t5-t4)
