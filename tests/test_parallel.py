#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module parallel.

This script should be submitted with the command:
sbatch test_parallel.job
"""

from lpa.input import parallel
from test_sets import *
import warnings

d = sets.Distribution('square', 2000, *rdd)
s = sets.Sample(625, 'circle', 1000, *rrdd, S=0+parallel.rank)

# export the pooled statistical analysis of the samples s of each core
parallel.export(d)
warnings.filterwarnings("ignore") # (read the warning in the module parallel)
parallel.export(s, expfmt='svg')
