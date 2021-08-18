#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module overlap.

This script should be submitted with the command:
sbatch test_parallel.job
"""

from lpa.input import parallel
from test_sets import *

d = sets.Distribution('circle', 1000, *rdd, c='IDBC', S=0)
s = sets.Sample(625, 'square', 1000, *rrdd, S=0)

# export the pooled statistical analysis of the samples s of each core
parallel.export(d)
parallel.export(s, exfmt='svg')
