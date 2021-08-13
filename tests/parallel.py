#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module overlap.

This script should be submitted with the command "sbatch parallel.job".
"""

from lpa.input import parallel
from lpa.input import sets
from lpa.input.models import RRDD

r = {'v': 'E', 'f': 2, 's': 200, 'r': 0}
s = sets.Sample(10000, 'circle', 1000, RRDD, r)

# export the pooled statistical analysis of the samples s of each core
parallel.export(s)
