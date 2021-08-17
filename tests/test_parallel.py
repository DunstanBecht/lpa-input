#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module overlap.

This script should be submitted with the command:
sbatch test_parallel.job
"""

from lpa.input import parallel
from test_sets import *

# export the pooled statistical analysis of the samples s of each core
parallel.export(d_rrdd, exfmt='svg')
