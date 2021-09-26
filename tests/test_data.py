#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module data.
"""

from lpa.input import data
from test_sets import *
import numpy as np

d = sets.Distribution('circle', 400, *rdd, S=0)
s = sets.Sample(10, 'square', 400, *rrdd, S=0)

# export the input data file of a distribution or a sample
data.export(d, expfmt='txt')
data.export(s, expdir='data', g=np.array([1,0,0]), expstm='stem')
data.export(s, expdir='data', g=np.array([1,0,0]), pbc=1)

input("OK")
