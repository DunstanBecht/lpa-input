#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module data.
"""

from lpa.input import data
from test_sets import *
import numpy as np

d = sets.Distribution('circle', 400, *rdd, S=0)
s = sets.Sample(10, 'circle', 400, *rrdd, S=0)

# export the input data file of a distribution or a sample
data.export(d, exfmt='txt')
data.export(s, exdir='data', g=np.array([1,0,0]), exstm='stem')

input("OK")
