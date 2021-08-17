#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module data.
"""

from lpa.input import data
from test_sets import *
import numpy as np

# export the input data file of a distribution or a sample
data.export(d_rdd, exfmt='txt')
data.export(s_rrdd, exdir='data', g=np.array([1,0,0]))
data.export(s_rcdd, exdir='data', exstm='stem')

input("OK")
