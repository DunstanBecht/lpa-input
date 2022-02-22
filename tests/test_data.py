#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module data.
"""

from lpa.input import data
from test_sets import *
import numpy as np

"""
The following lines instantiate a distribution and a sample.
"""
d = sets.Distribution('circle', 400, *rdd, S=0)
s = sets.Sample(10, 'square', 400, *rrdd, S=0)

"""
The following line exports the dislocation positions and Burgers
vectors of the distribution in a text file.
"""
data.export(d, expfmt='txt')

"""
The following line exports the dislocation positions and Burgers
vectors of the distributions contained in the sample to the folder
data/ and with a diffraction vector direction (100). The name of the
exported folder will be 'stem'.
"""
data.export(s, expdir='data', g=np.array([1,0,0]), expstm='stem')

"""
The following line exports the dislocation positions and Burgers
vectors of the distributions contained in the sample to the folder
data/ and with a diffraction vector direction (100). The XRD simulation
program will replicate the region of interest 1 time.
"""
data.export(s, expdir='data', g=np.array([1,0,0]), pbc=1)

input("OK")
