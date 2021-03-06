#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module maps.
"""

from lpa.input import maps
from test_sets import *

"""
The following lines instantiate some dislocation distributions.
"""
d1 = sets.Distribution('circle', 1000, *rdd, S=0)
d2 = sets.Distribution('square', 2000, *rrdd, S=0)
d3 = sets.Distribution('circle', 2000, *rcdd, S=0, c='ISD')
d4 = sets.Distribution('square', 2000, *rcdd, S=0, c='PBC1')
d5 = sets.Distribution('square', 2000, *rrdd, S=0, c='GBB1')

"""
The following lines export the dislocation maps of the previsously
instantiated distributions.
"""
maps.export(d1, expfmt='svg')
maps.export(d2, expfmt='svg')
maps.export(d3, expfmt='svg')
maps.export(d4, expfmt='svg')
maps.export(d5, expdir='maps/', supttl='suptitle', subttl='subtitle')

input("OK")
