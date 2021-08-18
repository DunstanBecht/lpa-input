#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module maps.
"""

from lpa.input import maps
from test_sets import *

d1 = sets.Distribution('circle', 1000, *rdd, S=0)
d2 = sets.Distribution('square', 2000, *rrdd, S=0)
d3 = sets.Distribution('circle', 1000, *rrdd, S=0, c='IDBC')

# export distribution maps
maps.export(d1, exfmt='svg')
maps.export(d2, exfmt='svg')
maps.export(d3, exdir='maps/', ttlsp='suptitle', ttlsb='subtitle')

input("OK")
