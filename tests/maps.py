#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module maps.
"""

from lpa.input import maps
from lpa.input import sets
from lpa.input.models import RDD, RRDD

rdd = (RDD, {'d': 5e13*1e-18, 'r': 0})
rrdde = (RRDD, {'v': 'E', 'f': 2, 's': 200, 'r': 0})

drdd = sets.Distribution('circle', 1000, *rdd)
drrdde = sets.Distribution('square', 2000, *rrdde)

# export distribution maps
maps.export(drdd, exdir='maps', exstm='stem')
maps.export(drrdde, exdir='maps', title='title')
