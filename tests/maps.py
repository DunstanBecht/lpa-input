#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module maps.
"""

from lpa.input import maps
from lpa.input import sets
from lpa.input.models import RDD

d = sets.Distribution('circle', 1000, RDD, {'d': 5e13*1e-18, 'r': 0})

# export a distribution map
maps.export(d, ep='maps')
