#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module data.
"""

from lpa.input import data
from lpa.input import sets
from lpa.input.models import RDD

d = sets.Distribution('circle', 1000, RDD, {'d': 1e15*1e-18}, c='IDBC')
s = sets.Sample(2, 'circle', 1000, RDD, {'d': 1e15*1e-18})

# export the input data file of a distribution or a sample
data.export(d, p='data')
data.export(s, p='data')
