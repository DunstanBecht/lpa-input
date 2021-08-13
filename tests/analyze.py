#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.input import analyze
from lpa.input import sets
from lpa.input.models import RDD

d = sets.Distribution('circle', 1000, RDD, {'d': 1e15*1e-18})
s = sets.Sample(2, 'circle', 1000, RDD, {'d': 1e15*1e-18})

# create a study range
r, iK = analyze.intervals(100, 1000)

# calculate statistical functions on a distribution or a sample
KKKK, GaGs = analyze.calculate(['KKKK', 'GaGs'], d, r)
gggg, KKKK = analyze.calculate(['gggg', 'KKKK'], s, r)

# export a complete statistical analysis of a distribution or a sample
analyze.export(d, p='analyses')
analyze.export(s, p='analyses')
