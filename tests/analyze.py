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
rd, iKd = analyze.intervals(d.i, d.s)
rs, iKs = analyze.intervals(s.i, s.s)

# calculate statistical functions on a distribution or a sample
KKKK, GaGs = analyze.calculate(['KKKK', 'GaGs'], d, rd)
gggg, MMMM = analyze.calculate(['gggg', 'MMMM'], s, rs)

# export a complete statistical analysis of a distribution or a sample
analyze.export(d, exdir='analyses')
analyze.export(s, exdir='analyses')
