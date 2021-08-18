#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.input import analyze
from test_sets import *

d = sets.Distribution('circle', 1000, *rdd, S=0)
s = sets.Sample(10, 'square', 2000, *rcdd, c='PBCG1', S=0)

# create a study range
rd, i = analyze.intervals(d.i, d.s)
rs, _ = analyze.intervals(s.i, s.s)

# calculate statistical functions on a distribution or a sample
KKKK, GaGs = analyze.calculate(['KKKK', 'GaGs'], d, rd)
gggg, MMMM = analyze.calculate(['gggg', 'MMMM'], s, rs)

# export a complete statistical analysis of a distribution or a sample
analyze.export(d, exdir='analyses/', exstm='stem')
analyze.export(s, exdir='analyses', title='title')

input("OK")
