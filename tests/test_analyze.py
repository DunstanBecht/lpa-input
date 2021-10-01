#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.input import analyze
from test_sets import *

d = sets.Distribution('circle', 1000, *rdd, c='ISD', S=0)
s = sets.Sample(10, 'square', 2000, *rrdd, S=0)
r = np.linspace(0, 2000, 100)

# calculate statistical functions on a distribution or a sample
KKKK, GaGs = analyze.calculate(['KKKK', 'GaGs'], d, r, ec='NEC')
gggg, MMMM = analyze.calculate(['gggg', 'MMMM'], s, r)

# export a complete statistical analysis of a distribution or a sample
analyze.export(d, expdir='analyses/', expstm='stem', savtxt=True)
analyze.export(s, expdir='analyses', figttl='title', intrad=r)
analyze.export(s, expdir='analyses', figttl='title', edgcon='WOA')
analyze.export(s, expdir='analyses', figttl='title', edgcon='PBC')
analyze.export(s, expdir='analyses', figttl='title', edgcon='GBB')

input("OK")
