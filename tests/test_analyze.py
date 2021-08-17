#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.input import analyze
from test_sets import *

# create a study range
rd, i = analyze.intervals(d_rdd.i, d_rdd.s)
rs, _ = analyze.intervals(s_rrdd.i, s_rrdd.s)

# calculate statistical functions on a distribution or a sample
KKKK, GaGs = analyze.calculate(['KKKK', 'GaGs'], d_rdd, rd)
gggg, MMMM = analyze.calculate(['gggg', 'MMMM'], s_rrdd, rs)

# export a complete statistical analysis of a distribution or a sample
analyze.export(d_rrdd, exdir='analyses/', exstm='stem')
analyze.export(s_rdd, exdir='analyses', title='title')

input("OK")
