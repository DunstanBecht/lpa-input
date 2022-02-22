#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module analyze.
"""

from lpa.input import analyze
from test_sets import *

"""
The following line generates a distribution of dislocations. The
dislocations are arranged in a circular region of interest with a
radius of 1000 nm according to the model RDD and with ISD boundary
conditions. A random seed is also set to 0 in order to reproduce the
same random draw at each execution.
"""
d = sets.Distribution('circle', 1000, *rdd, c='ISD', S=0)

"""
The following line generates a sample of 10 distributions of
dislocations. The dislocations are arranged in a square region of
interest with a side of 2000 nm according to the model RRDD and with
ISD boundary conditions. A random seed is also set to 0 in order to
reproduce the same random draw at each execution.
"""
s = sets.Sample(10, 'square', 2000, *rrdd, S=0)

"""
The following line defines the radius interval to study with the
statistical analysis functions.
"""
r = np.linspace(0, 2000, 100)

"""
The following line is used to calculate the Ripley’s K functions
and the antisymmetrical and symmetrical functions on the previously
instantiated distribution. The edge consideration chosen is NEC.
"""
KKKK, GaGs = analyze.calculate(['KKKK', 'GaGs'], d, r, ec='NEC')

"""
The following line is used to calculate the pair correlation functions
and the M functions on the previously instantiated sample of
distributions. The edge consideration chosen by default is NEC.
"""
gggg, MMMM = analyze.calculate(['gggg', 'MMMM'], s, r)

"""
The following line is used to export the plots of the Ripley’s K
functions, the pair correlation functions and the antisymmetrical and
symmetrical functions for the previously instantiated distribution.
The files are exported to the folder analyses/ and will be named from
the character string 'stem'. A backup of the calculated values will
be exported to a text file.
"""
analyze.export(d, expdir='analyses/', expstm='stem', savtxt=True)

"""
The following lines are used to export the plots of the Ripley’s K
functions, the pair correlation functions and the antisymmetrical and
symmetrical functions for the previously instantiated sample of
distributions and with different considerations at the edges. The files
are exported to the folder analyses/.
"""
analyze.export(s, expdir='analyses', figttl='title', intrad=r)
analyze.export(s, expdir='analyses', figttl='title', edgcon='WOA')
analyze.export(s, expdir='analyses', figttl='title', edgcon='PBC')
analyze.export(s, expdir='analyses', figttl='title', edgcon='GBB')

input("OK")
