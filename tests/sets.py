#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module sets.
"""

from lpa.input import sets
from lpa.input.models import RDD

# instantiate
d = sets.Distribution('circle', 1000, RDD, {'d': 1e15*1e-18}, c='IDBC')
s = sets.Sample(10, 'square', 2000, RDD, {'d': 1e15*1e-18}, t='edge')

# get file name proposition
print(d.stem(t=False))
print(s.stem(s=False), end="\n\n")

# get representation
print(repr(d))
print(repr(s), end="\n\n")

# evaluate a representation
eval('sets.'+repr(d))
eval('sets.'+repr(s))

# get plot title proposition
print(d.title())
print(s.title(), end="\n\n")

# print the object
print(d, end="\n\n")
print(str(s), end="\n\n")

# average over a sample
rho = s.average(lambda dist: len(dist)/dist.v)
print(str(rho))
