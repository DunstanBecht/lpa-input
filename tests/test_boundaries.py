#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module boundaries.
"""

from lpa.input import boundaries
import numpy as np

# generate image dislocations
p = np.array([[1, 0], [0, 0]])
b = np.array([1])
cp, cb = boundaries.IDBC(2, p, b, 'screw')
print(np.column_stack((cb, cp)), end="\n\n")

# generate replicated dislocations
p = np.array([[10, 10]])
b = np.array([1])
cp, cb = boundaries.PBCG(10, p, b, 1)
print(np.column_stack((cb, cp)), end="\n\n")

input("OK")
