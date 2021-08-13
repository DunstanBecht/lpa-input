#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module models.
"""

from lpa.input import models
import numpy as np

# choose or not a random seed
r1 = {'r': 0}
r2 = {}
models.seed(r1)
models.seed(r2)
print(r1)
print(r2, end="\n\n")

# generate ticks
s = 100
t = models.ticks('square', s, 50)
print(t, end="\n\n")

# generate evenly distributes positions and Burgers vector senses
x, y = models.even_positions(t, 2)
b = models.even_senses (t, 2)
print(np.concatenate((b, x, y)), end="\n\n")

# use RDD model
p, b = models.RDD('square', s, s**2, {'d': 0.0001})
print(np.column_stack((b, p)), end="\n\n")

# use RRDD model
p, b = models.RRDD('square', s, s**2, {'v': 'E', 'f': 2, 's': 100})
print(np.column_stack((b, p)), end="\n\n")

# use RCDD model
p, b = models.RCDD('square', s, s**2, {'v': 'R', 'd': 0.0005, 's': 50, 't': 10})
print(np.column_stack((b, p)), end="\n\n")
