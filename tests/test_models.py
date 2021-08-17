#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module models.
"""

from lpa.input import models
import numpy as np

G = np.random.default_rng(0)

# generate ticks
s = 100
t = models.ticks('square', s, 50)
print(t, end="\n\n")

# generate evenly distributes positions and Burgers vector senses
x, y = models.even_positions(t, 2)
b = models.even_senses (t, 2)
print(np.concatenate((b, x, y)), end="\n\n")

# use RDD model
r = {'d': 0.0001}
p, b = models.RDD('circle', s, 2*np.pi*s**2, r, G)
print(np.column_stack((b, p)), end="\n\n")

# use RRDD model
r = {'v': 'E', 'f': 2, 's': 100}
p, b = models.RRDD('square', s, s**2, r, G)
print(np.column_stack((b, p)), end="\n\n")

# use RCDD model
r = {'v': 'R', 'd': 0.0005, 's': 50, 't': 10}
p, b = models.RCDD('square', s, s**2, r, G)
print(np.column_stack((b, p)), end="\n\n")

input("OK")
