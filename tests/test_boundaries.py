#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module boundaries.
"""

from lpa.input import boundaries
import numpy as np

print("Image screw dislocation positions")
p = np.array([[1, 0], [0, 0]])
print(boundaries.image_positions(2, p))
print()

print("Displacements for replications")
print(boundaries.replication_displacements(1, 2))
print()

input("OK")
