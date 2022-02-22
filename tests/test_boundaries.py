#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module boundaries.
"""

from lpa.input import boundaries
import numpy as np

"""
The following lines display the positions of the dislocations resulting
of the ISD boundary conditions computed for 2 dislocation positions in
a circular region of interest of radius 2.
"""
print("Image screw dislocation positions")
p = np.array([[1, 0], [0, 0]])
print(boundaries.image_positions(2, p))
print()

"""
The following lines display the translations to operate on each
dislocation position in order to obtain PBC boundary conditions with
1 replication applied to a square region of interest of side 2.
"""
print("Displacements for replications")
print(boundaries.replication_displacements(1, 2))
print()

input("OK")
