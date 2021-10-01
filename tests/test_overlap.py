#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module overlap.
"""

from lpa.input import overlap
import numpy as np

print("Circle-circle overlapping")
rA = np.array([1, 1, 1000])
rB = 1
d = np.array([2, 0, 1000])
area = overlap.circle_circle(rA, rB, d, rA**2, rB**2, d**2)
print(area/np.pi/rB**2)
print()

print("Circle-square overlapping")
x = np.array([1, 2, 2])
y = np.array([1, 1, 2])
r = 1
s = 2
area = overlap.circle_square(x, y, r, r**2, s)
print(area/np.pi/r**2)
print()

input("OK")
