#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module overlap.
"""

from lpa.input import overlap
import numpy as np
import matplotlib.pyplot as plt

"""
The following lines display the relative overlapping area between two
circles in three different cases:
- Two circles of radius 1 whose centers are spaced by a distance 2.
- Two circles of radius 1 at the same position.
- A circle of radius 1000 and a circle of radius 1 centered on its edge.
"""
print("Circle-circle overlapping")
rA = np.array([1, 1, 1000])
rB = 1
d = np.array([2, 0, 1000])
area = overlap.circle_circle(rA, rB, d, rA**2, rB**2, d**2)
print(area/np.pi/rB**2)
print()

"""
The following lines display the relative overlapping area between two
a circle of radius 1 and a square of side 2 in three different cases:
- The circle is centered at (1, 1).
- The circle is centered at (2, 1).
- The circle is centered at (2, 2).
"""
print("Circle-square overlapping")
x = np.array([1, 2, 2])
y = np.array([1, 1, 2])
r = 1
s = 2
area = overlap.circle_square(x, y, r, r**2, s)
print(area/np.pi/r**2)
print()

"""
In the following we plot the mean circle-circle overlapping function
and the mean circle-square overlapping function.
"""

rng = np.random.default_rng(0)
x = np.sort(rng.random(10000))

print("Mean circle-circle overlapping")
x_cc = 2*x
cc1 = overlap.mean_circle_circle_analytic(x_cc, 1)
cc2 = overlap.mean_circle_circle_interpolation(x_cc, 1)
print(np.std(cc1-cc2)/np.mean(cc1))
print()

print("Mean circle-square overlapping")
x_cs = np.sqrt(2)*x
cs1 = overlap.mean_circle_square_analytic(x_cs, 1)
cs2 = overlap.mean_circle_square_interpolation(x_cs, 1)
print(np.std(cs1-cs2)/np.mean(cs1))
print()

plt.plot(x_cc, cc1, label="circle-circle")
plt.plot(x_cc, cc2, ':', label="circle-circle")
plt.plot(x_cc, cs1, label="circle-square")
plt.plot(x_cc, cs2, ':', label="circle-square")
plt.show()

input("OK")
