#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module models.
"""

from lpa.input import models
import numpy as np

"""
The following line defines the parameters for an RDD model with a
dislocation density of 5e13 m^-2.
"""
prm_rdd = {'d': 5e13*1e-18}

"""
The following line defines the parameters for an RRDD model with
2 dislocations in every subarea of 200 nm side. The variant E is
chosen.
"""
prm_rrdd = {'v': 'E', 's': 200, 'f': 2}

"""
The following line defines the parameters for an RCDD model with
a dislocation density of 5e13 m^-2. The cells have a side of 200 nm
and a wall thickness of 20 nm. The variant R is chosen.
"""
prm_rcdd = {'v': 'R', 'd': 5e13*1e-18, 's': 200, 't': 20}

if __name__ == "__main__":

    """
    The following line sets a random state.
    """
    G = np.random.default_rng(0)

    """
    The following lines display the list of values allowing to split
    the region of interest into sub-areas or cells.
    """
    print("Ticks")
    print(models.ticks('circle', 100, 50), "circle")
    print(models.ticks('square', 200, 50), "square")
    print()

    """
    In the following several sets of positions and Burgers vectors are
    generated for different models and parameters.
    """

    print("Evenly distributed positions and Burgers vectors")
    t = models.ticks('square', 2, 1)[:-1]
    p = models.even_positions(t, 2)
    b = models.even_senses (t, 2)
    print(np.column_stack((b, p)))
    print()

    print("RDD")
    p, b = models.RDD('circle', 100, 2*np.pi*100**2, prm_rdd, G)
    print(np.column_stack((b, p)))
    print()

    print("RRDD-E")
    p, b = models.RRDD('square', 400, 400**2, prm_rrdd, G)
    print(np.column_stack((b, p)))
    print()

    print("RCDD-R")
    p, b = models.RCDD('square', 400, 400**2, prm_rcdd, G)
    print(np.column_stack((b, p)))
    print()

    input("OK")
