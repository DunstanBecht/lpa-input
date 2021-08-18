#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module models.
"""

from lpa.input import models
import numpy as np

prm_rdd = {'d': 5e13*1e-18}
prm_rrdd = {'v': 'E', 's': 200, 'f': 2}
prm_rcdd = {'v': 'R', 'd': 5e13*1e-18, 's': 500, 't': 20}

if __name__ == "__main__":

    G = np.random.default_rng(0)

    # generate ticks
    print(models.ticks('circle', 100, 50))
    print(models.ticks('square', 200, 50))
    print()

    # generate evenly distributed positions and Burgers vector senses
    t = models.ticks('square', 1, 0.9)[:-1]
    x, y = models.even_positions(t, 2)
    b = models.even_senses (t, 2)
    print(np.column_stack((b, x, y)))
    print()

    # use RDD model
    p, b = models.RDD('circle', 100, 2*np.pi*100**2, prm_rdd, G)
    print(np.column_stack((b, p)))
    print()

    # use RRDD model
    p, b = models.RRDD('square', 400, 400**2, prm_rrdd, G)
    print(np.column_stack((b, p)))
    print()

    # use RCDD model
    p, b = models.RCDD('square', 400, 400**2, prm_rcdd, G)
    print(np.column_stack((b, p)))
    print()

    input("OK")
