#!/usr/bin/env python
# coding: utf-8

"""
Tools for managing the possible geometries of the region of interest.
"""

from . import *

@beartype
def nvolume(
    g: str,
    s: Scalar
) -> tuple:
    """
    Return the dimension of space n and the n-volume of the shape.

    Input:
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]

    Output:
        n (int): dimension of space
        v (Scalar): n-volume of the shape [nm^n]
    """
    if g == 'circle':
        return 2, np.pi*s**2
    elif g == 'square':
        return 2, s**2
    else:
        raise ValueError("unknown geometry: "+str(g))

@beartype
def mask(
    g: str,
    s: Scalar,
    p: VectorList,
) -> np.ndarray:
    """
    Return the mask filtering positions inside the region of interest.

    Input:
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]
        p (VectorList): positions of the dislocations [nm]

    Output:
        m (np.ndarray): mask of positions within the region of interest
    """
    if g == 'circle':
        return np.sum(np.square(p), axis=1) < s**2
    elif g == 'square':
        return (p[:,0]<s) & (p[:,1]<s) & (p[:,0]>0) & (p[:,1]>0)
    else:
        raise ValueError("unknown geometry: "+str(g))
