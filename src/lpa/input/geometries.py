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
) -> Tuple[Scalar, Scalar]:
    """
    Return the dimension of space n and the n-volume of the shape.

    Input:
        g: geometry of the region of interest
        s: size of the region of interest [nm]

    Output:
        n: dimension of space
        v: n-volume of the shape [nm^n]
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
        g: geometry of the region of interest
        s: size of the region of interest [nm]
        p: positions of the dislocations [nm]

    Output:
        m: mask of positions within the region of interest
    """
    if g == 'circle':
        return np.sum(np.square(p), axis=1) < s**2
    elif g == 'square':
        return (p[:,0]<s) & (p[:,1]<s) & (p[:,0]>0) & (p[:,1]>0)
    else:
        raise ValueError("unknown geometry: "+str(g))
