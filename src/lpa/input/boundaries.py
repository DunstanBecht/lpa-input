#!/usr/bin/env python
# coding: utf-8

"""
Tools for the generation of special conditions at the boundaries.
"""

from . import *

@beartype
def image_positions(
    s: Scalar,
    p: VectorList,
) -> VectorList:
    """
    Return the image screw positions for a circle of radius s.

    Adding the image dislocations allows to simulate correctly the
    complete diffraction profile (and not only for the asymptotic
    values of the Fourier variable). However, the formula used here
    only applies to dislocations of type 'screw'.

    Input:
        s (Scalar): radius of the region of interest [nm]
        p (VectorList): dislocation positions [nm]

    Output:
        cp (VectorList): image dislocation positions [nm]

    Input example:
        s = 2
        p = np.array([[1, 0], [-1, 0], [0, 0]])

    Output example:
        cp = np.array([[4, 0], [-4, 0]])

    Complexity:
        O( len(p) )
    """
    n = np.linalg.norm(p, axis=1) # distance to the origin
    m = n != 0 # mask to avoid division by zero
    v = p[m] # points different from (0, 0)
    t = np.arctan2(v[:,1], v[:,0]) # angles of the image dislocations
    r = s**2/n[m] # radii of the image dislocations
    cp = np.stack((r*np.cos(t), r*np.sin(t)), axis=1) # image positions
    return cp

@beartype
def replication_displacements(
    r: int,
    s: Scalar = 1,
) -> VectorList:
    """
    Return the displacements for the replication of a square of side s.

    Input:
        r (int): rank of replications of the region of interest
        s (Scalar): side of the region of interest

    Output:
        u (VectorList): displacements for the replications

    Input example:
        r = 1
        s = 1

    Output example:
        u = np.array([
            [-2  2],
            [ 2  2],
            [ 2 -2],
            [-2 -2],
            [-2  0],
            [ 0  2],
            [ 2  0],
            [ 0 -2],
        ])

    Complexity:
        O( r^2 )
    """
    u = [] # list containing displacements for the replications
    for i in range(1, r+1):
        for j in range(2*i):
            for k in (1, -1):
                u.append(np.array((-i*k, (i-j)*k)))
                u.append(np.array(((i-j)*k,  i*k)))
    return s*np.array(u)
