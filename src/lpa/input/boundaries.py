#!/usr/bin/env python
# coding: utf-8

"""
Tools for the generation of special conditions at the boundaries.
"""

from . import *

@beartype
def IDBC(
    s: Scalar,
    p: VectorList,
    b: ScalarList,
    t: str,
) -> Tuple[VectorList, ScalarList]:
    """
    Return the image dislocations of (p, b) for a circle of radius s.

    Adding the image dislocations allows to simulate correctly the
    complete diffraction profile (and not only for the asymptotic
    values of the Fourier variable). However, the formula used here
    only applies to dislocations of type 'screw'.

    Input:
        s: radius of the region of interest [nm]
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]
        t: dislocation type

    Output:
        cp: image dislocation positions [nm]
        cb: image dislocation Burgers vectors sense [1]

    Complexity:
        O( len(p) )
    """
    if t == 'screw':
        n = np.linalg.norm(p, axis=1) # distance to the origin
        m = n != 0 # mask to avoid division by zero
        v = p[m] # points different from (0, 0)
        t = np.arctan2(v[:,1], v[:,0]) # angles of the image dislocations
        r = s**2/n[m] # radii of the image dislocations
        cp = np.stack((r*np.cos(t), r*np.sin(t)), axis=1) # image positions
        cb = -b # image Burgers vector senses
        return cp, -cb
    else:
        raise ValueError("cannot calculate images of edge dislocations")

@beartype
def PBCG(
    s: Scalar,
    p: VectorList,
    b: ScalarList,
    r: int,
) -> Tuple[VectorList, ScalarList]:
    """
    Return the replicated dislocations of (p, b) for a square of side s.

    Input:
        s: side of the region of interest [nm]
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]
        r: number of replications at the boundaries

    Output:
        cp: replicated dislocation positions [nm]
        cb: replicated dislocation Burgers vectors sense [1]

    Complexity:
        O( len(p) * r^2 )
    """
    u = [] # list containing translated positions
    for i in range(1, r+1):
        for j in range(2*i):
            for k in (1, -1):
                u.append(np.array((-i*k, (i-j)*k)))
                u.append(np.array(((i-j)*k,  i*k)))
    cp = np.concatenate([p + s*d for d in u]) # replicated positions
    cb = np.tile(b, len(u)) # replicated Burgers vector senses
    return cp, cb
