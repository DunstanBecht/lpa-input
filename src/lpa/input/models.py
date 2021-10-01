#!/usr/bin/env python
# coding: utf-8

"""
Models for the random generation of dislocations.

The following fields can be used in the model parameters:
    v: variant (in RRDD or RCDD)
    d: density (in RDD, RRDD or RCDD)
    s: side of the subarea (in RRDD) or the cell (for RCDD)
    f: filling of the subareas (in RRDD or RCDD)
    t: wall thickness (in RCDD)
    l: dipole length (in RCDD)
    name: parameters set nickname (in RDD, RRDD or RCDD)
"""

import math
from . import *
from . import geometries

@beartype
def RDD(
    g: str,
    s: Scalar,
    a: Scalar,
    r: dict,
    G: np.random._generator.Generator,
) -> tuple:
    """
    Return the positions and Burgers vector generated with RDD model.

    "RDD" means: "Random Dislocation Distribution"

    The positions of the dislocations are uniformly distributed in the
    region of interest.

    The number of dislocations with a positive Burger vector is equal
    to the number of dislocations with a negative Burger vector.

    Input:
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]
        a (Scalar): area of region of interest [nm^2]
        r (dict): model parameters
        G (np.random._generator.Generator): random number generator

    Output:
        p (VectorList): dislocation positions [nm]
        b (ScalarList): dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'd' (Scalar): density of dislocations [nm^-2]

    Complexity:
        O( r['d'] * a )
    """
    nh = round(r['d']*a/2) # half of the total number of dislocations
    nt = 2*nh # total number of dislocations
    # dislocation Burgers vectors
    b = np.concatenate((np.ones(nh), -np.ones(nh)))
    # dislocation positions
    if g == 'circle':
        phi = 2*np.pi*G.random(nt)
        rad = s*np.sqrt(G.random(nt))
        p = np.stack((rad*np.cos(phi), rad*np.sin(phi)), axis=1)
    elif g == 'square':
        p = s*G.random([nt, 2])
    return p, b

@beartype
def ticks(
    g: str,
    l: Scalar,
    s: Scalar,
) -> ScalarList:
    """
    Standardize the construction of a grid of subareas or cells.

    Input:
        g (str): geometry
        l (Scalar): shape size [nm]
        s (Scalar): size of the step [nm]

    Output:
        t (ScalarList): ticks along an axis

    Input example:
        g = 'circle'
        l = 3
        s = 1

    Output example:
        t = np.array([-3, -2, -1, 0, 1, 2, 3])

    Complexity:
        O( l/s )
    """
    m = math.ceil(l/s)
    if np.abs(m*s-l)/l>0.0001:
        warnings.warn("the step is not a divisor of the length", Warning)
    if g == 'circle':
        return s*np.arange(-m, m+1) # [-2s, -s, 0, s, 2s] with m=2
    if g == 'square':
        return s*np.arange(m+1) # [0, s, 2s] with m=2

@beartype
def even_positions(
    t: ScalarList,
    f: int,
) -> ScalarListList:
    """
    Return the grid case lower left corner positions.

    Input:
        t (ScalarList): ticks along an axis
        f (int): number of points in each case

    Output:
        p (ScalarList): case lower left corner positions

    Input example:
        t = np.array([0, 1])
        f = 4

    Output example:
        p = np.array([
            [0, 0], [0, 0], [0, 0], [0, 0],
            [1, 0], [1, 0], [1, 0], [1, 0],
            [0, 1], [0, 1], [0, 1], [0, 1],
            [1, 1], [1, 1], [1, 1], [1, 1],
        ])

    """
    x = np.repeat(np.tile(t, len(t)), f)
    y = np.repeat(t, f*len(t))
    return np.stack((x,y), axis=1)

@beartype
def even_senses(
    t: ScalarList,
    f: int,
) -> ScalarList:
    """
    Return the evenly distributed Burgers senses for even_positions.

    Input:
        t (ScalarList): ticks along an axis
        f (int): number of points in each case

    Output:
        b (ScalarList): Burgers vector senses

    Input example:
        t = np.array([0, 1])
        f = 4

    Output example:
        b = [- - + + - - + + - - + + - - + +]
    """
    bp = np.ones(f//2, dtype=int)
    b = np.tile(np.concatenate((bp, -bp)), len(t)**2)
    return b

@beartype
def RRDD(
    g: str,
    s: Scalar,
    a: Scalar,
    r: dict,
    G: np.random._generator.Generator,
) -> tuple:
    """
    Return the positions and Burgers vector generated with RRDD model.

    "RRDD" means: "Restrictedly Random Dislocation Distribution"

    There is a fixed number of dislocations per subarea. These
    dislocations are uniformly distributed within the subareas.

    In the square geometry, the number of dislocations with a positive
    Burger vector is equal to the number of dislocations with a
    negative Burger vector. When the parameter 'v' is set to 'R', the
    Burgers vectors are randomly distributed. When it is set to 'E',
    the number of positive dislocations is equal to the number of
    negative dislocations in each subarea.

    Input:
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]
        a (Scalar): area of region of interest [nm^2]
        r (dict): model parameters
        G (np.random._generator.Generator): random number generator

    Output:
        p (VectorList): dislocation positions [nm]
        b (ScalarList): dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'v' (str): variant ('R' or 'E')
        's' (int): side of the subareas [nm]
        'f' (int): number of dislocations per subarea
        'd' (Scalar): density of dislocations [nm^-2]

    The parameter 'd' is considered only if the parameter 'f' is not
    specified.

    Complexity:
        O( r['f'] * s^2/r['s']^2 )
    """
    if 'f' in r:
        f = r['f'] # specified filling
    else:
        f = round(r['d']*r['s']**2) # deduced filling
    if f%2 != 0:
        raise ValueError("odd number of dislocations per subarea")
    t = ticks(g, s, r['s'])[:-1] # left or bottom location of subareas
    nt = f * len(t)**2 # total number of dislocations
    nh = nt//2 # half of the total number of dislocations
    # dislocation positions
    p = even_positions(t, f) + r['s']*G.random((nt, 2))
    # dislocation Burgers vector
    if r['v'] == 'R':
        b = G.permutation(np.concatenate((np.ones(nh), -np.ones(nh))))
    elif r['v'] == 'E':
        b = even_senses(t, f)
    # masking
    if g == 'circle':
        m = geometries.mask(g, s, p)
        p, b = p[m], b[m]
    return p, b

@beartype
def RCDD(
    g: str,
    s: Scalar,
    a: Scalar,
    r: dict,
    G: np.random._generator.Generator,
) -> tuple:
    """
    Return the positions and Burgers vector generated with RCDD model.

    "RCDD" means: "Random Cell Dislocation Distribution"

    There is a fixed number of dislocations per subarea. These
    dislocations are uniformly distributed within the cell walls.

    In the square geometry, the number of dislocations with a positive
    Burger vector is equal to the number of dislocations with a
    negative Burger vector. When the parameter 'v' is set to 'R', the
    Burgers vectors are randomly distributed. When it is set to 'E',
    the number of positive dislocations is equal to the number of
    negative dislocations in each cell. When it is set to 'D', the
    positive dislocations and negative dislocations are spaced by a
    distance defined with the parameter 'l'.

    Input:
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]
        a (Scalar): area of region of interest [nm^2]
        r (dict): model parameters
        G (np.random._generator.Generator): random number generator

    Output:
        p (VectorList): dislocation positions [nm]
        b (ScalarList): dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'v' (str): variant ('R', 'E' or 'D')
        's' (int): side of the subareas [nm]
        't' (int): thickness of the cell walls [nm]
        'l' (int): length of the dipoles [nm]
        'f' (int): number of dislocations per subarea
        'd' (Scalar): density of dislocations [nm^-2]

    Complexity:
        O( r['d'] * a )
    """
    if 'f' in r:
        f = r['f'] # specified filling
    else:
        f = round(r['d']*r['s']**2) # deduced filling
    f = 2*max(1, round(f/2))
    if r['t'] > r['s']/2:
        raise ValueError("inconsistent thickness and cell side")
    t = ticks(g, s, r['s'])[:-1] # left or bottom location of subareas
    nt = f * len(t)**2 # total number of dislocations
    nh = nt//2 # half of the total number of dislocations
    # dislocation positions
    if r['v'] == 'D':
        rp = G.random((nh, 2))
    else:
        rp = G.random((nt, 2))
    l1 = r['t']/2 # wall brick length 1
    l2 = r['s'] - l1 # wall brick length 2
    m = G.integers(4, size=len(rp))
    m2 = np.stack((m, m)).T
    rp = np.where(m2==0, np.array([l1,  0])+np.array([l2, l1])*rp, rp)
    rp = np.where(m2==1, np.array([ 0, l2])+np.array([l2, l1])*rp, rp)
    rp = np.where(m2==2, np.array([ 0,  0])+np.array([l1, l2])*rp, rp)
    rp = np.where(m2==3, np.array([l2, l1])+np.array([l1, l2])*rp, rp)
    if r['v'] == 'D':
        phi = 2*np.pi*G.random(nh)
        ux = r['l']/2 * np.cos(phi)
        uy = r['l']/2 * np.sin(phi)
        ux = np.stack((ux, -ux), axis=1).ravel()
        uy = np.stack((uy, -uy), axis=1).ravel()
        rp = np.repeat(rp, 2, axis=0) + np.stack((ux, uy), axis=1)
    p = even_positions(t, f) + rp
    # dislocation Burgers vector
    if r['v'] == 'R':
        b = G.permutation(np.concatenate((np.ones(nh), -np.ones(nh))))
    elif r['v'] == 'E':
        b = even_senses(t, f)
    elif r['v'] == 'D':
        b = np.tile((1, -1), nh)
    # masking
    if g == 'circle':
        m = geometries.mask(g, s, p)
        p, b = p[m], b[m]
    return p, b
