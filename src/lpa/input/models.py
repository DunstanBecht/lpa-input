#!/usr/bin/env python
# coding: utf-8

"""
Models for the random generation of dislocations.

The following parameters can be used in the model parameters:
    v: variant (RRDD, RCDD)
    d: density (RDD, RRDD, RCDD)
    s: side of the subarea (RRDD) or the cell (RCDD)
    f: filling of the subareas (RRDD)
    t: wall thickness (RCDD)
    l: dipole length (RCDD)
    name: parameters set nickname (RDD, RRDD, RCDD)
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
) -> Tuple[VectorList, ScalarList]:
    """
    Return the positions and Burgers vector generated with RDD model.

    "RDD" means: "Random Dislocation Distribution"

    Input:
        g: geometry of the region of interest
        s: size of the region of interest [nm]
        a: area of region of interest [nm^2]
        r: model parameters
        G: random number generator

    Output:
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'd' (Scalar): density of dislocations [nm^-2]

    Complexity:
        O( r['d'] * a )
    """
    # parameters
    rd = r['d'] # density of dislocations
    n = round(rd*a) # number of dislocations
    # dislocation positions
    if g == 'circle':
        theta = 2*np.pi*G.random(n)
        radius = s*np.sqrt(G.random(n))
        p = np.stack((radius*np.cos(theta), radius*np.sin(theta)), axis=1)
    elif g == 'square':
        p = s*G.random([n, 2])
    # dislocation Burgers vectors
    bp = np.ones(n//2, dtype=int) # positive Burgers vectors
    b = np.concatenate((bp, -bp)) # add negative Burgers vectors
    if n%2 == 1: # add a random burger vector
        b = np.concatenate((b, G.choice([1, -1], size=1)))
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
        g: geometry
        l: shape size [nm]
        s: size of the step [nm]

    Output:
        t: ticks along an axis

    Complexity:
        O( l/s )
    """
    m = math.ceil(l/s)
    if g == 'circle':
        return s*np.arange(-m, m+1) # [-2s, -s, 0, s, 2s] with m=2
    if g == 'square':
        return s*np.arange(m+1) # [0, s, 2s] with m=2

@beartype
def even_positions(
    t: VectorList,
    f: int,
) -> Tuple[ScalarList, ScalarList]:
    """
    Return the positions with a number of points f in each case.

    Input:
        t: ticks along an axis
        f: number of points in each case

    Output:
        x: x coordinates of points
        y: y coordinates of points

    Input example:
        t = np.array([0, 1])
        f = 4

    Output example:
        x = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]
        y = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1]
    """
    x = np.repeat(np.tile(t, len(t)), f)
    y = np.repeat(t, f*len(t))
    return x, y

@beartype
def even_senses(
    t: VectorList,
    f: int,
) -> ScalarList:
    """
    Return the evenly distributed Burgers senses for even_positions.

    Input:
        t: ticks
        f: number of points in each case

    Output:
        b: Burgers vector senses

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
) -> Tuple[VectorList, ScalarList]:
    """
    Return the positions and Burgers vector generated with RRDD model.

    "RRDD" means: "Restrictedly Random Dislocation Distribution"

    Input:
        g: geometry of the region of interest
        s: size of the region of interest [nm]
        a: area region of interest [nm^2]
        r: model parameters
        G: random number generator

    Output:
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'v' (str): variant ('R' or 'E')
        's' (int): side of the subareas [nm]
        'f' (int): number of dislocations per subarea
        'd' (Scalar): density of dislocations [nm^-2]

    The parameter 'd' is considered only if the parameter 'f' is not
    specified. When the parameter 'v' is set to 'R', the Burgers
    vectors are randomly distributed. When it is set to 'E', the
    Burgers vectors are evenly distributed in each subarea.

    Complexity:
        O( r['f'] * s^2/r['s']^2 )
    """
    # parameters
    rv = r['v'] # model variant
    rs = r['s'] # size of the subareas
    if 'f' in r:
        rf = r['f'] # specified filling
    else:
        rf = round(r['d']*rs**2) # deduced filling
    if rf%2 != 0:
        raise ValueError('odd number of dislocations per subarea')
    t = ticks(g, s, rs)[:-1] # left or bottom location of subareas
    n = rf * len(t)**2 # number of dislocations
    # dislocation positions
    x, y = even_positions(t, rf)
    p = np.stack((x,y), axis=1) + rs*G.random((n, 2))
    # dislocation Burgers vector
    if rv == 'R':
        b = G.choice([1, -1], size=n)
    elif rv == 'E':
        b = even_senses(t, rf)
    # masking
    m = geometries.mask(g, s, p)
    return p[m], b[m]

@beartype
def RCDD(
    g: str,
    s: Scalar,
    a: Scalar,
    r: dict,
    G: np.random._generator.Generator,
) -> Tuple[VectorList, ScalarList]:
    """
    Return the positions and Burgers vector generated with RCDD model.

    "RCDD" means: "Random Cell Dislocation Distribution"

    Input:
        g: geometry of the region of interest
        s: size of the region of interest [nm]
        a: area region of interest [nm^2]
        r: model parameters
        G: random number generator

    Output:
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'v' (str): variant ('R', 'E' or 'D')
        'd' (Scalar): density of dislocations [nm^-2]
        's' (int): side of the subareas [nm]
        't' (int): thickness of the cell walls [nm]
        'l' (int): length of the dipoles [nm]

    When the parameter 'v' is set to 'R', the positions and the Burgers
    vectors senses are randomly distributed. When 'v' is set to 'E',
    each cell contains a fixed number of dislocations and the Burgers
    vectors senses are evenly distributed in each cell. When 'v' is set
    to 'D', the Burgers vectors with positive sense and negative sense
    are spaced by a distance defined with 'l'.

    Complexity:
        O( r['d'] * a )
    """
    # parameters
    rv = r['v'] # variant
    rd = r['d'] # density of dislocations
    rs = r['s'] # side of the cells
    rt = r['t'] # thickness of the cell walls
    if rt > rs/2:
        raise ValueError("inconsistent thickness and cell side")
    if rv == 'D':
        rl = r['l'] # length of the dipoles
    t = ticks(g, s, rs)[:-1] # left or bottom location of subareas
    n = round(rd*rs**2*len(t)**2) # number of dislocations
    # dislocation or dipole positions
    if rv == 'D':
        c = n//2 # 1 random point for 1 dipole of 2 dislocations
    elif rv == 'R':
        c = n # 1 random point for 1 dislocation
    elif rv == 'E':
        rf = round(r['d']*rs**2) # deduced filling
        if rf%2 != 0:
            raise ValueError("odd number of dislocations per cell")
        c = rf * len(t)**2 # 1 random point for 1 dislocation
    # choice of cells
    if rv in ['R', 'D']:
        x = 1.0*G.choice(t, size=c)
        y = 1.0*G.choice(t, size=c)
    elif rv == 'E':
        x, y = even_positions(t, rf)
    # choice of position in the cells
    l1 = rt/2
    l2 = rs - l1
    m = G.integers(4, size=c)
    rx = G.random(c)
    ry = G.random(c)
    np.add(x, l1 + rx*l2, x, where=m==0)
    np.add(y,      ry*l1, y, where=m==0)
    np.add(x,      rx*l2, x, where=m==1)
    np.add(y, l2 + ry*l1, y, where=m==1)
    np.add(x,      rx*l1, x, where=m==2)
    np.add(y,      ry*l2, y, where=m==2)
    np.add(x, l2 + rx*l1, x, where=m==3)
    np.add(y, l1 + ry*l2, y, where=m==3)
    # Burgers vector
    if rv == 'R':
        b = G.choice([1, -1], size=c)
    elif rv == 'E':
        b = even_senses(t, rf)
    elif rv == 'D':
        b = np.concatenate((np.ones(c), -np.ones(c)))
        theta = 2*np.pi*G.random(c)
        ux = rl/2 * np.cos(theta)
        uy = rl/2 * np.sin(theta)
        x = np.concatenate((x+ux, x-ux))
        y = np.concatenate((y+uy, y-uy))
    p = np.stack((x,y), axis=1)
    # masking
    if g == 'circle':
        m = geometries.mask(g, s, p)
        p = p[m]
        b = b[m]
    return p, b
