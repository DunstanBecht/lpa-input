#!/usr/bin/env python
# coding: utf-8

"""
Models for the random generation of dislocations.

The following parameters can be used in the model parameters:
    r: random seed (rdd, rrdd, rcdd)
    v: variant (rrdd, rcdd)
    d: density (rdd, rrdd, rcdd)
    s: side of the subarea (rrdd) or the cell (rcdd)
    f: filling of the subareas (rrdd)
    l: dipole length (rcdd)
    name: parameters set nickname (rdd, rrdd, rcdd)
"""

import math
from . import *

@beartype
def seed(
    r: dict,
) -> None:
    """
    Set the specified seed or randomly select a seed if not specified.

    Letting the user choose a seed allows the replicability of the
    generation of a distribution. This function must be called at the
    start of the model definitions.

    Input:
        r: model parameters

    Complexity:
        O( 1 )
    """
    if not 'r' in r:
        r['r'] = np.random.randint(0, 2**32, dtype=np.int64)
    np.random.seed(r['r'])

@beartype
def rdd(
    g: str,
    s: Scalar,
    a: Scalar,
    r: dict,
) -> Tuple[VectorList, VectorList]:
    """
    Return the positions and Burgers vector generated with rdd model.

    "rdd" means: "random dislocation distribution"

    Input:
        g: geometry of the region of interest
        s: size of the region of interest [nm]
        a: area of region of interest [nm^2]
        r: model parameters

    Output:
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'r' (int): random seed
        'd' (Scalar): density of dislocations [nm^-2]

    When the parameter 'r' is not specified, the random seed is chosen
    randomly and added to the parameters dictionary.

    Complexity:
        O( r['d'] * a )
    """
    seed(r)
    rd = r['d']
    n = round(rd*a) # number of dislocations
    # dislocation positions
    if g == 'circle':
        theta = 2*np.pi*np.random.random(n)
        radius = s*np.sqrt(np.random.random(n))
        p = np.stack((radius*np.cos(theta), radius*np.sin(theta)), axis=1)
    elif g == 'square':
        p = s*np.random.random([n, 2])
    # dislocation Burgers vectors
    bp = np.ones(n//2, dtype=int)
    l = [bp, -bp]
    if n%2 == 1:
        l.append(np.random.choice([1, -1], size=1))
    b = np.concatenate(l)
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
        l: length to subdivide [nm]
        s: size of the step [nm]

    Output:
        t: ticks along an axis

    Complexity:
        O( l/s )
    """
    m = math.ceil(l/s)
    if g == 'circle':
        return s*np.arange(-m, m+1)
    if g == 'square':
        return s*np.arange(m+1)

@beartype
def rrdd(
    g: str,
    s: Scalar,
    a: Scalar,
    r: dict,
) -> Tuple[VectorList, VectorList]:
    """
    Return the positions and Burgers vector generated with rrdd model.

    "rrdd" means: "restrictedly random dislocation distribution"

    Input:
        g: geometry of the region of interest
        s: size of the region of interest [nm]
        a: area region of interest [nm^2]
        r: model parameters

    Output:
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'r': (int): random seed
        'v' (str): variant ('r' or 'e')
        's' (int): size of the subareas [nm]
        'f' (int): number of dislocations per subarea
        'd' (Scalar): density of dislocations [nm^-2]

    When the parameter 'r' is not specified, the random seed is chosen
    randomly and added to the parameters dictionary. The parameter 'd'
    is considered only if the parameter 'f' is not specified.

    When the parameter 'v' is set to 'r', the Burgers vectors are
    randomly distributed. When it is set to 'e', the Burgers vectors
    are evenly distributed in each subarea.

    Complexity:
        O( r['f'] * s^2/r['s']^2 )
    """
    seed(r)
    rv = r['v']
    rs = r['s']
    if 'filling' in r:
        rf = r['f']
    else:
        rf = round(r['d']*rs**2)
    # division of the space in subareas
    t = ticks(g, s, rs)[:-1] # left or bottom location of subareas
    if rf%2 != 0:
        raise ValueError('odd number of dislocations per subarea')
    n = rf * len(t)**2 # number of dislocations
    x = np.repeat(np.tile(t, len(t)), rf)
    # x=[0,0,0,1,1,1,0,0,0,1,1,1] with t=[0,1], f=3
    y = np.repeat(t, rf*len(t))
    # y=[0,0,0,0,0,0,1,1,1,1,1,1] with t=[0,1], f=3
    p = np.stack((x,y), axis=1) + rs*np.random.random((n, 2))
    # dislocation Burgers vector
    if rv == 'r':
        b = np.random.choice([1, -1], size=n)
    elif rv == 'e':
        bp = np.ones(rf//2, dtype=int)
        b = np.tile(np.concatenate((bp, -bp)), len(t)**2)
    if g == 'circle':
        mask = np.sum(np.square(p), axis=1) < s**2
        p = p[mask]
        b = b[mask]
    return p, b

@beartype
def rcdd(
    g: str,
    s: Scalar,
    a: Scalar,
    r: dict,
) -> Tuple[VectorList, VectorList]:
    """
    Return the positions and Burgers vector generated with rcdd model.

    "rcdd" means: "random cell dislocation distribution"

    Input:
        g: geometry of the region of interest
        s: size of the region of interest [nm]
        a: area region of interest [nm^2]
        r: model parameters

    Output:
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'r': (int): random seed
        'v' (str): 'r', 'e' or 'd'
        'd' (Scalar): density of dislocations [nm^-2]
        's' (int): size of the subareas [nm]
        't' (int): thickness of the border [nm]
        'l' (int): length of the dipoles [nm]

    When the parameter 'r' is not specified, the random seed is chosen
    randomly and added to the parameters dictionary.

    When the parameter 'v' is set to 'r', the Burgers vectors are
    randomly distributed. When 'v' is set to 'e', the Burgers vectors
    are evenly distributed in each cell. When 'v' is set to 'd', the
    Burgers vectors with positive sense and negative sense are spaced
    by a distance defined by 'l'.

    Complexity:
        O( r['d'] * a )
    """
    seed(r)
    rv = r['v']
    rd = r['d']
    rs = r['s']
    rt = r['t']
    if rv=='d' and r['l']>rt:
        raise ValueError('inconsistent sizes')
    # positions
    t = ticks(g, s, rs)[:-1] # left or bottom location of subareas
    n = round(rd*rs**2*len(t)**2) # number of dislocations
    if rv == 'd':
        c = n//2 # 1 random point for 1 dislocation
    else:
        c = n # 1 random point for 1 dipole of 2 dislocations
    if rv == 'd':
        l1 = (rt-r['l'])/2
    else:
        l1 = rt/2
    l2 = rs - l1
    x = 1.0*np.random.choice(t, size=c)
    y = 1.0*np.random.choice(t, size=c)
    m = np.random.randint(4, size=c)
    rx = np.random.random(c)
    ry = np.random.random(c)
    np.add(x, l1 + rx*l2, x, where=m==0)
    np.add(y,      ry*l1, y, where=m==0)
    np.add(x,      rx*l2, x, where=m==1)
    np.add(y, l2 + ry*l1, y, where=m==1)
    np.add(x,      rx*l1, x, where=m==2)
    np.add(y,      ry*l2, y, where=m==2)
    np.add(x, l2 + rx*l1, x, where=m==3)
    np.add(y, l1 + ry*l2, y, where=m==3)
    if rv == 'r':
        b = np.random.choice([1, -1], size=c)
    elif rv == 'e':
        l = [np.ones(c//2), -np.ones(c//2)]
        if c%2 == 1:
            l.append(np.random.choice([1, -1], size=1))
        b = np.concatenate(l)
    elif rv == 'd':
        b = np.concatenate((np.ones(c), -np.ones(c)))
        theta = 2*np.pi*np.random.random(c)
        ux = r['l']/2 * np.cos(theta)
        uy = r['l']/2 * np.sin(theta)
        x = np.concatenate((x+ux, x-ux))
        y = np.concatenate((y+uy, y-uy))
    p = np.stack((x,y), axis=1)
    if g == 'circle':
        mask = np.sum(np.square(p), axis=1) < s**2
        p = p[mask]
        b = b[mask]
    return p, b
