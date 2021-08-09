#!/usr/bin/env python
# coding: utf-8

"""
Models for the random generation of dislocations.
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
    if not 'seed' in r:
        r['seed'] = np.random.randint(0, 2**32, dtype=np.int64)
    np.random.seed(r['seed'])

@beartype
def urdd(
    g: str,
    s: Scalar,
    a: Scalar,
    r: dict,
) -> Tuple[VectorList, VectorList]:
    """
    Return the positions and Burgers vector generated with urdd model.

    "urdd" means: "uniformly random dislocation distribution"

    Input:
        g: geometry of the region of interest
        s: size of the region of interest [nm]
        a: area of region of interest [nm^2]
        r: model parameters

    Output:
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]

    The following parameters can be specified in r:
        'seed' (int): random seed
        'variant' (str): 'r' or 'e'
        'density' (Scalar): density of dislocations [nm^-2]

    When 'variant' equals 'r', the Burgers vectors are randomly
    distributed. When 'variant' equals 'e', the Burgers vectors are
    evenly distributed. When 'seed' is not specified, the random seed
    is chosen randomly.

    Complexity:
        O( r['density'] * a )
    """
    seed(r)
    # dislocation positions
    n = round(r['density']*a) # number of dislocations
    if g == 'circle':
        theta = 2*np.pi*np.random.random(n)
        radius = s*np.sqrt(np.random.random(n))
        p = np.stack((radius*np.cos(theta), radius*np.sin(theta)), axis=1)
    elif g == 'square':
        p = s*np.random.random([n, 2])
    # dislocation Burgers vectors
    if r['variant'] == 'r':
        b = np.random.choice([1, -1], size=n)
    elif r['variant'] == 'e':
        c = n//2
        l = [np.ones(c), -np.ones(c)]
        if n%2 == 1:
            l.append(np.random.choice([1, -1], size=1))
        b = np.concatenate(l)
    return p, b

@beartype
def ticks(
    g: str,
    l: int,
    s: int,
) -> ScalarList:
    """
    Standardize the construction of a grid of subareas.

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
        'seed': (int): random seed
        'variant' (str): 'r' or 'e'
        'subarea' (int): size of the subareas [nm]
        'filling' (int): number of dislocations per subarea

    When 'variant' equals 'r', the Burgers vectors are randomly
    distributed. When 'variant' equals 'e', the Burgers vectors are
    evenly distributed. When 'seed' is not specified, the random seed
    is chosen randomly.

    Complexity:
        O( r['filling'] * s^2/r['subarea']^2 )
    """
    seed(r)
    # division of the space in subareas
    t = ticks(g, s, r['subarea'])[:-1] # left or bottom location of subareas
    n = r['filling'] * len(t)**2 # number of dislocations
    x = np.repeat(np.tile(t, len(t)), r['filling'])
    # x=[0,0,0,1,1,1,0,0,0,1,1,1] with t=[0,1], f=3
    y = np.repeat(t, r['filling']*len(t))
    # y=[0,0,0,0,0,0,1,1,1,1,1,1] with t=[0,1], f=3
    p = np.stack((x,y), axis=1) + r['subarea']*np.random.random((n, 2))
    # dislocation Burgers vector
    if r['variant'] == 'r':
        b = np.random.choice([1, -1], size=n)
    elif r['variant'] == 'e':
        c = r['filling']//2
        bp, bm, l = np.ones(c), -np.ones(c), []
        for i in range(len(t)**2):
            ls = [bp, bm]
            if r['filling']%2 == 1:
                ls.append(np.random.choice([1, -1], size=1))
            l.append(np.concatenate(ls))
        b = np.concatenate(l)
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
        'seed': (int): random seed
        'variant' (str): 'r', 'e' or 'd'
        'density' (Scalar): density of dislocations [nm^-2]
        'subarea' (int): size of the subareas [nm]
        'border' (int): thickness of the border [nm]
        'dipole' (int): size of dislocation dipoles [nm]

    When 'variant' equals 'r', the Burgers vectors are randomly
    distributed. When 'variant' equals 'e', the Burgers vectors are
    evenly distributed. When 'variant' equals 'd' the Burgers vectors
    with sense + and sense - are spaced by a distance 'dipole'. When
    'seed' is not specified, the random seed is chosen randomly.

    Complexity:
        O( r['density'] * a )
    """
    seed(r)
    if r['variant']=='d' and r['dipole']>r['border']:
        raise ValueError('inconsistent sizes')
    # positions
    t = ticks(g, s, r['subarea'])[:-1] # left or bottom location of subareas
    n = round(r['density']*r['subarea']**2*len(t)**2) # number of dislocations
    if r['variant'] == 'd':
        c = n//2 # 1 random point for 1 dislocation
    else:
        c = n # 1 random point for 1 dipole of 2 dislocations
    if r['variant'] == 'd':
        l1 = (r['border']-r['dipole'])/2
    else:
        l1 = r['border']/2
    l2 = r['subarea'] - l1
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
    if r['variant'] == 'r':
        b = np.random.choice([1, -1], size=c)
    elif r['variant'] == 'e':
        l = [np.ones(c//2), -np.ones(c//2)]
        if c%2 == 1:
            l.append(np.random.choice([1, -1], size=1))
        b = np.concatenate(l)
    elif r['variant'] == 'd':
        b = np.concatenate((np.ones(c), -np.ones(c)))
        theta = 2*np.pi*np.random.random(c)
        ux = r['dipole']/2 * np.cos(theta)
        uy = r['dipole']/2 * np.sin(theta)
        x = np.concatenate((x+ux, x-ux))
        y = np.concatenate((y+uy, y-uy))
    p = np.stack((x,y), axis=1)
    if g == 'circle':
        mask = np.sum(np.square(p), axis=1) < s**2
    elif g == 'square':
        mask = (x<=s) & (x>=0) & (y<=s) & (y>=0)
    p = p[mask]
    b = b[mask]
    return p, b
