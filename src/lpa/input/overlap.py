#!/usr/bin/env python
# coding: utf-8

"""
Tools for calculating the overlapping of two geometric objects.
"""

import scipy.integrate
from . import *

@beartype
def circular_segment(
    r: Union[Scalar, ScalarList],
    d: Union[Scalar, ScalarList],
) -> Union[Scalar, ScalarList]:
    """
    Return the circular segment area.

    Input:
        r (Scalar|ScalarList): circle radius/ii
        d (Scalar|ScalarList): distance(s) of the chord from the center

    Output:
        o (Scalar|ScalarList): circular segment area
    """
    return r**2*np.arccos(d/r) - d*np.sqrt(r**2-d**2)

@beartype
def circle_circle(
    rA: Union[Scalar, ScalarList],
    rB: Union[Scalar, ScalarList],
    d: Union[Scalar, ScalarList],
    r2A: Union[Scalar, ScalarList],
    r2B: Union[Scalar, ScalarList],
    d2: Union[Scalar, ScalarList],
) -> Union[Scalar, ScalarList]:
    """
    Return the overlapping area of two circles.

    All input parameters can be either an array or a scalar. If one of
    them is an array, the result will be an array of the same size.

    Input:
        rA (Scalar|ScalarList): radius/ii of the circle A
        rB (Scalar|ScalarList): radius/ii of the circle B
        d (Scalar|ScalarList): distance(s) between the circle centers
        r2A (Scalar|ScalarList): squared radius/ii of the circle A
        r2B (Scalar|ScalarList): squared radius/ii of the circle B
        d2 (Scalar|ScalarList): squared distance(s) between the centers

    Output:
        o (Scalar|ScalarList): overlapping area/s

    Complexity:
        O( max(rA.size, rB.size, d.size) )
    """
    m0 = rA + rB <= d # mask: zero intersection
    mA = d + rA <= rB # mask: A is inside B
    mB = d + rB <= rA # mask: B is inside A
    m = np.logical_not(m0) & np.logical_not(mA) & np.logical_not(mB)
    # m is true for non-trivial cases
    o = np.subtract(
        np.add(
            np.multiply(
                r2A,
                np.arccos(
                    np.divide(
                        np.add(d2, np.subtract(r2A, r2B, where=m), where=m),
                        np.multiply(np.multiply(2, d, where=m), rA, where=m),
                        where=m),
                    where=m),
                where=m),
            np.multiply(
                r2B,
                np.arccos(
                    np.divide(
                        np.add(d2, np.subtract(r2B, r2A, where=m), where=m),
                        np.multiply(np.multiply(2, d, where=m), rB, where=m),
                        where=m),
                    where=m),
                where=m),
            where=m),
        np.divide(
            np.sqrt(
                np.multiply(
                    np.multiply(
                        np.subtract(np.add(rA, rB, where=m), d, where=m),
                        np.add(np.add(rA, rB, where=m), d, where=m),
                        where=m),
                    np.multiply(
                        np.add(np.subtract(rA, rB, where=m), d, where=m),
                        np.add(np.subtract(rB, rA, where=m), d, where=m),
                        where=m),
                    where=m),
                where=m),
            2,
            where=m),
        where=m)
    o = np.where(mA, np.pi*r2A, o)
    o = np.where(mB, np.pi*r2B, o)
    o = np.where(m0, 0, o)
    return o

@beartype
def circle_square(
    x: Union[Scalar, ScalarList],
    y: Union[Scalar, ScalarList],
    r: Union[Scalar, ScalarList],
    r2: Union[Scalar, ScalarList],
    s: Union[Scalar, ScalarList],
) -> Union[Scalar, ScalarList]:
    """
    Return the overlapping area of a circle and a square.

    All input parameters can be either an array or a scalar. If one of
    them is an array, the result will be an array of the same size.

    Input:
        x (Scalar|ScalarList): x coordinate(s) of the circle center
        y (Scalar|ScalarList): y coordinate(s) of the circle center
        r (Scalar|ScalarList): circle radius/ii
        r2 (Scalar|ScalarList): squared circle radius/ii
        s (Scalar|ScalarList): square side(s)

    Output:
        o (Scalar|ScalarList): overlapping area/s

    Complexity:
        O( max(x.size, r.size, s.size) )
    """
    e = [] # list of areas of the circle outside the square for each frame
    dA = s - x # width of the right quadrant
    dB = s - y # height of the upper quadrant
    dC = x # width of the left quadrant
    dD = y # height of the lower quadrant
    # loop on the four quadrants and add their contributions to E
    for d1, d2 in [[dA, dB], [dB, dC], [dC, dD], [dD, dA]]:
        mq = d1**2 + d2**2 > r2
        m1 = mq & (d1 < r)
        m2 = mq & (d2 < r)
        mr = np.logical_not(mq)
        # area of the circle outside the square in the current quadrant
        eq = np.where(mr, np.subtract(np.pi*r2/4, d1*d2, where=mr), 0)
        for a in [[m1, d1], [m2, d2]]:
            np.add(
                np.divide(
                    np.subtract(
                        np.multiply(
                            r2,
                            np.arccos(
                                np.divide(a[1], r, where=a[0]),
                                where=a[0]),
                            where=a[0]),
                        np.multiply(
                            a[1],
                            np.sqrt(
                                np.subtract(
                                    r2,
                                    np.multiply(a[1], a[1], where=a[0]),
                                    where=a[0]),
                                where=a[0]),
                            where=a[0]),
                        where=a[0]),
                    2,
                    where=a[0]),
                eq,
                eq,
                where=a[0])
        e.append(eq)
    return np.pi*r2 - sum(e)

@np.vectorize
@beartype
def mean_circle_circle_analytic(
    r: Scalar,
    R: Scalar,
) -> Scalar:
    """
    Return the mean overlapping area of two circles.

    All input parameters can be either an array or a scalar. If one of
    them is an array, the result will be an array of the same size.

    Input:
        r (Scalar|ScalarList): circle 1 radius/ii
        R (Scalar|ScalarList): circle 2 radius/ii

    Output:
        o (Scalar|ScalarList): mean overlapping area/s

    Complexity:
        O( r.size )
    """
    r2 = r**2
    R2 = R**2
    f = lambda x: circle_circle(r, R, R*np.sqrt(x), r2, R2, R2*x)
    return scipy.integrate.quad(f, 0, 1)[0]

@np.vectorize
@beartype
def mean_circle_square_analytic(
    r: Scalar,
    s: Scalar,
) -> Scalar:
    """
    Return the mean overlapping area of a circle and a square.

    All input parameters can be either an array or a scalar. If one of
    them is an array, the result will be an array of the same size.

    Input:
        r (Scalar|ScalarList): circle radius/ii
        s (Scalar|ScalarList): square side(s)

    Output:
        o (Scalar|ScalarList): mean overlapping area/s

    Complexity:
        O( r.size )
    """
    f12 = lambda x, phi: circular_segment(r, x*np.cos(phi))*x/(2*s**2)
    f3a = lambda x, phi: (np.pi*r**2/4 - x*np.cos(phi)*x*np.sin(phi))*x/(s**2)
    f3b = lambda x, y: (np.pi*r**2/4 - x*y)/s**2
    phi1 = np.arccos(s/max(r, s))
    phi2 = np.arctan(s/min(r, s))
    phi3 = np.arcsin(s/max(r, s))
    x1 = lambda phi: min(r, s)/np.cos(phi)
    x2 = lambda phi: s/np.sin(phi)
    x3 = lambda y: np.tan(phi1)*y
    if r == 0:
        return 0
    elif r <= np.sqrt(2)*s:
        E12 = (scipy.integrate.dblquad(f12, phi1, phi2, r, x1)[0]
            + scipy.integrate.dblquad(f12, phi2, phi3, r, x2)[0])
        E3 = (scipy.integrate.dblquad(f3a, phi1, phi3, 0, r)[0]
            + 2*scipy.integrate.dblquad(f3b, 0, s, 0, x3)[0])
        return np.pi*r**2 - 4*(2*E12+E3)
    else:
        return s**2

@beartype
def mean_circle_circle_simulation(
    r: Union[Scalar, ScalarList],
    R: Union[Scalar, ScalarList],
    n: int = 1000000,
    G: np.random._generator.Generator = np.random.default_rng(0),
) -> Union[Scalar, ScalarList]:
    """
    Return the mean overlapping area of two circles.

    Input parameters r and R can be either an array or a scalar. If one
    of them is an array, the result will be an array of the same size.

    Input:
        r (Scalar|ScalarList): circle 1 radius/ii
        R (Scalar|ScalarList): circle 2 radius/ii
        n (int): number of tested positions
        G (np.random._generator.Generator): random number generator

    Output:
        o (Scalar|ScalarList): mean overlapping area/s

    Complexity:
        O( r.size )
    """
    d = R*np.sqrt(G.random(n))
    d2 = d**2
    m = lambda r, R: circle_circle(r, R, d, r**2, R**2, d2).mean()
    return np.vectorize(m)(r, R)

@beartype
def mean_circle_square_simulation(
    r: Union[Scalar, ScalarList],
    s: Union[Scalar, ScalarList],
    n: int = 1000000,
    G: np.random._generator.Generator = np.random.default_rng(0),
) -> Union[Scalar, ScalarList]:
    """
    Return the mean overlapping area of a circle and a square.

    All input parameters can be either an array or a scalar. If one of
    them is an array, the result will be an array of the same size.

    Input:
        r (Scalar|ScalarList): circle radius/ii
        s (Scalar|ScalarList): square side(s)
        n (int): number of tested positions
        G (np.random._generator.Generator): random number generator

    Output:
        o (Scalar|ScalarList): mean overlapping area/s

    Complexity:
        O( r.size )
    """
    x = G.random(n)*s
    y = G.random(n)*s
    m = lambda ri, si: circle_square(x, y, ri, ri**2, si).mean()
    return np.vectorize(m)(r, s)
