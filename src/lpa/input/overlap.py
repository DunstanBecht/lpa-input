#!/usr/bin/env python
# coding: utf-8

"""
Tools for calculating the overlapping of two geometric objects.
"""

from . import *

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
        d2 (Scalar|ScalarList): squared distance(s) between the circle centers

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
    return np.where(mA, np.pi*r2A, np.where(mB, np.pi*r2B, o))

@beartype
def circle_square(
    x: Union[Scalar, ScalarList],
    y: Union[Scalar, ScalarList],
    r: Union[Scalar, ScalarList],
    r2: Union[Scalar, ScalarList],
    s: Union[Scalar, ScalarList],
) -> Union[float, ScalarList]:
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
