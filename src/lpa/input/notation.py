#!/usr/bin/env python
# coding: utf-8

"""
Tools to standardize notations.
"""

from . import *

@beartype
def number(
    x: Scalar,
    t: str = 'file',
    w: int = 4,
) -> str:
    """
    Return a notation of x that can be used in a file name.

    Input:
        x: number
        t: type of use ('file', 'console', 'title')
        w: width of the output string

    Output:
        s: string of x
    """
    if x<10**w and x>=1 or x == 0: # diplay x with fixed comma
        s = str(round(x))
        if t == 'file':
            s = s.zfill(w) # fill with zeros
        elif t == 'console':
            s = format(s, '>'+str(w)) # fill with spaces
    elif t == 'title': # LaTeX scientific notation
        s = format(x, '1.'+str(max(w-6, 0))+'e')
        m, e = s.split("e") # separate mantissa and exponent
        s = m+" \\times 10^{"+str(int(e))+"}"
    else: # scientific notation with 'e'
        s = format(x, '1.0e').replace("+", "")
        s = s.replace("e0", "e").replace("e-0", "e-")
        d = w-len(s) # number of missing characters to reach the ideal width
        if d>0:
            if t=='file' or d==1:
                if "e-" in s:
                    s = s.replace("e-", "e-0")
                else:
                    s = s.replace("e", "e0")
            else:
                fmt = '1.'+str(d-1)+'e'
                s = format(x, fmt)
                s = s.replace("+", "").replace("e0", "e").replace("e-0", "e-")
    return s

@beartype
def parameters(
    r: dict,
    t: str = 'file',
) -> str:
    """
    Return a string descibing the model parameters.

    Input:
        r: model parameters
        t: type of use ('file', 'console', 'title')

    Output:
        n: model parameters description
    """
    if 'name' in r: # return the name chosen for the parameter set
        return r['name']
    if t == 'file':
        sep = "_" # separator
        weq = lambda a, b: a+b # writing of equalities
    elif t == 'console':
        sep = ", " # separator
        weq = lambda a, b: a+"="+b # writing of equalities
    elif t == 'title':
        sep = ", \ " #separator
        weq = lambda a, b: a+r" = "+b # writing of equalities
    n = ""
    if 'v' in r:
        n += "-"+r['v'] # model version
    if t == 'title':
        n += r" $ \left("
    elif t == 'console':
        n += " ("
    else:
        n += sep
    if 'd' in r:
        n += weq('d', number(r['d'], t)) # density
    if 's' in r:
        n += sep+weq("s", number(r['s'], t, 4)) # subarea or cell side
    if 'f' in r:
        n += sep+weq("f", str(r['f'])) # filling
    if 't' in r:
        n += sep+weq("t", number(r['t'], t, 3)) # wall thickness
    if 'l' in r:
        n += sep+weq("l", number(r['l'], t, 3)) # dipole length
    if t == 'title':
        n += r"\right) $ "
    elif t == 'console':
        n += ")"
    return n
