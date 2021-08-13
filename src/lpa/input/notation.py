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
    # user-defined name for the set of parameters
    if 'name' in r:
        return r['name']
    n = ""
    # writing of equalities
    if t == 'file':
        weq = lambda a, b: a+b
    elif t == 'console':
        weq = lambda a, b: a+"="+b
    elif t == 'title':
        weq = lambda a, b: a+r" = "+b
    # version
    if 'v' in r:
        n += "-"+r['v']
    # parameter list
    p = []
    if 'd' in r:
        p.append(weq('d', number(r['d'], t))) # density
    if 's' in r:
        p.append(weq("s", number(r['s'], t, 4))) # subarea or cell side
    if 'f' in r:
        p.append(weq("f", str(r['f']))) # filling
    if 't' in r:
        p.append(weq("t", number(r['t'], t, 3))) # wall thickness
    if 'l' in r:
        p.append(weq("l", number(r['l'], t, 3))) # dipole length
    if 'r' in r:
        p.append(weq("r", str(r['r'])))
    # concatenate
    if len(p) > 0:
        if t == 'title':
            n += r" $ \left( "+", \ ".join(p)+r" \right) $ "
        elif t == 'console':
            n += " ("+" ".join(p)+")"
        else:
            n += "_"+"_".join(p)
    return n
