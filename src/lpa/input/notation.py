#!/usr/bin/env python
# coding: utf-8

"""
Tools to standardize notations.

The following abbreviations define the types of use of strings:
    stm: the string is used in the stem of a file name
    ttl: the string is a title and contains LaTeX code
    csl: the string is displayed in a console
"""

from . import *

@beartype
def number(
    x: Scalar,
    u: str = 'csl',
    w: int = 4,
) -> str:
    """
    Return a notation of x for the type of use u.

    Input:
        x: number
        u: type of use ('stm', 'ttl', 'csl')
        w: width of the output string

    Output:
        s: notation of x
    """
    if x<10**w and x>=1 or x == 0: # diplay x with fixed comma
        s = str(round(x))
        if u == 'stm':
            s = s.zfill(w) # fill with zeros
        elif u == 'csl':
            s = format(s, '>'+str(w)) # fill with spaces
    elif u == 'ttl': # LaTeX scientific notation
        s = format(x, '1.'+str(max(w-6, 0))+'e')
        m, e = s.split("e") # separate mantissa and exponent
        s = m+" \\times 10^{"+str(int(e))+"}"
    elif u in ('stm', 'csl'): # scientific notation with 'e'
        s = format(x, '1.0e').replace("+", "")
        s = s.replace("e0", "e").replace("e-0", "e-")
        d = w-len(s) # number of missing characters to reach the ideal width
        if d>0:
            if u=='stm' or d==1:
                if "e-" in s:
                    s = s.replace("e-", "e-0")
                else:
                    s = s.replace("e", "e0")
            else:
                fmt = '1.'+str(d-1)+'e'
                s = format(x, fmt)
                s = s.replace("+", "").replace("e0", "e").replace("e-0", "e-")
    else:
        raise ValueError("unknown usage type: "+str(u))
    return s

@beartype
def parameters(
    r: dict,
    u: str = 'csl',
    s: bool = True,
) -> str:
    """
    Return the notation descibing the model parameters set r.

    Input:
        r: model parameters
        u: type of use ('stm', 'ttl', 'csl')
        s: display the seed

    Output:
        n: model parameters description
    """
    # user-defined name for the set of parameters
    if 'name' in r:
        return r['name']
    n = ""
    # version
    if 'v' in r:
        n += "-"+r['v']
    # parameter list
    p = []
    if 'd' in r:
        p.append(('d', number(r['d'], u))) # density
    if 's' in r:
        p.append(("s", number(r['s'], u, 4))) # subarea or cell side
    if 'f' in r:
        p.append(("f", str(r['f']))) # filling
    if 't' in r:
        p.append(("t", number(r['t'], u, 3))) # wall thickness
    if 'l' in r:
        p.append(("l", number(r['l'], u, 3))) # dipole length
    if 'r' in r and s:
        p.append(("r", str(r['r']))) # random seed
    # concatenate
    if len(p) > 0:
        if u == 'stm':
            d = [a+b for a, b in p]
            n += "_"+"_".join(d)
        elif u == 'ttl':
            d = [a+" = "+b for a, b in p]
            n += r" $ \left( "+", \ ".join(d)+r" \right) $ "
        elif u == 'csl':
            d = [a+"="+b for a, b in p]
            n += " ("+" ".join(d)+")"
        else:
            raise ValueError("unknown usage type: "+str(u))
    return n
