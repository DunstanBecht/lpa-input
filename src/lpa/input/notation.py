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
) -> str:
    """
    Return a notation of x that can be used in a file name.

    Input:
        x: number
        t: type of use ('file', 'console', 'title')

    Output:
        s: string of x
    """
    if x<10000 and x>1:
        s = str(round(x))
        if t == 'file':
            s = s.zfill(4)
        elif t == 'console':
            s = format(s, '>4')
        return s
    if t == 'file':
        s = format(x, '1.0e').replace("+", "").replace("0", "")
    elif t == 'console':
        s = format(x, '1.2e').replace("+", "")
    elif t == 'title':
        e = format(x, '1.0e').split("e") # mantissa and exponent
        s = e[0]+"\\times 10^{"+str(int(e[1]))+"}" # LaTeX
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
    if 'name' in r:
        return r['name']
    if t == 'file':
        sep = "_"
        eq = lambda a, b: a+b
    elif t == 'console':
        sep = " "
        eq = lambda a, b: a+"="+b
    elif t == 'title':
        sep = ", \ "
        eq = lambda a, b: a+r" = "+b
    n = ""
    if 'v' in r:
        n += "-"+r['v']
    if len(r) != 0:
        if t == 'title':
            n += r" $ \left("
        else:
            n += sep
        if 'd' in r:
            n += eq('d', number(r['d'], t))
        if 's' in r:
            n += sep+eq("s", number(r['s'], t))
        if 'f' in r:
            n += sep+eq("f", str(r['f']))
        if 't' in r:
            n += sep+eq("t", number(r['t'], t))
        if 'l' in r:
            n += sep+eq("l", number(r['l'], t))
        if t == 'title':
            n += r"\right) $ "
    return n
