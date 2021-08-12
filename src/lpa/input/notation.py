#!/usr/bin/env python
# coding: utf-8

"""
Tools to standardize notations.
"""

from . import *

@beartype
def number(
    x: Scalar,
) -> str:
    """
    Return a notation of x that can be used in a file name.

    Input:
        x: number

    Output:
        s: string of x
    """
    s = format(round(x), '1.0e').replace("+", "")
    return s.replace("0", "").strip("e")


@beartype
def parameters(
    r: dict,
) -> str:
    """
    Return a string descibing the model parameters.

    Input:
        r: model parameters

    Output:
        n: model parameters description
    """
    if 'name' in r:
        return r['name']
    n = ""
    if 'v' in r:
        n += "-"+r['v']
    if 's' in r:
        n += "_s"+number(r['s'])+"nm"
    if 'f' in r:
        n += "_f"+str(r['f'])
    if 't' in r:
        n += "_t"+number(r['t'])+"nm"
    if 'l' in r:
        n += "_l"+number(r['l'])+"nm"
    return n
