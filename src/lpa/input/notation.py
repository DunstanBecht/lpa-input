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
    """
    return format(round(x), '1.0e').replace("+", "").replace("0", "")


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
    #if 'name' in r:
    #    return r['name']
    n = ""
    if 'v' in r:
        n += r['v']
    if 's' in r:
        n += "_s"+number(r['s'])+"nm"
    if 'f' in r:
        n += "_f"+str(r['f'])
    if 'l' in r:
        n += "_l"+number(r['l'])+"nm"
    return n
