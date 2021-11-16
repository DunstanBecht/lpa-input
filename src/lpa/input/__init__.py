#!/usr/bin/env python
# coding: utf-8

"""
Tools for the generation and analysis of dislocation distributions.

A dislocation is characterized by a position and a Burgers vector. A
distribution is a set of dislocations in a material. The dislocations
of a distribution are generated randomly according to a model. A sample
is a set of distributions generated with the same model.

A distribution can be exported as a map showing the location of its
dislocations. It can also be exported in a standardized form for input
to an X-ray diffraction simulation program. Finally, it is possible to
make a spatial analysis of a model by means of statistical functions
on a sample of distributions.

In this package, Burgers vectors are represented only by their sense
(positive: +1, or negative: -1). Their direction is only determined in
the header of the exported input data files. It is parallel to the
dislocation line for the 'screw' type or perpendicular for the 'edge'
type. The superposition principle simplifies the analysis by allowing
to study only one direction at a time.

To obtain relevant statistical results, the spatial analysis must be
averaged over a large number of distributions. To reduce the
computation time it is possible to parallelize the calculations. The
MPI library can be used for this purpose. A sample of distributions is
generated and analyzed on each core, then the results are averaged
over the cores.
"""

__author__ = "Dunstan Becht"
__version__ = "1.2.4"

import os
import sys
import warnings
from typing import Union, Optional, Any

if sys.version_info[0]>=3 and sys.version_info[1]>=9:
    from collections.abc import Callable
else:
    from typing import Callable

if sys.version_info[0]>=3 and sys.version_info[1]>=8:
    from typing import get_args, get_origin
else:
    get_args = lambda x: x.__args__
    get_origin = lambda x: x.__origin__ if hasattr(x, '__origin__') else x

import numpy as np
from beartype import beartype

# scalar and vectors
Scalar = Union[int, np.integer, float, np.floating]
Vector = np.ndarray # shape: (n,)

# sets
ScalarList = np.ndarray # shape: (...,)
ScalarListList = np.ndarray # shape: (..., ...)
VectorList = np.ndarray # shape: (..., n)

# generic analysis function output
AnalysisOutput = Union[tuple, Scalar, np.ndarray]

# edge correction functions
CorrectionFunction = Callable[[Vector, ScalarList, ScalarList], ScalarList]

# model generation functions
GenerationFunction = Callable[[str, Scalar, Scalar, dict], tuple]

def getkwa(
    key: str,
    kwa: dict,
    typ: type,
    dft: Any,
) -> Any:
    """
    Retrieve the keyword argument and check its type.

    Input:
        key (str): keyword argument name
        kwa (dict): keyword arguments dictionary
        typ (type): keyword argument type
        dft (Any): keyword argument default value

    Output:
        val (Any): value of the keyword argument
    """
    val = kwa.pop(key, dft) # get the argument or its default value
    if get_origin(typ) == Union: # if the type is a union of types
        typ = get_args(typ) # get the accepted types
    if not isinstance(val, typ): # if the type is not correct
        msg = (f"keyword argument {key!r} must be of type {typ.__name__!r} "
               f"but {val!r} is of type {type(val).__name__!r}")
        raise TypeError(msg)
    return val

def endkwa(
    kwa: dict,
) -> None:
    """
    Check that no bad keyword arguments have been passed.

    Input:
        kwa (dict): keyword arguments dictionary
    """
    if len(kwa) > 0: # if there are still keyword arguments after recovery
        msg = "unexpected keyword argument: {!r}={!r}".format(*kwa.popitem())
        raise TypeError(msg)
