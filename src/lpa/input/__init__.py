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
to an X-ray diffraction simulation program. Finally it is possible to
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
__version__ = "0.9.1"

import sys
import numpy as np
from typing import Union, Optional, NewType, Any
from beartype import beartype

if sys.version_info[0]>=3 and sys.version_info[1]>=9:
    from collections.abc import Callable
    List = list
    Tuple = tuple
else:
    from typing import Callable
    from typing import List
    from typing import Tuple

Scalar = Union[int, float, np.intc]
Vector = np.ndarray # shape: (n,)
ScalarList = np.ndarray # shape: (...,)
ScalarListList = np.ndarray # shape: (..., ...)
VectorList = np.ndarray # shape: (..., n)
CorrectionFunction = Callable[[Vector, ScalarList, ScalarList], ScalarList]
GenerationFunction = Callable[[str, int], Tuple[VectorList, VectorList]]
AnalysisOutput = Union[Tuple, Union[Scalar, np.ndarray]]
