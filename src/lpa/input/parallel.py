#!/usr/bin/env python
# coding: utf-8

"""
Tools for parallelizing the spatial analysis of distributions.
"""

from mpi4py import MPI
name = MPI.Get_processor_name()
comm = MPI.COMM_WORLD
rank = comm.Get_rank() # number of the processor executing this script
size = comm.Get_size() # number of processors
root = 0 # master processor
import warnings
from . import *
from . import sets
from . import analyze

@beartype
def average_on_cores(
    w: AnalysisOutput,
    b: bool = False
) -> AnalysisOutput:
    """
    Return the average value of w ​​over the cores.

    When b is True, all cores get the average value. When b is False,
    the workers get None and the root gets the average value.

    Input:
        w (AnalysisOutput): worker value
        b (bool): broadcast the result to all cores

    Output:
        m (AnalysisOutput): averaged value of w over the cores

    Input example:
        On master core:
            w = 1
            b = False
        On worker core:
            w = 0
            b = False

    Output example:
        On master core:
            m = 0.5
        On worker core:
            m = None
    """
    if not isinstance(w, np.ndarray):
        w = np.array(w)
    if rank == root: # if the script is executed by the master
        m = np.zeros_like(w) # prepare the buffer for data reception
    else: # if the script is executed by a worker
        m = None
    comm.Reduce([w, MPI.DOUBLE], [m, MPI.DOUBLE], op=MPI.SUM, root=root)
    if rank == root:
        m = m/size # average the value
    if b:
        m = comm.bcast(m, root=root) # broadcast to the workers
    return m

@beartype
def export(
    o: Union[sets.Distribution, sets.Sample],
    **kwargs,
) -> None:
    """
    Export a complete pooled analysis of the object o of each core.

    Function similar to the export function of the analyze module.

    Input:
        o (Distribution|Sample): distribution or sample to analyze on the core
      **expdir (str): export directory
      **expfmt (str): export format
      **expstm (str): export stem
      **title (str): figure title
    """
    if rank==root and not o.S is None:
        msg = ("chosen random seed detected on "+o.name('stm')+" "
            + "(The use of a seed is not recommended. In order for "
            + "the distributions or samples not to be identical from "
            + "one core to another, the seed must not be identical "
            + "from one core to another.")
        warnings.warn(msg, Warning)
    i = average_on_cores(o.i, True) # averaged inter dislocation distance
    r, iK = analyze.intervals(i, o.s) # intervals to display
    f = ['KKKK', 'gggg', 'GaGs'] # functions to calculate
    worker = analyze.calculate(f, o, r) # function results
    master = [average_on_cores(worker[i]) for i in range(len(worker))]
    if rank == root:
        if isinstance(o, sets.Distribution):
            c = str(size) # number of distributions analyzed
        else:
            c = str(len(o)*size) # number of distributions analyzed
        # optional parameters
        expdir = getkwa('expdir', kwargs, str, '')
        expfmt = getkwa('expfmt', kwargs, str, 'pdf')
        expstm = getkwa('expstm', kwargs, str, c+"_"+o.name('dmgsS', c='stm'))
        title = getkwa('title', kwargs, str, c+" "+o.name('mgsd', c='ttl'))
        endkwa(kwargs)
        # export
        args = (expdir, expfmt, expstm, title)
        analyze.plot_KKKK(r[:iK], master[f.index('KKKK')].T[:iK].T, *args)
        analyze.plot_gggg(r, master[f.index('gggg')], *args)
        analyze.plot_GaGs(r, master[f.index('GaGs')], *args)
