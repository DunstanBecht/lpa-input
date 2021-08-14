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

from . import *
from . import sets
from . import analyze

def average_on_cores(
    w,
    b: bool = False
):
    """
    Return the average value of w ​​over the cores.

    When b is True, all cores get the average value. When b is False,
    the workers get None and the root gets the average value.

    Input:
        w: worker value
        b: broadcast the result to all cores

    Output:
        m: averaged value of w over the cores
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

def export(
    o: Union[sets.Distribution, sets.Sample],
    p: str = "",
) -> None:
    """
    Export a complete pooled analysis of the object o of each core.

    Function similar to the export function of the analyze module.

    Input:
        o: distribution or sample to analyze on the core
        p: path where to export the file
        n: name of the exported file
        t: title of the plot
    """
    if p!="" and p[-1]!="/":
        p += "/"
    i = average_on_cores(o.i, True) # averaged inter dislocation distance
    r, iK = analyze.intervals(i, o.s) # intervals to display
    f = ['KKKK', 'gggg', 'GaGs'] # functions to calculate
    worker = analyze.calculate(f, o, r) # function results
    master = [average_on_cores(worker[i]) for i in range(len(worker))]
    if rank == root:
        if isinstance(o, sets.Distribution):
            c = str(size) # number of distributions analyzed
            tt = c+" "+o.plotTitle() # plots title
            id = c+"_"+o.fileName() # plots file name
        else:
            c = str(len(o)*size) # number of distributions analyzed
            tt = c+" "+o[0].plotTitle() # plots title
            id = c+"_"+o[0].fileName() # plots file name
        analyze.plot_KKKK(r[:iK], master[f.index('KKKK')].T[:iK].T, p, id, tt)
        analyze.plot_gggg(r, master[f.index('gggg')], p, id, tt)
        analyze.plot_GaGs(r, master[f.index('GaGs')], p, id, tt)
