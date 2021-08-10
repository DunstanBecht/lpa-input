#!/usr/bin/env python
# coding: utf-8

"""
Tools for parallelizing the spatial analysis of distributions.
"""

from mpi4py import MPI
name = MPI.Get_processor_name()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
root = 0

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
    if rank == root:
        m = np.zeros_like(w)
    else:
        m = None
    comm.Reduce([w, MPI.DOUBLE], [m, MPI.DOUBLE], op=MPI.SUM, root=root)
    if rank == root:
        m = m/size
    if b:
        m = comm.bcast(m, root=root)
    return m

def export(
    o: Union[sets.Distribution, sets.Sample],
    p: str = "",
    n: Optional[str] = None,
    t: Optional[str] = None,
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
    i = average_on_cores(o.i, True) # averaged inter dislocation distance
    r, iK = analyze.intervals(i, o.s)
    f = ['KKKK', 'gggg', 'GaGs']
    worker = analyze.calculate(f, o, r)
    master = [average_on_cores(worker[i]) for i in range(len(worker))]
    if rank == root:
        if isinstance(o, sets.Distribution):
            c = str(size)
            t = c+" "+o.plotTitle()
            n = c+"_"+o.fileName()
        else:
            c = str(len(o)*size)
            t = c+" "+o[0].plotTitle()
            n = c+"_"+o[0].fileName()
        analyze.plot_KKKK(r[:iK], master[f.index('KKKK')].T[:iK].T, p, n, t)
        analyze.plot_gggg(r, master[f.index('gggg')], p, n, t)
        analyze.plot_GaGs(r, master[f.index('GaGs')], p, n, t)
