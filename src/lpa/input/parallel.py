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

@beartype
def average_on_cores(
    w: AnalysisOutput,
    b: bool = False
) -> Optional[AnalysisOutput]:
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
) -> Optional[tuple]:
    """
    Export a complete pooled analysis of the object o of each core.

    Function similar to the export function of the analyze module.

    Input:
        o (Distribution|Sample): object to analyze on the core
      **expdir (str): export directory (default: '')
      **expfmt (str): export format (default: 'pdf')
      **expstm (str): export stem (default: o.name(fstm, c='stm'))
      **figttl (str): figure title (default: o.name(fttl, c='ttl'))
      **edgcon (str): edge consideration (default: 'NEC')
      **intrad (ScalarList): interval of radii [nm] (default: ROI size)
      **savtxt (bool): save data to text files (default: False)

    Output:
        r (ScalarList): radius of the neighborhoods [nm]
        KKKK (ScalarListList): stacked K++, K-+, K+-, K-- values [nm^n]
        gggg (ScalarListList): stacked g++, g-+, g+-, g-- values [1]
        GaGs (ScalarListList): stacked Ga and Gs values [m^-1]
    """
    if rank==root and not o.S is None:
        msg = (f"chosen random seed detected on {o.name('stm')} "
               f"(The use of a seed is not recommended. In order for "
               f"the distributions or samples not to be identical "
               f"from one core to another, the seed must not be "
               f"identical from one core to another.)")
        warnings.warn(msg, Warning)
    rmax = o.s
    if o.g == 'circle':
        rmax *= 2
    edgcon = getkwa('edgcon', kwargs, str, 'NEC')
    intrad = getkwa('intrad', kwargs, ScalarList, np.linspace(0, rmax, 200))
    i = average_on_cores(o.i, True) # averaged inter dislocation distance
    fun = ('KKKK', 'gggg', 'GaGs') # functions to calculate
    worker = analyze.calculate(fun, o, intrad, ec=edgcon) # function results
    master = [average_on_cores(worker[i]) for i in range(len(worker))]
    if rank == root:
        if isinstance(o, sets.Distribution):
            c = str(size) # number of distributions analyzed
        else:
            c = str(len(o)*size) # number of distributions analyzed
        # optional parameters
        expstm = getkwa('expstm', kwargs, str, c+"_"+o.name('dmgsS', c='stm'))
        figttl = getkwa('title', kwargs, str, c+" "+o.name('mgsd', c='ttl'))
        savtxt = getkwa('savtxt', kwargs, bool, False)
        kwargs['figttl'] = figttl
        kwargs['edgcon'] = edgcon
        # export
        KKKK, gggg, GaGs = master
        KKKKstm = expstm+"_KKKK_"+edgcon
        ggggstm = expstm+"_gggg_"+edgcon
        GaGsstm = expstm+"_GaGs_"+edgcon
        analyze.plot_KKKK(intrad, KKKK, **kwargs, expstm=KKKKstm)
        analyze.plot_gggg(intrad, gggg, **kwargs, expstm=ggggstm)
        analyze.plot_GaGs(intrad, GaGs, **kwargs, expstm=GaGsstm)
        if savtxt:
            pth = getkwa('expdir', kwargs, str, '')
            np.savetxt(os.path.join(pth, expstm+"_radii.txt"), intrad)
            np.savetxt(os.path.join(pth, KKKKstm+'.txt'), KKKK)
            np.savetxt(os.path.join(pth, ggggstm+'.txt'), gggg)
            np.savetxt(os.path.join(pth, GaGsstm+'.txt'), GaGs)
        return intrad, KKKK, gggg, GaGs
