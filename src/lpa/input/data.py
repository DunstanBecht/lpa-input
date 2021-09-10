#!/usr/bin/env python
# coding: utf-8

"""
Standardized files for input to an X-ray diffraction simulation program.
"""

import os
from . import *
from . import sets

@beartype
def indices(
    v: Vector
) -> str:
    """
    Return the formated string of Miller indices of v.

    Input:
        v (Vector): Miller indices

    Output:
        s (str): formated string of Miller indices
    """
    sep, fmt = " ", '2.0f'
    return sep.join([format(c, fmt) for c in v])

@beartype
def contrast_factor(
    t: str,
    g: Vector,
    l: Vector,
    b: Vector,
    nu: Scalar,
) -> Scalar:
    """
    Return the contrast factor in elastically isotropic crystals.

    Input:
        t (str): dislocation type
        g (Vector): diffraction vector direction (hkl)
        l (Vector): dislocation line vector direction [uvw]
        b (Vector): Burgers vector direction [uvw]
        nu (Scalar): Poisson's number [1]

    Output:
        C (Scalar): dislocation contrast factor [1]

    Input example:
        t = 'screw'
        g = np.array([2, 0, 0])
        l = np.array([1, 1, 0])
        b = np.array([1, 1, 0])
        nu = 0.345

    Output example:
        C = 0.25
    """
    ng = np.linalg.norm(g) # diffraction vector norm to normalize
    nl = np.linalg.norm(l) # dislocation line vector norm to normalize
    psi = np.arccos(np.dot(g, l)/(ng*nl)) # angle between g and l
    if t == 'screw':
        C = np.sin(psi)**2 * np.cos(psi)**2
    elif t == 'edge':
        pg = g - (np.dot(g, l)/nl**2)*l # projection of g
        pb = b - (np.dot(b, l)/nl**2)*l # projection of b
        npg = np.linalg.norm(pg) # norm of the projection of g
        npb = np.linalg.norm(pb) # norm of the projection of b
        gamma = np.arccos(np.dot(pg, pb)/(npg*npb)) # angle betw. projections
        C = (np.sin(psi)**4
            / (8*(1-nu)**2)
            * (1-4*nu+8*nu**2+4*(1-2*nu)*np.cos(gamma)**2))
    return C

dft_l = {'edge': np.array([ 1, -1, -2]), 'screw': np.array([ 1,  1,  0])}
dft_L = {'edge': np.array([ 1,  1,  0]), 'screw': np.array([-1,  1,  0])}

@beartype
def export_distribution(
    d: sets.Distribution,
    i: Scalar,
    **kwargs,
) -> None:
    """
    Export the dislocations of d to a standardized input data file.

    The inter-dislocation distance must be specified because it may be
    different from d.i when averaged over several distributions.

    Input:
        d (Distribution): distribution to be exported
        i (Scalar): inter dislocation distance [nm]
      **g (Vector): diffraction vector direction (hkl) (default: [2, 0, 0])
      **b (Vector): Burgers vector direction [uvw] (default: [2, 0, 0])
      **l (Vector): line vector direction [uvw] (default: depends on d. type)
      **L (Vector): direction along Lx [uvw] (defaut: depends on d. type)
      **a (Scalar): cell side [nm] (default: 0.40494)
      **a3 (Scalar): step size along Lx [nm] (default: depends on i)
      **nu (Scalar): Poisson's number [1] (default: 0.345)
      **expdir (str): export directory (default: '')
      **expfmt (str): export format (default: 'dat')
      **expstm (str): export stem (default: d.name())

    Complexity:
        O( len(d) )
    """
    # optional parameters
    g = getkwa('g', kwargs, Vector, np.array([2, 0, 0]))
    b = getkwa('b', kwargs, Vector, np.array([1, 1, 0]))
    l = getkwa('l', kwargs, Vector, dft_l[d.t])
    L = getkwa('L', kwargs, Vector, dft_L[d.t])
    a = getkwa('a', kwargs, Scalar, 0.40494)
    a3 = getkwa('a3', kwargs, Scalar, max(2, i/12))
    nu = getkwa('nu', kwargs, Scalar, 0.345)
    expdir = getkwa('expdir', kwargs, str, '')
    expfmt = getkwa('expfmt', kwargs, str, 'dat')
    expstm = getkwa('expstm', kwargs, str, d.name(c='stm'))
    endkwa(kwargs)
    # parameters
    if d.t == 'screw' and np.linalg.norm(np.cross(l, b)) != 0:
        raise ValueError("screw type but l and b not parallel")
    elif d.t == 'edge' and np.dot(l, b) != 0:
        raise ValueError("edge type but l and b not perpendicular")
    C = contrast_factor(d.t, g, l, b, nu)
    if d.g == 'circle':
        str_g = "Cylinder radius"
    if d.g == 'square':
        if d.c and d.c[:4]=='PBCR':
            str_g = "Square_"+d.c[4:]+" side"
        else:
            str_g = "Square side"
    # write
    with open(os.path.join(expdir, expstm+"."+expfmt), "w") as f:
        h = ("# please keep the structure of this file unchanged\n"
            + indices(l)+" # z: direction of 'l' (line vector) [uvw]\n"
            + indices(L)+" # x: direction of 'L' (Fourier variable) [uvw]\n"
            + indices(b)+" # b: Burgers vector direction [uvw]\n"
            + indices(g)+" # g: diffraction vector direction (hkl)\n"
            + format(C,      '8.6f')+" # C: contrast coefficient [1]\n"
            + format(a,      '8.6f')+" # a: cell parameter [nm]\n"
            + format(d.s,    '8.0f')+" # s: "+str_g+" [nm]\n"
            + format(a3,     '8.1f')+" # a3: step size of 'L' along x [nm]\n"
            + format(nu,     '8.3f')+" # nu: Poisson's number [1]\n"
            + format(len(d), '8.0f')+" # number of dislocations\n"
            + "# Burgers vector and dislocation (x,y) coordinates\n")
        f.write(h)
        fmt = "%2.0f %22.15E %22.15E"
        np.savetxt(f, np.stack((d.b, d.p[:,0], d.p[:,1])).T, fmt=fmt)

@beartype
def export_sample(
    s: sets.Sample,
    **kwargs,
) -> None:
    """
    Export a standardized input data file for each distributions of s.

    Input:
        s (Sample): sample of distributions to be exported
        <see export_distribution function for keyword arguments>

    Complexity:
        O( len(s) * complexity(export_distribution) )
    """
    # optional parameters
    expdir = kwargs.pop('expdir', '') # export directory
    expfmt = kwargs.pop('expfmt', 'dat') # export format
    expstm = kwargs.pop('expstm', s.name(c='stm')) # export stem
    # export
    stmdir = os.path.join(expdir, expstm) # folder where to export the files
    w = len(str(len(s))) # number of characters in file names
    if os.path.exists(stmdir):
        for f in os.listdir(stmdir):
            os.remove(os.path.join(stmdir, f)) # delete the existing
    else:
        os.mkdir(stmdir)
    for i in range(len(s)):
        export_distribution(
            s[i],
            s.i,
            expdir=stmdir,
            expfmt=expfmt,
            expstm=str(i+1).zfill(w),
            **kwargs,
        )

@beartype
def export(
    o: Union[sets.Distribution, sets.Sample],
    **kwargs,
) -> None:
    """
    Versatile function for exporting a standardized data file.

    Input:
        o (Distribution|Sample): distribution or sample to export
        <see export_distribution function for keyword arguments>

    Complexity:
        O( complexity(export_distribution) ) if o is a distribution
        O( complexity(export_sample) ) if o is a sample
    """
    if isinstance(o, sets.Distribution):
        export_distribution(o, o.i,**kwargs)
    else:
        export_sample(o, **kwargs)
