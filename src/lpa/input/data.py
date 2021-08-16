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
        v: Miller indices

    Output:
        s: formated string of Miller indices
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
        t: dislocation type
        g: diffraction vector direction (hkl)
        l: dislocation line vector direction [uvw]
        b: Burgers vector direction [uvw]
        nu: Poisson's number [1]

    Output:
        C: dislocation contrast factor [1]

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

@beartype
def export_distribution(
    d: sets.Distribution,
    i: Scalar,
    g: Vector,
    b: Vector,
    exdir: str = "",
    exfmt: str = "dat",
    exstm: Optional[str] = None,
) -> None:
    """
    Export the dislocations of d to a standardized input data file.

    The inter-dislocation distance must be specified because it may be
    different from d.i when averaged over several distributions.

    Input:
        d: distribution to be exported
        i: inter dislocation distance [nm]
        g: diffraction vector direction (hkl)
        b: Burgers vector direction [uvw]
        exdir: export directory
        exfmt: export format
        exstm: export stem

    Complexity:
        O( len(d) )
    """
    if exstm is None:
        exstm = d.name(c='stm', s=True)
    # parameters
    m = 12 # number of points along Lx
    a = 0.40494 # cell side [nm]
    a3 = max(2, i/m) # step size along Lx [nm]
    if d.g == 'circle':
        a3 = min(a3, 3*i/m)
    if d.t == 'screw':
        l = np.array([ 1,  1,  0]) # dislocation line vector direction [uvw]
        L = np.array([-1,  1,  0]) # direction along Lx [uvw]
        if np.linalg.norm(np.cross(l, b)) != 0:
            raise ValueError("screw type but l and b not parallel")
    elif d.t == 'edge':
        l = np.array([ 1, -1, -2]) # dislocation line vector direction [uvw]
        L = np.array([ 1,  1,  0]) # direction along Lx [uvw]
        if np.dot(l, b) != 0:
            raise ValueError("edge type but l and b not perpendicular")
    nu = 0.345 # Poisson's number
    C = contrast_factor(d.t, g, l, b, nu)
    if d.g == 'circle':
        str_g = "Cylinder radius"
    if d.g == 'square':
        if d.c and d.c[:4]=='PBCR':
            str_g = "Square_"+d.c[4:]+" side"
        else:
            str_g = "Square side"
    # write
    with open(os.path.join(exdir, exstm+"."+exfmt), "w") as f:
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
        fmt = "%2.0f %27.20E %27.20E"
        np.savetxt(f, np.stack((d.b, d.p[:,0], d.p[:,1])).T, fmt=fmt)

@beartype
def export_sample(
    s: sets.Sample,
    g: Vector,
    b: Vector,
    exdir: str = "",
    exfmt: str = "dat",
    exstm: Optional[str] = None,
) -> None:
    """
    Export a standardized input data file for each distributions of s.

    Input:
        s: sample of distributions to be exported
        g: diffraction vector direction (hkl)
        b: Burgers vector direction [uvw]
        exdir: export directory
        exfmt: export format
        exstm: export stem

    Complexity:
        O( len(s) * complexity(export_distribution) )
    """
    if exstm is None:
        exstm = s.name(c='stm', s=True)
    stmdir = os.path.join(exdir, exstm) # folder where to export the files
    # export
    w = len(str(len(s))) # number of characters in file names
    if os.path.exists(stmdir):
        for f in os.listdir(stmdir):
            os.remove(os.path.join(stmdir, f)) # delete the existing
    else:
        os.mkdir(stmdir)
    for i in range(len(s)):
        export_distribution(s[i], s.i, g, b, stmdir, exfmt, str(i+1).zfill(w))

@beartype
def export(
    o: Union[sets.Distribution, sets.Sample],
    g: Vector = np.array([2, 0, 0]),
    b: Vector = np.array([1, 1, 0]),
    exdir: str = "",
    exfmt: str = "dat",
    exstm: Optional[str] = None,
) -> None:
    """
    Convenient function for exporting a standardized data file.

    Input:
        o: distribution or sample of distributions to export
        g: diffraction vector direction (hkl)
        b: Burgers vector direction [uvw]
        exdir: export directory
        exfmt: export format
        exstm: export stem

    Complexity:
        O( complexity(export_distribution) ) if o is a distribution
        O( complexity(export_sample) ) if o is a sample
    """
    if isinstance(o, sets.Distribution):
        export_distribution(o, o.i, g, b, exdir, exfmt, exstm)
    else:
        export_sample(o, g, b, exdir, exfmt, exstm)
