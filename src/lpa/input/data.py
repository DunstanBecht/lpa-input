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
    """
    ng = np.linalg.norm(g)
    nl = np.linalg.norm(l)
    psi = np.arccos(np.dot(g, l)/(ng*nl))
    if t == 'screw':
        C = np.sin(psi)**2 * np.cos(psi)**2
    elif t == 'edge':
        pg = g - (np.dot(g, l)/nl**2)*l
        pb = b - (np.dot(b, l)/nl**2)*l
        npg = np.linalg.norm(pg)
        npb = np.linalg.norm(pb)
        gamma = np.arccos(np.dot(pg, pb)/(npg*npb))
        C = (np.sin(psi)**4
            / (8*(1-nu)**2)
            * (1-4*nu+8*nu**2+4*(1-2*nu)*np.cos(gamma)**2))
    return C

@beartype
def export_distribution(
    d: sets.Distribution,
    i: Scalar,
    p: str = "",
    n: Optional[str] = None,
    g: Vector = np.array([2, 0, 0]),
    b: Vector = np.array([1, 1, 0]),
) -> None:
    """
    Export the dislocations of d to a standardized input data file.

    The inter-dislocation distance must be specified because it may be
    different from d.i when averaged over several distributions.

    Input:
        d: distribution to be exported
        i: inter dislocation distance [nm]
        p: path where to place the exported file
        n: name of the exported file
        g: diffraction vector direction (hkl)
        b: Burgers vector direction [uvw]

    Complexity:
        O( len(d) )
    """
    if p!="" and p[-1]!="/":
        p += "/"
    if n == None:
        n = d.fileName()
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
        if d.c!=None and 'pbcr' in d.c:
            str_g = "Square_"+d.c[4:]+" side"
        else:
            str_g = "Square side"

    # export
    with open(p+n+".dat", "w") as f:
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
    p: str = "",
    n: Optional[str] = None,
) -> None:
    """
    Export a standardized input data file for each distributions of s.

    Input:
        s: sample of distributions to be exported
        p: path where to place the exported file
        n: name of the exported directory

    Complexity:
        O( len(s) * complexity(export_distribution) )
    """
    if p!="" and p[-1]!="/":
        p += "/"
    e = p+s.fileName()
    # export
    m = len(str(len(s)))
    @beartype
    def fmt(i: int) -> str:
        str_i = str(i+1)
        return "0"*(m-len(str_i))+str_i
    if os.path.exists(e):
        for f in os.listdir(e):
            os.remove(e+"/"+f)
    else:
        os.mkdir(e)
    for i in range(len(s)):
        export_distribution(s[i], s.i, e+"/", fmt(i))

@beartype
def export(
    o: Union[sets.Distribution, sets.Sample],
    p: str = "",
    n: Optional[str] = None,
) -> None:
    """
    Convenient function for exporting a standardized data file.

    Input:
        o: distribution or sample of distributions to export
        p: path where to place the exported file
        n: name of the exported file or directory

    Complexity:
        O( complexity(export_distribution) ) if o is a distribution
        O( complexity(export_sample) ) if o is a sample
    """
    if isinstance(o, sets.Distribution):
        export_distribution(o, p, n)
    else:
        export_sample(o, p, n)
