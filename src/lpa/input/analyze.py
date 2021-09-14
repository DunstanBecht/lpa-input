#!/usr/bin/env python
# coding: utf-8

"""
Tools for spatial analysis of dislocation distributions.
"""

import matplotlib.pyplot as plt
import scipy.special
from . import *
from . import sets
from . import geometries

@beartype
def N(
    a: Vector,
    B: VectorList,
    r2: ScalarList,
) -> ScalarList:
    """
    Return the number of points of B in the neighborhoods of a.

    The function does not use the notion of dislocation. It allows only
    a spatial analysis of a distribution of points in space. If a point
    of B coincides with the point a then it is not counted.

    Input:
        a (Vector): point around which the neighborhoods are formed
        B (VectorList): points that can be counted when in the neighborhoods
        r2 (ScalarList): squared neighborhood radii in ascending order

    Output:
        n (ScalarList): list of the number of points for each radius value

    Input example:
        a = np.array([a_x, a_y])
        B = np.array([[b_0_x, b_0_y], [b_1_x, b_1_y], ...])
        r2 = np.array([r_0^2, r_1^2, r_2^2, ...])

    Output example:
        n = np.array([N(r_0), (r_1), N(r_2), ...])

    Complexity:
        O( len(B)*log(len(B)) + len(r2) )
    """
    d2 = np.sum(np.square(np.subtract(B, a)), axis=1) # squared distances to a
    d2 = np.sort(d2[(d2<=r2[-1]) & (d2>0)]) # sorted squared distances to a
    j, dn = 0, np.zeros(len(r2), dtype=np.int64)
    for i in range(len(d2)):
        while d2[i]>r2[j]:
            j += 1
        dn[j] += 1
    return np.cumsum(dn)

@beartype
def M(
    A: VectorList,
    B: VectorList,
    w: CorrectionFunction,
    r: ScalarList,
    r2: ScalarList,
) -> ScalarList:
    """
    Return the value of N corrrected and averaged on the points of A.

    The function does not use the notion of dislocation. It allows only
    a spatial analysis of a distribution of points in space.

    Input:
        A (VectorList): points around which the neighborhoods are formed
        B (VectorList): points that can be counted when in the neighborhoods
        w (CorrectionFunction): weighting function for correction at the edges
        r (ScalarList): neighborhood radii in ascending order
        r2 (ScalarList): squared neighborhood radii in ascending order

    Output:
        m (ScalarList): average number of points for each radius value

    Input example:
        A = np.array([[a_0_x, a_0_y], [a_1_x, a_1_y], ...])
        B = np.array([[b_0_x, b_0_y], [b_1_x, b_1_y], ...])
        w = lambda a, r, r2: np.ones(len(r))
        r = np.array([r_0, r_1, r_2, ...])
        r2 = np.array([r_0^2, r_1^2, r_2^2, ...])

    Output example:
        m = np.array([M(r_0), M(r_1), M(r_2)])

    Complexity:
        O( len(A) * (complexity(w)+complexity(N)) )
    """
    s = np.zeros(len(r))
    for i in range(len(A)):
        s += w(A[i], r, r2) * N(A[i], B, r2)
    return s/len(A)

@beartype
def MMMM_cp_cm(
    d: sets.Distribution,
    r: ScalarList,
    r2: ScalarList,
) -> tuple:
    """
    Return M++, M-+, M+-, M-- and the number of + and - dislocations.

    M++ is the average number of dislocations with positive Burgers
    vector sense around dislocations with positive sense. The
    dislocation in the center of the neighborhood is not counted.

    M-+ is the average number of dislocations with positive Burgers
    vector sense around dislocations with negative sense.

    M+- is the average number of dislocations with negative Burgers
    vector sense around dislocations with positive sense.

    M-- is the average number of dislocations with negative Burgers
    vector sense around dislocations with negative sense. The
    dislocation in the center of the neighborhood is not counted.

    Input:
        d (Distribution): distribution of dislocations to analyze
        r (ScalarList): neighborhood radii in ascending order [nm]
        r2 (ScalarList): squared neighborhood radii in ascending order [nm^2]

    Output:
        MMMM (ScalarListList): stacked values of M++, M-+, M+-, M-- [1]
        cp (Scalar): number of dislocations with Burgers vector sense + [1]
        cm (Scalar): number of dislocations with Burgers vector sense - [1]

    Output example:
        MMMM = np.array([
            [M++(r_0), M++(r_1), M++(r_2), ...],
            [M-+(r_0), M-+(r_1), M-+(r_2), ...],
            [M+-(r_0), M+-(r_1), M+-(r_2), ...],
            [M--(r_0), M--(r_1), M--(r_2), ...],
        ])
        cp = 2
        cm = 3

    Complexity:
        O( complexity(M) )
    """
    Pp = d.p[d.b>0] # positive Burgers vector sense dislocation positions
    Pm = d.p[d.b<0] # negative Burgers vector sense dislocation positions
    if d.c:
        Pp = Pp[geometries.mask(d.g, d.s, Pp)]
        Pm = Pm[geometries.mask(d.g, d.s, Pm)]
    Mpp = M(Pp, Pp, d.w, r, r2) # M++
    Mmp = M(Pm, Pp, d.w, r, r2) # M-+
    Mpm = M(Pp, Pm, d.w, r, r2) # M+-
    Mmm = M(Pm, Pm, d.w, r, r2) # M--
    MMMM = np.stack((Mpp, Mmp, Mpm, Mmm))
    cp = len(Pp)
    cm = len(Pm)
    return MMMM, cp, cm

@beartype
def KKKK_dp_dm(
    MMMM: ScalarListList,
    cp: Scalar,
    cm: Scalar,
    v: Scalar,
) -> tuple:
    """
    Return K++, K-+, K+-, K-- and densities of + and - dislocations.

    Implementation of the Ripley’s K function.

    This function is incremental. It uses pre-calculated MMMM values
    to avoid unnecessary repetition of calculations common to spatial
    analysis functions

    Input:
        MMMM (ScalarListList): stacked values of M++, M-+, M+-, M-- [1]
        cp (Scalar): number of dislocations with Burgers vector sense + [1]
        cm (Scalar): number of dislocations with Burgers vector sense - [1]
        v (Scalar): n-volume of the region of interest [nm^n]

    Output:
        KKKK (ScalarListList): stacked values of K++, K-+, K+-, K-- [nm^n]
        dp (Scalar): density of dislocations with Burgers vector + [nm^-n]
        dm (Scalar): density of dislocations with Burgers vector - [nm^-n]

    Output example:
        KKKK = np.array([
            [K++(r_0), K++(r_1), K++(r_2), ...],
            [K-+(r_0), K-+(r_1), K-+(r_2), ...],
            [K+-(r_0), K+-(r_1), K+-(r_2), ...],
            [K--(r_0), K--(r_1), K--(r_2), ...],
        ])
        dp = 0.002
        dp = 0.003

    Complexity:
        O( len(MMMM) )
    """
    dp, dm = cp/v, cm/v
    KKKK = np.concatenate((MMMM[0:2]/dp, MMMM[2:4]/dm))
    return KKKK, dp, dm

@beartype
def gggg_dV_dKKKK(
    KKKK: ScalarListList,
    dr: ScalarList,
    r: ScalarList,
    n: int,
) -> tuple:
    """
    Return g++, g-+, g+-, g-- and calculated differentials.

    Implementation of the pair correlation function.

    This function is incremental. It uses pre-calculated KKKK values
    to avoid unnecessary repetition of calculations common to spatial
    analysis functions

    Input:
        KKKK (ScalarListList): stacked values of K++, K-+, K+-, K-- [nm^n]
        dr (ScalarList): differentials of the radius [nm]
        r (ScalarList): neighborhood radii in ascending order [nm]
        n (int): dimension of the space of the region of interest

    Ouput:
        gggg (ScalarListList): stacked values of g++, g-+, g+-, g-- [1]
        dV (ScalarList): differentials of the neighborhood n-volumes [nm^n]
        dKKKK (ScalarList): differential of KKKK [nm^n]

    Output example:
        gggg = np.array([
            [g++(r_0), g++(r_1), g++(r_2), ...],
            [g-+(r_0), g-+(r_1), g-+(r_2), ...],
            [g+-(r_0), g+-(r_1), g+-(r_2), ...],
            [g--(r_0), g--(r_1), g--(r_2), ...],
        ])
        dV = np.array([pi*r_0*dr_0, pi*r_1*dr_1, ...])
        dKKKK = np.array([
            [K++(r_1)-K++(r_0), (K++(r_2)-K++(r_0))/2, ...],
            [K-+(r_1)-K-+(r_0), (K-+(r_2)-K-+(r_0))/2, ...],
            [K+-(r_1)-K+-(r_0), (K+-(r_2)-K+-(r_0))/2, ...],
            [K--(r_1)-K--(r_0), (K--(r_2)-K--(r_0))/2, ...],
        ])

    Complexity:
        O( len(r) )
    """
    dV = (np.pi**(n/2)/scipy.special.gamma(n/2+1))*n*r**(n-1)*dr # diff. of V
    dV[dV==0] = np.nan # mask zero values
    dKKKK = np.gradient(KKKK, axis=1) # differential of KKKK
    gggg = dKKKK/dV
    return gggg, dV, dKKKK

@beartype
def GaGs_dMMMM(
    MMMM: ScalarListList,
    dr: ScalarList,
    cp: Scalar,
    cm: Scalar,
) -> tuple:
    """
    Return Ga and Gs.

    This function is incremental. It uses pre-calculated MMMM values
    to avoid unnecessary repetition of calculations common to spatial
    analysis functions

    Input:
        MMMM (ScalarListList): stacked values of M++, M-+, M+-, M-- [1]
        dr (ScalarList): differentials of the radius [nm]
        cp (Scalar): number of dislocations with Burgers vector sense + [1]
        cm (Scalar): number of dislocations with Burgers vector sense - [1]

    Output:
        GaGs (ScalarListList): stacked values of Ga and Gs [nm^-1]
        dMMMM (ScalarList): differential of MMMM [1]

    Output example:
        GaGs = np.array([
            [Ga(r_0), Ga(r_1), Ga(r_2), ...],
            [Gs(r_0), Gs(r_1), Gs(r_2), ...],
        ])
        dMMMM = np.array([
            [M++(r_1)-M++(r_0), (M++(r_2)-M++(r_0))/2, ...],
            [M-+(r_1)-M-+(r_0), (M-+(r_2)-M-+(r_0))/2, ...],
            [M+-(r_1)-M+-(r_0), (M+-(r_2)-M+-(r_0))/2, ...],
            [M--(r_1)-M--(r_0), (M--(r_2)-M--(r_0))/2, ...],
        ])

    Complexity:
        O( len(MMMM) )
    """
    dMMMM = np.gradient(MMMM, axis=1)
    dMppdr, dMmpdr, dMpmdr, dMmmdr = dMMMM/dr
    Ga = cp*(dMppdr-dMpmdr) + cm*(dMmmdr-dMmpdr)
    Gs = cp*(dMppdr+dMpmdr) + cm*(dMmmdr+dMmpdr)
    return np.stack((Ga, Gs)), dMMMM

@beartype
def calculate(
    q: list,
    o: Union[sets.Distribution, sets.Sample],
    r: ScalarList,
    **kwargs,
) -> AnalysisOutput:
    """
    Return the values of the quantities in q calculated incrementally.

    The following quantities can be requested:
        'MMMM': stacked values of M++, M-+, M+-, M-- [1]
        'KKKK': stacked values of K++, K-+, K+-, K-- [nm^n]
        'gggg': stacked values of g++, g-+, g+-, g-- [1]
        'GaGs': stacked values of Ga and Gs [1]
        'cp': number of dislocations with Burgers vector sense + [1]
        'cm': number of dislocations with Burgers vector sense - [1]
        'dp': density of dislocations with Burgers vector sense + [nm^-n]
        'dm': density of dislocations with Burgers vector sense - [nm^-n]
        'dV': differentials of the neighborhood volume [nm^n]
        'dKKKK': differentials of K++, K-+, K+-, K-- [nm^n]
        'dMMMM': differentials of M++, M-+, M+-, M-- [1]

    Input:
        q (list): name of the quantities to calculate
        o (Distribution|Sample): distribution or sample to analyze
        r (ScalarList): neighborhood radii in ascending order [nm]
      **r2 (ScalarList): squared neighborhood radii [nm^2]
      **dr (ScalarList): differentials of the neighborhood radii [nm]

    Output:
        v (AnalysisOutput): values of the quantities in requested order

    Exemple:
        KKKK, GaGs = calculate(['KKKK', 'GaGs'], distribution, r)

    Complexity:
        O( complexity(MMMM_Pp_Pm) ) if o is a distribution
        O( complexity(MMMM_Pp_Pm) * len(o) ) if o is a sample
    """
    # optional parameters
    r2 = getkwa('r2', kwargs, ScalarList, np.square(r))
    dr = getkwa('dr', kwargs, ScalarList, np.gradient(r))
    endkwa(kwargs)
    # create an incremental analysis function to be applied to a distribution
    @beartype
    def calculate_on_distribution(
        d: sets.Distribution,
    ) -> AnalysisOutput:
        """
        Auxiliary function applied to a distribution only.
        """
        s = {}
        clcgggg = 'gggg' in q or 'dp' in q or 'dm' in q or 'dKKKK' in q
        clcKKKK = clcgggg or 'KKKK' in q or 'cp' in q or 'cm' in q
        clcGaGs = 'GaGs' in q or 'dMMMM' in q
        s['MMMM'], s['cp'], s['cm'] = MMMM_cp_cm(
            d,
            r,
            r2,
        )
        if clcKKKK:
            s['KKKK'], s['dp'], s['dm'] = KKKK_dp_dm(
                s['MMMM'],
                s['cp'],
                s['cm'],
                d.v,
            )
        if clcgggg or clcGaGs:
            s['dr'] = np.gradient(r)
        if clcgggg:
            s['gggg'], s['dV'], s['dKKKK'] =  gggg_dV_dKKKK(
                s['KKKK'],
                s['dr'],
                r,
                d.n,
            )
        if clcGaGs:
            s['GaGs'], s['dMMMM'] = GaGs_dMMMM(
                s['MMMM'],
                s['dr'],
                s['cp'],
                s['cm'],
            )
        return tuple([s[k] for k in q])
    # return the result of the function applied to a distribution or a sample
    if isinstance(o, sets.Distribution):
        return calculate_on_distribution(o)
    return o.average(calculate_on_distribution)

@beartype
def plot_KKKK(
    r: ScalarList,
    KKKK: ScalarListList,
    expdir: str,
    expfmt: str,
    expstm: str,
    title: str,
) -> None:
    """
    Export the figure showing functions K++, K-+, K+-, K--.

    Input:
        r (ScalarList): radius of the neighborhoods [nm]
        KKKK (ScalarListList): stacked values of K++, K-+, K+-, K-- [nm^n]
        expdir (str): export directory
        expfmt (str): export format
        expstm (str): export stem
        title (str): figure title

    Complexity:
        O( len(r) )
    """
    # K of a uniform distribution
    r_compare = np.linspace(r[0], r[-1], 10)
    k_compare = np.pi*r_compare**2
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.1)
    fig.suptitle(title, fontsize=16)
    # ax1
    ax1.plot(r, KKKK[0], label=r"$K_{++}(r)$")
    ax1.plot(r, KKKK[1], label=r"$K_{-+}(r)$")
    ax1.plot(r_compare, k_compare, "*", label=r"$ \pi r^2 $", color='black')
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r"$r \ (nm)$")
    # ax2
    ax2.plot(r, KKKK[3], label=r"$K_{--}(r)$")
    ax2.plot(r, KKKK[2], label=r"$K_{+-}(r)$")
    ax2.plot(r_compare, k_compare, "*", label=r"$ \pi r^2 $", color='black')
    ax2.legend()
    ax2.grid()
    ax2.set_xlabel(r"$r \ (nm)$")
    ax2.set_ylabel(r"$(nm^2)$")
    # export
    plt.savefig(os.path.join(expdir, expstm+"_KKKK."+expfmt), format=expfmt)
    plt.close('all')

@beartype
def plot_gggg(
    r: ScalarList,
    gggg: ScalarListList,
    expdir: str,
    expfmt: str,
    expstm: str,
    title: str,
) -> None:
    """
    Export the figure showing functions g++, g-+, g+-, g--.

    Input:
        r (ScalarList): radius of the neighborhoods [nm]
        gggg (ScalarListList): stacked values of g++, g-+, g+-, g-- [1]
        expdir (str): export directory
        expfmt (str): export format
        expstm (str): export stem
        title (str): figure title

    Complexity:
        O( len(r) )
    """
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.1)
    fig.suptitle(title, fontsize=16)
    ymin = min(np.nanmin(gggg[gggg != -np.inf])*0.95, 0.9)
    ymax = max(np.nanmax(gggg[gggg != np.inf])*1.05, 1.1)
    # ax1
    ax1.plot(r, gggg[0], label=r"$g_{++}(r)$")
    ax1.plot(r, gggg[1], label=r"$g_{-+}(r)$")
    ax1.hlines(1, r[0], r[-1], label=r"$1$", color='black')
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r"$r \ (nm)$")
    ax1.set_ylim(ymin, ymax)
    # ax2
    ax2.plot(r, gggg[3], label=r"$g_{--}(r)$")
    ax2.plot(r, gggg[2], label=r"$g_{+-}(r)$")
    ax2.hlines(1, r[0], r[-1], label=r"$1$", color='black')
    ax2.legend()
    ax2.grid()
    ax2.set_xlabel(r"$r \ (nm)$")
    ax2.set_ylim(ymin, ymax)
    # export
    plt.savefig(os.path.join(expdir, expstm+"_gggg."+expfmt), format=expfmt)
    plt.close('all')

@beartype
def plot_GaGs(
    r: ScalarList,
    GaGs: ScalarListList,
    expdir: str,
    expfmt: str,
    expstm: str,
    title: str,
) -> None:
    """
    Export the figure showing functions Ga and Gs.

    Input:
        r (ScalarList): radius of the neighborhoods [nm]
        GaGs (ScalarListList): stacked values of Ga and Gs [1]
        expdir (str): export directory
        expfmt (str): export format
        expstm (str): export stem
        title (str): figure title

    Complexity:
        O( len(r) )
    """
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.1)
    fig.suptitle(title, fontsize=16)
    # ax1
    ax1.plot(r, GaGs[0], label=r"$G_A(r)$")
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r"$r \ (nm)$")
    # ax2
    ax2.plot(r, GaGs[1], label=r"$G_S(r)$")
    ax2.legend()
    ax2.grid()
    ax2.set_xlabel(r"$r \ (nm)$")
    # export
    plt.savefig(os.path.join(expdir, expstm+"_GaGs."+expfmt), format=expfmt)
    plt.close('all')

@beartype
def intervals(
    i: Scalar,
    s: Scalar,
    n: int = 100,
) -> tuple:
    """
    Suggest default radius interval for the analysis of a distribution.

    It is proposed to conduct the analysis for a radius varying from 0
    to half the characteristic size of the region of interest. For the
    K-function plot it is proposed to vary the radius from 0 to 4 times
    the inter-dislocation distance.

    Input:
        i (Scalar): inter dislocation distance [nm]
        s (Scalar): shape size [nm]
        n (int): number of steps between a and b

    Output:
        r (ScalarList): full range of radii for analysis [nm]
        iK (int): index of the maximum value for the plot of K
    """
    a, b, c = 0, 4*i, s/2
    return np.linspace(a, max(c, b), max(int((c-a)*n/(b-a)), n)), n-1

@beartype
def export(
    o: Union[sets.Distribution, sets.Sample],
    **kwargs,
) -> None:
    """
    Export a complete analysis of the object o.

    Input:
        o: distribution or sample of distributions to analyze
      **expdir: export directory (default: '')
      **expfmt: export format
      **expstm: export stem
      **title: figure title

    Complexity:
        O( complexity(calculate) )
    """
    # optional parameters
    fstm, fttl = 'dmgsS', 'mgsd'
    if isinstance(o, sets.Sample):
        fstm, fttl = "n"+fstm, "n"+fttl
    expdir = getkwa('expdir', kwargs, str, '')
    expfmt = getkwa('expfmt', kwargs, str, 'pdf')
    expstm = getkwa('expstm', kwargs, str, o.name(fstm, c='stm'))
    title = getkwa('title', kwargs, str, o.name(fttl, c='ttl'))
    endkwa(kwargs)
    # export
    r, iK = intervals(o.i, o.s) # range of the study
    KKKK, gggg, GaGs = calculate(['KKKK', 'gggg', 'GaGs'], o, r)
    args = (expdir, expfmt, expstm, title)
    plot_KKKK(r[:iK], KKKK.T[:iK].T, *args)
    plot_gggg(r, gggg, *args)
    plot_GaGs(r, GaGs, *args)
