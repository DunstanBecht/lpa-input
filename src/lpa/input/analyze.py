#!/usr/bin/env python
# coding: utf-8

"""
Tools for spatial analysis of dislocation distributions.
"""

import matplotlib.pyplot as plt
import scipy.special
from . import *
from . import __version__
from . import sets
from . import geometries
from . import boundaries

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
        B (VectorList): points observed and potentially counted
        r2 (ScalarList): squared neighborhood radii in ascending order

    Output:
        n (ScalarList): number of points for each radius value

    Input example:
        a = np.array([a_x, a_y])
        B = np.array([[b_0_x, b_0_y], [b_1_x, b_1_y], ...])
        r2 = np.array([r_0^2, r_1^2, r_2^2, ...])

    Output example:
        n = np.array([N(a, B, r_0), N(a, B, r_1), N(a, B, r_2), ...])

    Complexity:
        O( len(B)*log(len(B)) + len(r2) )
    """
    d2 = np.sum(np.square(np.subtract(B, a)), axis=1) # squared distances to a
    d2 = np.sort(d2[(d2<=r2[-1]) & (d2>0)]) # sorted and filtered sqrd. dist.
    j = 0 # index to browse the neighborhood radii
    dn = np.zeros(len(r2), dtype=np.int64) # differential of n
    for i in range(len(d2)): # browse the observed points
        while d2[i]>r2[j]: # as long as the point is beyond the current radius
            j += 1 # take the following radius
        dn[j] += 1 # when the point is in the neighborhood, it is counted
    return np.cumsum(dn) # cumulative sum of differentials

@beartype
def M(
    A: VectorList,
    B: VectorList,
    w: CorrectionFunction,
    r: ScalarList,
    r2: ScalarList,
) -> ScalarList:
    """
    Return the value of N(a, B, r2) averaged over the points a in A.

    The function does not use the notion of dislocation. It allows only
    a spatial analysis of a distribution of points in space.

    Input:
        A (VectorList): points of the centers of the neighborhoods
        B (VectorList): points observed and potentially counted
        w (CorrectionFunction): weighting function for edge correction
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
        m = np.array([M(A, B, r_0), M(A, B, r_1), M(A, B, r_2)])

    Complexity:
        O( len(A) * (complexity_of(w)+complexity_of(N)) )
    """
    sumN = np.zeros(len(r)) # sum of the results of N
    sumw = np.zeros(len(r)) # weight sum
    for i in range(len(A)): # browse the center points of the neighborhoods
        sumN += N(A[i], B, r2)
        sumw += w(A[i], r, r2)
    return sumN/sumw

@beartype
def MMMM_cp_cm(
    d: sets.Distribution,
    r: ScalarList,
    r2: ScalarList,
    ec: str,
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

    When the parameter ec equals 'NEC' no edge correction is performed.
    When it equals 'WOA' the results are weighted by edge correction
    coefficients. When it equals 'PBC' or 'GBB', corresponding boundary
    conditions are applied with a rank of replication determined by the
    maximum value of the neighborhood radius.

    Input:
        d (Distribution): distribution of dislocations to analyze
        r (ScalarList): neighborhood radii in ascending order [nm]
        r2 (ScalarList): squared neighborhood radii [nm^2]
        ec (str): edge consideration

    Output:
        MMMM (ScalarListList): stacked M++, M-+, M+-, M-- values [1]
        cp (Scalar): number of dislocations with sense + [1]
        cm (Scalar): number of dislocations with sense - [1]

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
        O( complexity_of(M) )
    """
    # masks
    mskroi = geometries.mask(d.g, d.s, d.p) # the dislocation is in the ROI
    mskbsp = d.b>0 # the dislocation has a Burgers vector sense +
    mskbsm = d.b<0 # the dislocation has a Burgers vector sense -
    # neighbourhood centers
    Pp1 = d.p[mskbsp & mskroi] # positions is the ROI with sense +
    Pm1 = d.p[mskbsm & mskroi] # positions is the ROI with sense -
    # observed position
    if ec == 'NEC': # no edge correction
        Pp2 = d.p[mskbsp] # positions with sense +
        Pm2 = d.p[mskbsm] # positions with sense +
    elif ec == 'WOA': # weighting by overlapping area
        Pp2 = Pp1 # observed positions with sense +
        Pm2 = Pm1 # observed positions with sense -
    elif ec=='PBC' or ec=='GBB': # apply boundary conditions
        rep = int(np.ceil(np.max(r)/d.s)) # replication rank
        cp, cb = d.conditions(ec+str(rep)) # outer dislocations
        P2 = np.concatenate((d.p, cp)) # observed positions
        B2 = np.concatenate((d.b, cb)) # observed Burgers vector senses
        Pp2 = P2[B2>0] # observed positions with sense +
        Pm2 = P2[B2<0] # observed positions with sense -
        if ec=='PBC' and r[-1]>=d.s:
            msg = ("discontinuities of M++ and M-- when using PBC "
                   "(This is due to the periodic presence of the same "
                   "dislocation at the same place in the replicated "
                   "regions. The dislocation is suddenly counted "
                   "multiple times at radius values corresponding to "
                   "multiples of d.s, sqrt(2)*d.s etc...)")
            warnings.warn(msg, Warning)
    else:
        raise ValueError(f"invalid edge consideration: {ec}")
    # weighting
    if ec == 'WOA': # no edge correction
        w = d.w # weighting function
    else:
        w = lambda a, r, r2: 1 # weighting function
    # calculate
    Mpp = M(Pp1, Pp2, w, r, r2) # M++
    Mmp = M(Pm1, Pp2, w, r, r2) # M-+
    Mpm = M(Pp1, Pm2, w, r, r2) # M+-
    Mmm = M(Pm1, Pm2, w, r, r2) # M--
    MMMM = np.stack((Mpp, Mmp, Mpm, Mmm)) # stacked M++, M-+, M+-, M--
    cp = len(Pp1) # number of dislocations with sense +
    cm = len(Pm1) # number of dislocations with sense -
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

    Implementation of the Ripley???s K function.

    This function is incremental. It uses pre-calculated MMMM values
    to avoid unnecessary repetition of calculations common to spatial
    analysis functions

    Input:
        MMMM (ScalarListList): stacked M++, M-+, M+-, M-- values [1]
        cp (Scalar): number of dislocations with sense + [1]
        cm (Scalar): number of dislocations with sense - [1]
        v (Scalar): n-volume of the region of interest [nm^n]

    Output:
        KKKK (ScalarListList): stacked K++, K-+, K+-, K-- values [nm^n]
        dp (Scalar): density of dislocations with sense + [nm^-n]
        dm (Scalar): density of dislocations with sense - [nm^-n]

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
    dp, dm = cp/v, cm/v # densities of + and - dislocations
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
    Return g++, g-+, g+-, g-- and differentials of V and KKKK.

    Implementation of the pair correlation function.

    This function is incremental. It uses pre-calculated KKKK values
    to avoid unnecessary repetition of calculations common to spatial
    analysis functions

    Input:
        KKKK (ScalarListList): stacked K++, K-+, K+-, K-- values [nm^n]
        dr (ScalarList): differentials of the radius [nm]
        r (ScalarList): neighborhood radii in ascending order [nm]
        n (int): dimension of the space of the region of interest

    Ouput:
        gggg (ScalarListList): stacked g++, g-+, g+-, g-- values [1]
        dV (ScalarList): differentials of the neighborhood [nm^n]
        dKKKK (ScalarList): differentials of KKKK [nm^n]

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
    dV[0] = np.nan # mask first value because of zero or wrong gradient
    dV[-1] = np.nan # mask last value because of wrong gradient
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

    Implementation of the antisymmetrical and symmetrical functions.

    This function is incremental. It uses pre-calculated MMMM values
    to avoid unnecessary repetition of calculations common to spatial
    analysis functions

    Input:
        MMMM (ScalarListList): stacked M++, M-+, M+-, M-- values [1]
        dr (ScalarList): differentials of the radius [nm]
        cp (Scalar): number of dislocations with sense + [1]
        cm (Scalar): number of dislocations with sense - [1]

    Output:
        GaGs (ScalarListList): stacked values of Ga and Gs [nm^-1]
        dMMMM (ScalarList): differentials of MMMM [1]

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
    dr[0] = np.nan # mask first value because of zero or wrong gradient
    dr[-1] = np.nan # mask last value because of wrong gradient
    dMMMM = np.gradient(MMMM, axis=1) # differentials of MMMM
    dMppdr, dMmpdr, dMpmdr, dMmmdr = dMMMM/dr # derivatives of MMMM
    Ga = cp*(dMppdr-dMpmdr) + cm*(dMmmdr-dMmpdr) # antisymmetrical function
    Gs = cp*(dMppdr+dMpmdr) + cm*(dMmmdr+dMmpdr) # symmetrical function
    return np.stack((Ga, Gs)), dMMMM

@beartype
def calculate(
    q: Union[list, tuple],
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
        'GaGs': stacked values of Ga and Gs [nm^-1]
        'cp': number of dislocations with sense + [1]
        'cm': number of dislocations with sense - [1]
        'dp': density of dislocations with sense + [nm^-n]
        'dm': density of dislocations with sense - [nm^-n]
        'dV': differentials of the neighborhood volume [nm^n]
        'dKKKK': differentials of K++, K-+, K+-, K-- [nm^n]
        'dMMMM': differentials of M++, M-+, M+-, M-- [1]

    Input:
        q (list): name of the quantities to calculate
        o (Distribution|Sample): object to analyze
        r (ScalarList): neighborhood radii in ascending order [nm]
      **r2 (ScalarList): squared neighborhood radii [nm^2]
      **dr (ScalarList): differentials of the neighborhood radii [nm]
      **ec (str): edge consideration (default: 'NEC')

    Output:
        v (AnalysisOutput): values of the quantities in requested order

    Exemple:
        KKKK, GaGs = calculate(['KKKK', 'GaGs'], distribution, r)

    Complexity:
        O( complexity_of(MMMM_Pp_Pm) ) if o is a distribution
        O( complexity_of(MMMM_Pp_Pm) * len(o) ) if o is a sample
    """
    # optional parameters
    r2 = getkwa('r2', kwargs, ScalarList, np.square(r))
    dr = getkwa('dr', kwargs, ScalarList, np.gradient(r))
    ec = getkwa('ec', kwargs, str, 'NEC')
    endkwa(kwargs)
    # create an incremental analysis function to be applied to a distribution
    @beartype
    def calculate_on_distribution(
        d: sets.Distribution,
    ) -> AnalysisOutput:
        """
        Auxiliary function applied to a distribution only.
        """
        s = {} # results temporary storage
        # flags determining what to calculate
        clcgggg = 'gggg' in q or 'dp' in q or 'dm' in q or 'dKKKK' in q
        clcKKKK = clcgggg or 'KKKK' in q or 'cp' in q or 'cm' in q
        clcGaGs = 'GaGs' in q or 'dMMMM' in q
        # calculate
        s['MMMM'], s['cp'], s['cm'] = MMMM_cp_cm(
            d,
            r,
            r2,
            ec,
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
    **kwargs,
) -> None:
    """
    Export the figure showing functions K++, K-+, K+-, K--.

    Input:
        r (ScalarList): radius of the neighborhoods [nm]
        KKKK (ScalarListList): stacked K++, K-+, K+-, K-- values [nm^n]
      **expdir (str): export directory (default: '')
      **expfmt (str): export format (default: 'pdf')
      **expstm (str): export stem (default: 'KKKK')
      **figttl (str): figure title (default: "")
      **edgcon (str): edge consideration (default: "")

    Complexity:
        O( len(r) )
    """
    # optional parameters
    expdir = getkwa('expdir', kwargs, str, '')
    expfmt = getkwa('expfmt', kwargs, str, 'pdf')
    expstm = getkwa('expstm', kwargs, str, 'KKKK')
    figttl = getkwa('figttl', kwargs, str, "")
    edgcon = getkwa('edgcon', kwargs, str, "")
    endkwa(kwargs)
    # K of a uniform distribution
    r_compare = np.linspace(r[0], r[-1], 10)
    k_compare = np.pi*r_compare**2
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.95, bottom=0.1)
    fig.suptitle(figttl, fontsize=16)
    # ax1
    ax1.plot(r, KKKK[0], label=fr"$K_{{++}}^{{ {edgcon} }}(r)$")
    ax1.plot(r, KKKK[1], label=fr"$K_{{-+}}^{{ {edgcon} }}(r)$")
    ax1.plot(r_compare, k_compare, "*", label=r"$ \pi r^2 $", color='black')
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r"$r \ (nm)$")
    ax1.set_ylabel(r"$(nm^2)$")
    # ax2
    ax2.plot(r, KKKK[3], label=fr"$K_{{--}}^{{ {edgcon} }}(r)$")
    ax2.plot(r, KKKK[2], label=fr"$K_{{+-}}^{{ {edgcon} }}(r)$")
    ax2.plot(r_compare, k_compare, "*", label=r"$ \pi r^2 $", color='black')
    ax2.legend()
    ax2.grid()
    ax2.set_xlabel(r"$r \ (nm)$")
    ax2.set_ylabel(r"$(nm^2)$")
    ax2.text(
        1.05,
        0.5,
        f"lpa-input ({__version__})",
        rotation=90,
        fontfamily='monospace',
        ha='left',
        va='center',
        transform=ax2.transAxes,
    )
    # export
    plt.savefig(os.path.join(expdir, expstm+"."+expfmt), format=expfmt)
    plt.close('all')

@beartype
def plot_gggg(
    r: ScalarList,
    gggg: ScalarListList,
    **kwargs,
) -> None:
    """
    Export the figure showing functions g++, g-+, g+-, g--.

    Input:
        r (ScalarList): radius of the neighborhoods [nm]
        gggg (ScalarListList): stacked g++, g-+, g+-, g-- values [1]
      **expdir (str): export directory (default: '')
      **expfmt (str): export format (default: 'pdf')
      **expstm (str): export stem (default: 'gggg')
      **figttl (str): figure title (default: "")
      **edgcon (str): edge consideration (default: "")

    Complexity:
        O( len(r) )
    """
    # optional parameters
    expdir = getkwa('expdir', kwargs, str, '')
    expfmt = getkwa('expfmt', kwargs, str, 'pdf')
    expstm = getkwa('expstm', kwargs, str, 'gggg')
    figttl = getkwa('figttl', kwargs, str, "")
    edgcon = getkwa('edgcon', kwargs, str, "")
    endkwa(kwargs)
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.95, bottom=0.1)
    fig.suptitle(figttl, fontsize=16)
    masked = np.ma.masked_invalid(gggg)
    margin = max(np.ptp(masked)*0.05, 0.1)
    ymin = np.min(masked) - margin
    ymax = np.max(masked) + margin
    # ax1
    ax1.plot(r, gggg[0], label=fr"$g_{{++}}^{{ {edgcon} }}(r)$")
    ax1.plot(r, gggg[1], label=fr"$g_{{-+}}^{{ {edgcon} }}(r)$")
    ax1.hlines(1, r[0], r[-1], label=r"$1$", color='black')
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r"$r \ (nm)$")
    ax1.set_ylim(ymin, ymax)
    # ax2
    ax2.plot(r, gggg[3], label=fr"$g_{{--}}^{{ {edgcon} }}(r)$")
    ax2.plot(r, gggg[2], label=fr"$g_{{+-}}^{{ {edgcon} }}(r)$")
    ax2.hlines(1, r[0], r[-1], label=r"$1$", color='black')
    ax2.legend()
    ax2.grid()
    ax2.set_xlabel(r"$r \ (nm)$")
    ax2.set_ylim(ymin, ymax)
    ax2.text(
        1.05,
        0.5,
        f"lpa-input ({__version__})",
        rotation=90,
        fontfamily='monospace',
        ha='left',
        va='center',
        transform=ax2.transAxes,
    )
    # export
    plt.savefig(os.path.join(expdir, expstm+"."+expfmt), format=expfmt)
    plt.close('all')

@beartype
def plot_GaGs(
    r: ScalarList,
    GaGs: ScalarListList,
    **kwargs,
) -> None:
    """
    Export the figure showing functions Ga and Gs.

    Input:
        r (ScalarList): radius of the neighborhoods [nm]
        GaGs (ScalarListList): stacked values of Ga and Gs [m^-1]
      **expdir (str): export directory (default: '')
      **expfmt (str): export format (default: 'pdf')
      **expstm (str): export stem (default: 'GaGs')
      **figttl (str): figure title (default: "")
      **edgcon (str): edge consideration (default: "")

    Complexity:
        O( len(r) )
    """
    # optional parameters
    expdir = getkwa('expdir', kwargs, str, '')
    expfmt = getkwa('expfmt', kwargs, str, 'pdf')
    expstm = getkwa('expstm', kwargs, str, 'GaGs')
    figttl = getkwa('figttl', kwargs, str, "")
    edgcon = getkwa('edgcon', kwargs, str, "")
    endkwa(kwargs)
    # fig
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    fig.subplots_adjust(left=0.06, right=0.95, bottom=0.1)
    fig.suptitle(figttl, fontsize=16)
    masked = np.ma.masked_invalid(GaGs)
    # ax1
    marg1 = max(np.ptp(masked[0])*0.05, 0.1)
    ax1.plot(r, GaGs[0], label=fr"$G_A^{{ {edgcon} }}(r)$")
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel(r"$r \ (nm)$")
    ax1.set_ylabel(r"$(nm^{-1})$")
    ax1.set_ylim(np.min(masked[0])-marg1, np.max(masked[0])+marg1)
    # ax2
    marg2 = max(np.ptp(masked[1])*0.05, 0.1)
    ax2.plot(r, GaGs[1], label=fr"$G_S^{{ {edgcon} }}(r)$")
    ax2.legend()
    ax2.grid()
    ax2.set_xlabel(r"$r \ (nm)$")
    ax2.set_ylabel(r"$(nm^{-1})$")
    ax2.set_ylim(np.min(masked[1])-marg2, np.max(masked[1])+marg2)
    ax2.text(
        1.05,
        0.5,
        f"lpa-input ({__version__})",
        rotation=90,
        fontfamily='monospace',
        ha='left',
        va='center',
        transform=ax2.transAxes,
    )
    # export
    plt.savefig(os.path.join(expdir, expstm+"."+expfmt), format=expfmt)
    plt.close('all')

@beartype
def export(
    o: Union[sets.Distribution, sets.Sample],
    **kwargs,
) -> tuple:
    """
    Export a complete analysis of the object o.

    Input:
        o (Distribution|Sample): object to analyze
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

    Complexity:
        O( complexity_of(calculate) )
    """
    # optional parameters
    fstm, fttl = 'dmgscS', 'mgscd'
    if isinstance(o, sets.Sample):
        fstm, fttl = "n"+fstm, "n"+fttl
    rmax = o.s
    if o.g == 'circle':
        rmax *= 2
    expstm = getkwa('expstm', kwargs, str, o.name(fstm, c='stm'))
    figttl = getkwa('figttl', kwargs, str, o.name(fttl, c='ttl'))
    edgcon = getkwa('edgcon', kwargs, str, 'NEC')
    intrad = getkwa('intrad', kwargs, ScalarList, np.linspace(0, rmax, 200))
    savtxt = getkwa('savtxt', kwargs, bool, False)
    kwargs['figttl'] = figttl
    kwargs['edgcon'] = edgcon
    # export
    fun = ('KKKK', 'gggg', 'GaGs') # functions to calculate
    KKKK, gggg, GaGs = calculate(fun, o, intrad, ec=edgcon) # function results
    KKKKstm = expstm+"_KKKK_"+edgcon
    ggggstm = expstm+"_gggg_"+edgcon
    GaGsstm = expstm+"_GaGs_"+edgcon
    plot_KKKK(intrad, KKKK, **kwargs, expstm=KKKKstm)
    plot_gggg(intrad, gggg, **kwargs, expstm=ggggstm)
    plot_GaGs(intrad, GaGs, **kwargs, expstm=GaGsstm)
    if savtxt:
        pth = kwargs.get('expdir', '')
        np.savetxt(os.path.join(pth, expstm+"_radii.txt"), intrad)
        np.savetxt(os.path.join(pth, KKKKstm+'.txt'), KKKK)
        np.savetxt(os.path.join(pth, ggggstm+'.txt'), gggg)
        np.savetxt(os.path.join(pth, GaGsstm+'.txt'), GaGs)
    return intrad, KKKK, gggg, GaGs
