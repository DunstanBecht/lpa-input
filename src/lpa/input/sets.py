#!/usr/bin/env python
# coding: utf-8

"""
Tools for managing distributions and samples of distributions.
"""

from . import *
from . import models
from . import overlap

@beartype
def images(
    s: Scalar,
    p: VectorList,
    b: ScalarList,
) -> Tuple[VectorList, ScalarList]:
    """
    Return the image dislocations for a circle of radius s.

    Input:
        s: radius of the region of interest [nm]
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sense [1]

    Output:
        cp: image dislocation positions [nm]
        cb: image dislocation Burgers vectors sense [1]

    Complexity:
        O( len(p) )
    """
    n = np.linalg.norm(p, axis=1)
    m = n != 0
    v = p[m]
    t = np.arctan2(v[:,1], v[:,0])
    r = s**2/n[m]
    return np.stack((r*np.cos(t), r*np.sin(t)), axis=1), -b

@beartype
def replications(
    s: Scalar,
    p: VectorList,
    b: ScalarList,
    r: int,
) -> Tuple[VectorList, ScalarList]:
    """
    Return replications dislocations for a square of side s.

    Input:
        s: radius of the region of interest [nm]
        p: dislocation positions [nm]
        b: dislocation Burgers vectors sens [1]
        r: number of replications at the boundaries

    Output:
        cp: replicated dislocation positions [nm]
        cb: replicated dislocation Burgers vectors sense [1]

    Complexity:
        O( len(p) * r^2 )
    """
    u = []
    for i in range(1, r+1):
        for j in range(2*i):
            for k in [1, -1]:
                u.append(np.array((-i*k, (i-j)*k)))
                u.append(np.array(((i-j)*k,  i*k)))
    return np.concatenate([p + s*d for d in u]), np.tile(b, len(u))

class Distribution:
    """
    Represent a distribution of dislocations.

    When the geometry is of type 'circle' the characteristic size is
    its radius. When it is of type 'square' the characteristic size
    is its side.

    The dislocations can be of type 'screw' or type 'edge'. If the type
    is not specified, type 'screw' is chosen by default.

    For circle geometry, if the parameter c is set to 'idbc', image
    dislocations are added to the distribution.

    For square geometry, the parameter c can be set to 'pbcg<r>' or
    'pbcr<r>' where <r> is the number of replicates of the region of
    interest around the boundaries. With 'pbcg<r>' the replicated
    dislocations are added to the distribution. With 'pbcr<r>' the
    replicated dislocations are not added to the distribution but the
    X-ray diffraction simulation program is warned to replicate the
    region of interest.

    Attributes:
        g (str): geometry of the region of interest
        s (int): size of the region of interest [nm]
        n (int): dimension of space of the region of interest
        v (Scalar): n-dimensional volume of the region of interest [nm^n]
        m (GenerationFunction): function of the random model used
        r (dict): parameters of the random model
        t (str): dislocations type
        p (VectorList): dislocation positions [nm]
        b (ScalarList): dislocation Burgers vectors sense [1]
        d (Scalar): density of dislocations [nm^-n]
        i (Scalar): inter dislocation distance [nm]
        c (str|None): conditions at the boundaries
    """

    @beartype
    def __init__(self,
        g: str,
        s: int,
        m: str,
        r: dict,
        t: str = 'screw',
        c: Optional[str] = None,
    ) -> None:
        """
        Initialize the distribution.

        Input:
            g: geometry of the region of interest
            s: size of the region of interest [nm]
            m: name of the random model used to generate the dislocations
            r: parameters of the random model
            t: dislocations type
            c: conditions at the boundaries

        Complexity:
            O( complexity(m) )
        """
        # shape
        self.g = g
        if s <= 0:
            raise ValueError("s must be strictly positive")
        self.s = s # [nm]
        if self.g == 'circle': # centered at the origin
            self.n, self.v = 2, np.pi*self.s**2 # [nm^2]
        elif self.g == 'square': # corner at the origin
            self.n, self.v = 2, self.s**2 # [nm^2]
        else:
            raise Exception('unknown geometry')
        # dislocations
        self.m = eval('models.' + m)
        self.r = r.copy()
        self.t = t
        self.p, self.b = self.m(self.g, self.s, self.v, self.r)
        self.d = len(self)/self.v # [nm^-n]
        self.i = 1/np.sqrt(self.d) # [nm]
        self.c = c
        if not (self.c is None or c[:4]=='pbcr' and self.g=='square'):
            if self.g=='circle' and c=='idbc':
                cp, cb = images(self.s, self.p, self.b)
            elif self.g=='square' and c[:4]=='pbcg':
                cp, cb = replications(self.s, self.p, self.b, int(c[4:]))
            else:
                raise Exception('invalid boundary conditions')
            self.p = np.concatenate((self.p, cp))
            self.b = np.concatenate((self.b, cb))
        # printable attributes
        def fmt(x):
            return format(round(x), '1.0e').replace("+", "").replace("0", "")
        self.str_s = fmt(self.s) # [nm]
        self.str_n = str(self.n)
        self.str_v = fmt(self.v) # [nm^n]
        self.str_m = m
        self.str_d = fmt(self.d*1e9**self.n) # [m^-n]
        self.str_i = format(self.i, '1.1f') # [nm]

    @beartype
    def __repr__(self) -> str:
        """Return the representation of the distribution."""
        args = [self.g, self.s, self.str_m, self.r, self.t, self.c]
        return "Distribution("+", ".join([repr(a) for a in args])+")"

    @beartype
    def __str__(self) -> str:
        """Return a str version of the distribution."""
        r = ", ".join([k+"="+str(self.r[k]) for k in self.r])
        s = ("Distribution: "+self.fileName()
            + "\n- geometry: "+self.g
            + "\n- size: "+self.str_s+" nm"
            + "\n- n-volume: "+self.str_v+" nm^"+self.str_n
            + "\n- dislocations type: "+self.t
            + "\n- model: "+self.str_m+" ("+r+")"
            + "\n- population: "+str(len(self))+" dislocations"
            + "\n- dislocation density: "+self.str_d+" m^-"+self.str_n
            + "\n- inter dislocation distance: "+self.str_i+" nm"
            + "\n- boundary conditions: "+str(self.c)
            + "\n- b+: "+str(len(self.b[self.b>0]))+" dislocations"
            + "\n- b-: "+str(len(self.b[self.b<0]))+" dislocations")
        return s

    @beartype
    def __len__(self) -> int:
        """Return the number of dislocations in the distribution."""
        return len(self.b)

    @beartype
    def fileName(self,
        t: bool = True,
    ) -> str:
        """
        Return a name that can be used to export the distribution.

        Input:
            t: add dislocation type information

        Output:
            n: file name of the distribution
        """
        n = (self.g+"_"+self.str_s+"nm"+"_"
            + self.str_m+self.r['variant']+"_"
            + self.str_d+"m-"+self.str_n)
        if t:
            n += "_"+self.t
        if not self.c is None:
            n += "_"+self.c
        return n

    @beartype
    def plotTitle(self) -> str:
        """Return a name that can be used as a title for a plot."""
        s = "in a "+self.g
        m = self.str_m.upper()
        if 'variant' in self.r:
            m += self.r['variant'].upper()
        e = self.str_d.split("e") # mantissa and exponent of the density
        d = e[0]+"\\times 10^{"+str(int(e[1]))+"} m^{-2}" # density in LaTeX
        d = r"$ \left( \rho \sim "+d+r" \right) $"
        t = m+" "+d+" "+s
        if not self.c is None:
            t += " with "+self.c.upper()
        return t

    @beartype
    def w(self,
        a: Vector,
        r: ScalarList,
        r2: ScalarList,
    ) -> ScalarList:
        """
        Return the edge correction coefficients of the region of interest.

        Input:
            a: point around which neighborhoods are formed [nm]
            r: radius of the neighborhoods [nm]
            r2: squared radius of neighborhoods [nm^2]

        Output:
            w: weighting coefficients for each value of radius around a

        Complexity:
            O( len(r) )
        """
        if self.g == 'circle':
            d2 = np.sum(np.square(a))
            d = np.sqrt(d2)
            vo = overlap.circle_circle(r, self.s, d, r2, self.s**2, d2)
            vv = np.pi*r2
        elif self.g == 'square':
            vo = overlap.circle_square(a[0], a[1], r, r2, self.s)
            vv = np.pi*r2
        return 1/np.divide(vo, vv, np.ones(vo.size), where=r2>0)

class Sample:
    """
    Represent a sample of distributions.

    Except for l, the attributes have the same meaning as for the
    Distribution class. The density and the inter dislocation
    distance are averaged over all distributions.

    Attributes:
        l (Tuple[Distribution, ...]): sampled distributions
        g (str): geometry of the region of interest
        s (int): size of the region of interest [nm]
        n (int): dimension of space of the region of interest
        v (Scalar): n-dimensional volume of the region of interest [nm^n]
        m (str): name of the random model function
        r (dict): parameters of the random model
        t (str): dislocations type
        d (Scalar): averaged density of dislocations [nm^-n]
        i (Scalar): inter dislocation distance [nm]
        c (str|None): conditions at the boundaries
    """

    @beartype
    def __init__(self,
        n: int,
        g: str,
        s: int,
        m: str,
        r: dict,
        t: str = 'screw',
        c: Optional[str] = None,
    ) -> None:
        """
        Initialize the sample of distributions.

        Input:
            n: number of distributions to generate
            g: geometry of the region of interest
            s: size of the region of interest [nm]
            m: name of the random model function
            r: parameters of the random model
            t: dislocations type
            c: conditions at the boundaries

        Complexity:
            O( c * complexity(Distribution) )
        """
        if n <= 0:
            raise ValueError("n must be strictly positive")
        self.l = tuple([Distribution(g, s, m, r, t, c) for i in range(n)])
        self.g = g
        self.s = s
        self.n = self[0].n
        self.v = self[0].v
        self.m = m
        self.r = r
        self.t = t
        self.d = self.average(lambda d : d.d)
        self.i = self.average(lambda d : d.i)
        self.c = c
        def fmt(x):
            return format(round(x), '1.0e').replace("+", "")
        self.str_s = fmt(self.s) # [nm]
        self.str_n = str(self.n)
        self.str_v = fmt(self.v) # [nm^n]
        self.str_d = fmt(self.d*1e9**self.n) # [m^-n]
        self.str_i = format(self.i, '1.1f') # [nm]

    @beartype
    def __repr__(self) -> str:
        """Return the representation of the sample of distributions."""
        args = [len(self), self.g, self.s, self.m, self.r, self.t, self.c]
        return "Sample("+", ".join([repr(a) for a in args])+")"

    @beartype
    def __str__(self) -> str:
        """Return a str version of the sample of distributions."""
        r = ", ".join([k+"="+str(self.r[k]) for k in self.r])
        s = ("Sample: "+self.fileName()
            + "\n- population: "+str(len(self))+" distributions"
            + "\n- geometry: "+self.g
            + "\n- size: "+self.str_s+" nm"
            + "\n- n-volume: "+self.str_v+" nm^"+self.str_n
            + "\n- dislocations type: "+self.t
            + "\n- model: "+self.m+" ("+r+")"
            + "\n- dislocation density: "+self.str_d+" m^-"+self.str_n
            + "\n- inter dislocation distance: "+self.str_i+" nm"
            + "\n- boundary conditions: "+str(self.c))
        return s

    @beartype
    def __getitem__(self, k: int) -> Distribution:
        """Return the k th distribution of the sample."""
        return self.l[k]

    @beartype
    def __len__(self) -> int:
        """Return the number of distributions in the sample."""
        return len(self.l)

    @beartype
    def fileName(self, **kwargs) -> str:
        """Return a name that can be used to export the sample."""
        return str(len(self))+'_'+self[0].fileName(**kwargs)

    @beartype
    def plotTitle(self) -> str:
        """Return a name that can be used as a title for a plot."""
        return "sample of "+str(len(self))+" "+self[0].plotTitle()

    @beartype
    def average(self,
        f: Callable[[Distribution, Any], AnalysisOutput],
        *args,
    ) -> AnalysisOutput:
        """
        Average the function f over the sample of distributions.

        The first parameter of f must be the distribution to analyze.
        The function f must return a summable object or a tuple of
        summable objects. In the case of a tuple the items are averaged
        one by one according to their position in the tuple.

        Input:
            f: function that can be applied to a distribution
            *args: additional arguments to pass to function f

        Outpur:
            r: result of f applied to each distribution and averaged

        Complexity:
            O( len(self) * complexity(f) )
        """
        r = f(self[0], *args)
        if isinstance(r, tuple):
            r = list(r)
            n = len(r)
            for i in range(1, len(self)):
                ri = f(self[i], *args)
                for j in range(n):
                    r[j] += ri[j]
            for j in range(n):
                r[j] = r[j]/len(self)
            return tuple(r)
        else:
            for i in range(1, len(self)):
                r += f(self[i], *args)
            return r/len(self)
