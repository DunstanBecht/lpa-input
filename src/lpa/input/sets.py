#!/usr/bin/env python
# coding: utf-8

"""
Tools for representing distributions and samples of distributions.
"""

from . import *
from . import models
from . import overlap
from . import notation
from . import boundaries

class Distribution:
    """
    Represent a distribution of dislocations.

    When the geometry is of type 'circle' the characteristic size is
    its radius. When it is of type 'square' the characteristic size is
    its side.

    The dislocations can be of type 'screw' or type 'edge'. If the type
    is not specified, type 'screw' is chosen by default.

    For circle geometry, if the parameter c is set to 'IDBC', image
    dislocations are added to the distribution.

    For square geometry, the parameter c can be set to 'PBCG<r>' or
    'PBCR<r>' where <r> is the number of replicates of the region of
    interest around the boundaries. With 'PBCG<r>' the replicated
    dislocations are added to the distribution. With 'PBCR<r>' the
    replicated dislocations are not added to the distribution but the
    X-ray diffraction simulation program is warned to replicate the
    region of interest.

    Attributes:
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]
        n (int): dimension of space of the region of interest
        v (Scalar): n-volume of the region of interest [nm^n]
        m (GenerationFunction): function of the random model used
        r (dict): parameters of the random model
        t (str): dislocations type
        p (VectorList): dislocation positions [nm]
        b (ScalarList): dislocation Burgers vectors sense [1]
        d (Scalar): density of dislocations [nm^-n]
        i (Scalar): inter dislocation distance [nm]
        c (str|None): conditions at the boundaries
    """

    geometries = ('square', 'circle') # available geometries

    @beartype
    def __init__(self,
        g: str,
        s: Scalar,
        m: GenerationFunction,
        r: dict,
        t: str = 'screw',
        c: Optional[str] = None,
    ) -> None:
        """
        Initialize the distribution.

        Input:
            g: geometry of the region of interest
            s: size of the region of interest [nm]
            m: model generating function
            r: parameters of the random model
            t: dislocations type
            c: conditions at the boundaries

        Complexity:
            O( complexity(m) )
        """
        # shape size
        if s <= 0:
            raise ValueError("incorrect size: "+str(s))
        self.s = s # characteristic size of the region of interest [nm]
        # shape geometry
        if not g in Distribution.geometries:
            raise ValueError("unknown geometry: "+str(g))
        self.g = g # geometry of the region of interest
        # shape n-volume
        if self.g == 'circle': # centered at the origin
            self.n = 2 # dimension of space
            self.v = np.pi*self.s**2 # surface [nm^2]
        elif self.g == 'square': # bottom left corner at the origin
            self.n = 2 # dimension of space
            self.v = self.s**2 # surface [nm^2]
        # model
        self.m = m # model generating function
        self.r = r.copy() # model parameters
        # dislocations
        self.t = t # dislocation type
        self.p, self.b = self.m(self.g, self.s, self.v, self.r) # generate
        self.d = len(self)/self.v # density of dislocations [nm^-n]
        if self.d == 0:
            raise ValueError("the model gives a density equal to 0")
        self.i = 1/np.sqrt(self.d) # inter dislocation distance [nm]
        # boundary conditions
        self.c = c # boundary conditions code
        if not (self.c is None or c[:4]=='PBCR' and self.g=='square'):
            if self.g=='circle' and c=='IDBC':
                cp, cb = boundaries.IDBC(self.s, self.p, self.b, self.t)
            elif self.g=='square' and c[:4]=='PBCG':
                cp, cb = boundaries.PBCG(self.s, self.p, self.b, int(c[4:]))
            else:
                raise Exception("invalid boundary conditions: "+str(c))
            self.p = np.concatenate((self.p, cp))
            self.b = np.concatenate((self.b, cb))

    @beartype
    def __repr__(self) -> str:
        """
        Return the eval-able representation of the distribution.

        Output:
            r: representation of the distribution
        """
        args = [
            repr(self.g), # geometry
            str(self.s), # size
            self.m.__name__, # model generating function
            repr(self.r), # model parameters
            repr(self.t), # dislocations type
            repr(self.c), # bondary conditions
        ]
        return "Distribution("+", ".join(args)+")"

    @beartype
    def __str__(self) -> str:
        """
        Return a textual description of the distribution.

        Output:
            s: distribution information
        """
        s = (
            "Distribution: "+self.stem()
            + "\n- geometry: "+self.g
            + "\n- size: "+notation.number(self.s)+" nm"
            + "\n- n-volume: "+notation.number(self.v)+" nm^"+str(self.n)
            + "\n- dislocations type: "+self.t
            + "\n- model: "+self.m.__name__+notation.parameters(self.r)
            + "\n- population: "+str(len(self))+" dislocations"
            + "\n- dislocation density: "+notation.number(self.d*1e9**self.n)
            + " m^-"+str(self.n)
            + "\n- inter dislocation distance: "+str(round(self.i))+" nm"
            + "\n- boundary conditions: "+str(self.c)
            + "\n- b+: "+str(len(self.b[self.b>0]))+" dislocations"
            + "\n- b-: "+str(len(self.b[self.b<0]))+" dislocations"
        )
        return s

    @beartype
    def __len__(self) -> int:
        """
        Return the number of dislocations in the distribution.

        Output:
            n: number of Burgers vectors
        """
        return len(self.b)

    @beartype
    def stem(self,
        t: bool = True,
        s: bool = True,
    ) -> str:
        """
        Return a string that can be used in a file name.

        Input:
            t: add dislocation type
            s: add seed

        Output:
            stm: name identifying the distribution
        """
        d = notation.number(self.d*1e9**self.n, 'stm') # density
        stm = (
            d+"m-"+str(self.n) # density
            + "_"+self.g # geometry
            + "_"+notation.number(self.s, 'stm', w=4)+"nm" # size
            + "_"+self.m.__name__ # model
            + notation.parameters(self.r, 'stm', s=s) # model parameters
        )
        if t:
            stm += "_"+self.t # add dislocation type information
        if self.c:
            stm += "_"+self.c # add bondary conditions information
        return stm

    @beartype
    def title(self,
        t: bool = True,
        s: bool = True,
    ) -> str:
        """
        Return a name that can be used as a title in a plot.

        Input:
            t: add dislocation type
            s: add seed

        Output:
            ttl: title containing LaTeX code
        """
        d = notation.number(self.d*1e9**self.n, 'ttl') # density
        ttl = (
            self.m.__name__ # model
            + notation.parameters(self.r, 'ttl', s=s) # model parameters
            + r" $ \rho \sim "+d+r" m^{-2} $" # density
        )
        if t:
            ttl += " "+self.t # add dislocation type information
        if self.c and self.c[:4]!='PBCR':
            ttl += " "+self.c # add bondary conditions information
        return ttl

    @beartype
    def w(self,
        a: Vector,
        r: ScalarList,
        r2: ScalarList,
    ) -> ScalarList:
        """
        Return the edge correction coefficients of the shape.

        Input:
            a: position around which neighborhoods are formed [nm]
            r: radius of the neighborhoods [nm]
            r2: squared radius of neighborhoods [nm^2]

        Output:
            w: weighting coefficients for each value of radius around a

        Input example:
            a = np.array([0, 0])
            r = np.array([1, 2, 3, ...])
            r2 = np.array([1, 4, 9, ...])

        Output example:
            w = np.array([1, 1, 0.9, ...])

        Complexity:
            O( len(r) )
        """
        # vo: n-volumes of the overlappings (for each value of r)
        # vv: n-volumes of the vicinities of a (for each value of r)
        if self.g == 'circle':
            d2 = np.sum(np.square(a)) # squared distance to the origin of a
            d = np.sqrt(d2) # distance to the origin of a
            vo = overlap.circle_circle(r, self.s, d, r2, self.s**2, d2)
            vv = np.pi*r2
        elif self.g == 'square':
            vo = overlap.circle_square(a[0], a[1], r, r2, self.s)
            vv = np.pi*r2
        w = 1/np.divide(vo, vv, np.ones(vo.size), where=r2>0) # ratio
        return w

class Sample:
    """
    Represent a sample of distributions.

    Except for l, the attributes have the same meaning as for the
    Distribution class. The density and the inter dislocation distance
    are averaged over all distributions.

    Attributes:
        l (Tuple[Distribution, ...]): sampled distributions
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]
        n (int): dimension of space of the region of interest
        v (Scalar): n-volume of the region of interest [nm^n]
        m (GenerationFunction): function of the random model used
        r (dict): parameters of the random model
        t (str): dislocations type
        d (Scalar): averaged density of dislocations [nm^-n]
        i (Scalar): averaged inter dislocation distance [nm]
        c (str|None): conditions at the boundaries
    """

    @beartype
    def __init__(self,
        n: int,
        g: str,
        s: Scalar,
        m: GenerationFunction,
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
            m: model generating function
            r: parameters of the random model
            t: dislocations type
            c: conditions at the boundaries

        Complexity:
            O( c * complexity(Distribution) )
        """
        if n <= 0:
            raise ValueError("incorrect number of distribution: "+str(n))
        self.l = tuple([Distribution(g, s, m, r, t, c) for i in range(n)])
        self.g = g # geometry of the region of interest
        self.s = s # size of the region of interest
        self.n = self[0].n # dimension of space
        self.v = self[0].v # n-volume of the region of interest
        self.m = m # model function
        self.r = r # model parameters
        self.t = t # dislocation type
        self.d = self.average(lambda d : d.d) # averaged dislocation density
        self.i = self.average(lambda d : d.i) # averaged inter disl. dist.
        self.c = c # boundary conditions

    @beartype
    def __repr__(self) -> str:
        """
        Return the eval-able representation of the sample of distributions.

        Output:
            r: representation of the sample of distributions
        """
        args = [
            str(len(self)), # number of distributions
            repr(self.g), # geometry
            str(self.s), # size
            self.m.__name__, # model function
            repr(self.r), # model parameters
            repr(self.t), # dislocation type
            repr(self.c), # boundary conditions
        ]
        return "Sample("+", ".join(args)+")"

    @beartype
    def __str__(self) -> str:
        """
        Return a textual description of the distribution.

        Output:
            s: distribution information
        """
        s = (
            "Sample: "+self.stem()
            + "\n- population: "+str(len(self))+" distributions"
            + "\n- geometry: "+self.g
            + "\n- size: "+notation.number(self.s)+" nm"
            + "\n- n-volume: "+notation.number(self.v)+" nm^"+str(self.n)
            + "\n- dislocations type: "+self.t
            + "\n- model: "+self.m.__name__+notation.parameters(self.r)
            + "\n- dislocation density: "+notation.number(self.d*1e9**self.n)
            + " m^-"+str(self.n)
            + "\n- inter dislocation distance: "+str(round(self.i))+" nm"
            + "\n- boundary conditions: "+str(self.c)
        )
        return s

    @beartype
    def __getitem__(self, k: int) -> Distribution:
        """
        Return the k-th distribution stored in the sample.

        Input:
            k: index of the distribution

        Output:
            d: k-th distribution
        """
        return self.l[k]

    @beartype
    def __len__(self) -> int:
        """
        Return the number of distributions stored in the sample.

        Output:
            n: number of stored distributions
        """
        return len(self.l)

    @beartype
    def stem(self, **kwargs) -> str:
        """
        Return a string that can be used in a file name.

        The arguments are passed to the equivalent function in the
        Distribution class.

        Output:
            stm: name identifying the sample
        """
        return str(len(self))+'_'+self[0].stem(**kwargs)

    @beartype
    def title(self, **kwargs) -> str:
        """
        Return a name that can be used as a title in a plot.

        Output:
            ttl: title containing LaTeX code
        """
        return str(len(self))+" "+self[0].title(**kwargs)

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

        Output:
            r: result of f applied to each distribution and averaged

        Complexity:
            O( len(self) * complexity(f) )
        """
        r = f(self[0], *args) # result of f on the first distribution
        if isinstance(r, tuple):
            r = list(r)
            n = len(r)
            for i in range(1, len(self)):
                ri = f(self[i], *args) # result of f on the ith distribution
                for j in range(n):
                    r[j] += ri[j]
            for j in range(n):
                r[j] = r[j]/len(self)
            return tuple(r)
        else:
            for i in range(1, len(self)):
                r += f(self[i], *args)
            return r/len(self)
