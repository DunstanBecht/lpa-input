#!/usr/bin/env python
# coding: utf-8

"""
Tools for representing distributions and samples of distributions.
"""

from . import *
from . import overlap
from . import notation
from . import boundaries
from . import geometries

class Distribution:
    """
    Represent a distribution of dislocations.

    When the geometry is of type 'circle' the characteristic size s is
    its radius. When it is of type 'square' s is its side.

    The dislocations can be of type 'screw' or type 'edge'. If the type
    is not specified, type 'screw' is chosen by default.

    For circle geometries, the parameter c can be set to 'ISD'. For
    square geometries, it can be set to 'PBC<r>' or 'GBB<r>' where <r>
    is the rank of replications of the region of interest beyond the
    boundaries. The dislocations outside the area of ​​interest are then
    added to the distribution.

    Attributes:
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]
        n (int): dimension of space in the region of interest
        v (Scalar): n-volume of the region of interest [nm^n]
        m (GenerationFunction): distribution model generation function
        r (dict): distribution model parameters
        t (str): dislocation type
        p (VectorList): dislocation positions [nm]
        b (ScalarList): dislocation Burgers vector senses [1]
        d (Scalar): dislocation density [nm^-n]
        i (Scalar): inter dislocation distance [nm]
        c (NoneType|str): boundary conditions
        S (NoneType|int): random seed
        G (np.random._generator.Generator): random number generator
    """

    @beartype
    def __init__(self,
        g: str,
        s: Scalar,
        m: GenerationFunction,
        r: dict,
        t: str = 'screw',
        c: Optional[str] = None,
        S: Optional[int] = None,
        G: Optional[np.random._generator.Generator] = None,
    ) -> None:
        """
        Initialize the distribution.

        Input:
            g (str): geometry of the region of interest
            s (Scalar): size of the region of interest [nm]
            m (GenerationFunction): distribution model generation function
            r (dict): distribution model parameters
            t (str): dislocation type
            c (NoneType|str): boundary conditions
            S (NoneType|int): random seed
            G (NoneType|np.random._generator.Generator): random gen.

        Complexity:
            O( complexity_of(m) )
        """
        # shape
        if s <= 0:
            raise ValueError(f"incorrect size: {s}")
        self.s = s # characteristic size of the region of interest [nm]
        self.g = g # geometry of the region of interest
        self.n, self.v = geometries.nvolume(g, s) # dimension and n-volume
        # model
        self.m = m # model generation function
        self.r = r.copy() # model parameters
        # dislocations
        self.t = t # dislocation type
        if not G is None:
            self.G = G
            self.S = None
        else:
            self.G = np.random.default_rng(S) # random generator
            self.S = S # random seed
        self.p, self.b = self.m(self.g, self.s, self.v, self.r, self.G)
        self.d = len(self)/self.v # density of dislocations [nm^-n]
        if self.d == 0:
            raise ValueError("the model gives a density equal to 0")
        self.i = 1/np.sqrt(self.d) # inter dislocation distance [nm]
        # boundary conditions
        if c:
            cp, cb = self.conditions(c)
            self.p = np.concatenate((self.p, cp))
            self.b = np.concatenate((self.b, cb))
        self.c = c # boundary conditions name

    @beartype
    def __repr__(self) -> str:
        """Return the eval-able representation of the distribution."""
        args = (
            repr(self.g), # geometry
            str(self.s), # size
            self.m.__qualname__, # model generation function
            repr(self.r), # model parameters
            repr(self.t), # dislocation type
            repr(self.c), # bondary conditions
            repr(self.S), # random seed
        )
        return f"{self.__class__.__name__}({', '.join(args)})"

    @beartype
    def __str__(self) -> str:
        """Return a regular description of the distribution."""
        return self.name()

    @beartype
    def __len__(self) -> int:
        """Return the number of dislocations in the distribution."""
        return len(self.b)

    @beartype
    def conditions(self,
        c: str,
    ) -> tuple:
        """
        Return the dislocations for boundary conditions c.

        Input:
            c (str): boundary conditions name

        Output:
            cp (VectorList): outer dislocation positions [nm]
            cb (ScalarList): outer dislocation Burgers vectors sense [1]
        """
        if hasattr(self, 'c') and self.c:
            raise ValueError("boundary conditions already applied")
        if c=='ISD' and self.g=='circle' and self.t=='screw':
            cp = boundaries.image_positions(self.s, self.p)
            cb = - self.b
        elif 'PBC' in c and self.g=='square':
            u = boundaries.replication_displacements(int(c[3:]), self.s)
            cp = np.concatenate([self.p + u[i] for i in range(len(u))])
            cb = np.tile(self.b, len(u))
        elif 'GBB' in c and self.g=='square':
            u = boundaries.replication_displacements(int(c[3:]), self.s)
            cp = np.empty((0,2), dtype=self.p.dtype)
            cb = np.array([], dtype=self.b.dtype)
            for i in range(len(u)):
                p, b = self.m(self.g, self.s, self.v, self.r, self.G)
                cp = np.concatenate((cp, p + u[i]))
                cb = np.concatenate((cb, b))
        else:
            raise Exception(f"invalid boundary conditions: {c}")
        return cp, cb

    @beartype
    def name(self,
        f: str = 'dgsmtcS',
        c: str = 'csl',
    ) -> str:
        """
        Return a string describing the distribution in context c.

        Input:
            f (str): format
            c (str): context

        Output:
            csl|stm|ttl (str): string describing the distribution

        The format is defined by an arrangement of the characters:
            d: dislocation density [m^-2]
            g: geometry of the region of interest
            s: size of the region of interest [nm]
            m: model and its parameters
            t: dislocation type
            c: boundary conditions
            n: number of dislocations
            S: random seed
        """
        rho = notation.number(self.d*1e9**self.n, c, 7)
        rho += notation.unit(f"m^{{-{self.n}}}", c)
        v = {
            'd': notation.equality(r"\rho", rho, c),
            'g': self.g,
            's': notation.quantity(self.s, r"nm", c, w=4),
            'm': self.m.__name__+notation.parameters(self.r, c),
            't': self.t,
            'c': self.c,
            'n': str(int(len(self))),
            'S': notation.equality(r"S", str(self.S), c)
                if not self.S is None else None,
        }
        return notation.fmt(f, v, c)

    @beartype
    def w(self,
        a: Vector,
        r: ScalarList,
        r2: ScalarList,
    ) -> ScalarList:
        """
        Return the edge correction coefficients of the shape.

        Input:
            a (Vector): position of the neighborhoods center [nm]
            r (ScalarList): radius of the neighborhoods [nm]
            r2 (ScalarList): squared radius of neighborhoods [nm^2]

        Output:
            w (ScalarList): weighting coefficients for each radius

        Input example:
            a = np.array([x, y])
            r = np.array([r_0, r_1, r_2, ...])
            r2 = np.array([r_0^2, r_1^2, r_2^2, ...])

        Output example:
            w = np.array([w_r_0, w_r_1, w_r_2, ...])

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
        w = np.divide(vo, vv, np.ones(vo.size), where=r2>0) # ratio
        return w

class Sample:
    """
    Represent a sample of distributions.

    Most of the attributes have a common meaning with those of class
    Distribution. The density and the inter dislocation distance are
    averaged over all the generated distributions.

    Attributes:
        l (tuple): generated distributions
        g (str): geometry of the region of interest
        s (Scalar): size of the region of interest [nm]
        n (int): dimension of space in the region of interest
        v (Scalar): n-volume of the region of interest [nm^n]
        m (GenerationFunction): distribution model generation function
        r (dict): distribution model parameters
        t (str): dislocation type
        d (Scalar): averaged dislocation density [nm^-n]
        i (Scalar): averaged inter dislocation distance [nm]
        c (NoneType|str): boundary conditions
        S (NoneType|int): random seed
        G (np.random._generator.Generator): random generator
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
        S: Optional[int] = None,
    ) -> None:
        """
        Initialize the sample of distributions.

        Input:
            n (int): number of distributions to generate
            g (str): geometry of the region of interest
            s (Scalar): size of the region of interest [nm]
            m (GenerationFunction): model generation function
            r (dict): model parameters
            t (str): dislocation type
            c (NoneType|str): boundary conditions
            S (NoneType|str): random seed

        Complexity:
            O( c * complexity_of(Distribution) )
        """

        self.G = np.random.default_rng(S) # random generator
        self.S = S # random seed
        args = (g, s, m, r, t, c, None, self.G) # distribution arguments
        if n <= 0:
            raise ValueError(f"incorrect number of distribution: {n}")
        self.l = tuple([Distribution(*args) for i in range(n)])
        self.g = g # geometry of the region of interest
        self.s = s # size of the region of interest [nm]
        self.n = self[0].n # dimension of space
        self.v = self[0].v # n-volume of the region of interest [nm^n]
        self.m = m # model function
        self.r = r # model parameters
        self.t = t # dislocation type
        self.d = self.average(lambda d : d.d) # averaged dislocation density
        self.i = self.average(lambda d : d.i) # averaged inter disl. dist.
        self.c = c # boundary conditions

    @beartype
    def __repr__(self) -> str:
        """Return the eval-able representation of the sample."""
        args = (
            str(len(self)), # number of distributions
            repr(self.g), # geometry
            str(self.s), # size
            self.m.__qualname__, # model function
            repr(self.r), # model parameters
            repr(self.t), # dislocation type
            repr(self.c), # boundary conditions
            repr(self.S), # random seed
        )
        return f"{self.__class__.__name__}({', '.join(args)})"

    @beartype
    def __str__(self) -> str:
        """Return a regular description of the sample."""
        return self.name()

    @beartype
    def __getitem__(self, k: int) -> Distribution:
        """Return the k-th distribution stored in the sample."""
        return self.l[k]

    @beartype
    def __len__(self) -> int:
        """Return the number of distributions stored in the sample."""
        return len(self.l)

    @beartype
    def name(self,
        f: str = 'ndgsmtcS',
        c: str = 'csl',
    ) -> str:
        """
        Return a string describing the sample in context c.

        Input:
            f (str): format
            c (str): context

        Output:
            csl|stm|ttl (str): string describing the sample

        The format is defined by an arrangement of the characters:
            d: dislocation density [m^-2]
            g: geometry of the region of interest
            s: size of the region of interest [nm]
            m: model and its parameters
            t: dislocation type
            c: boundary conditions
            n: number of distributions generated
            S: random seed
        """
        rho = notation.number(self.d*1e9**self.n, c, 7)
        rho += notation.unit(f"m^{{-{self.n}}}", c)
        v = {
            'd': notation.equality(r"\rho", rho, c),
            'g': self.g,
            's': notation.quantity(self.s, r"nm", c, w=4),
            'm': self.m.__name__+notation.parameters(self.r, c),
            't': self.t,
            'c': self.c,
            'n': str(int(len(self))),
            'S': notation.equality(r"S", str(self.S), c)
                if not self.S is None else None,
        }
        return notation.fmt(f, v, c)

    @beartype
    def average(self,
        f: Callable[[Distribution, Any], AnalysisOutput],
        *args,
    ) -> AnalysisOutput:
        """
        Average the function f over the set of generated distributions.

        The first parameter of f must be the distribution to analyze.
        The function f must return a summable object or a tuple of
        summable objects. In the case of tuple, the items are averaged
        one by one according to their position in the tuple.

        Input:
            f (Callable): function applied to a distribution
            *args: additional arguments to pass to function f

        Output:
            r (AnalysisOutput): result of f averaged over the distrib.

        Complexity:
            O( len(self) * complexity_of(f) )
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
