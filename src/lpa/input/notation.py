#!/usr/bin/env python
# coding: utf-8

"""
Tools to standardize notations.

The following abbreviations define the contexts:
    stm: the string is used in the stem of a file name
    ttl: the string is a title and contains LaTeX code
    csl: the string is displayed in a console
"""

from . import *

@beartype
def number(
    x: Scalar,
    c: str = 'csl',
    w: int = 4,
) -> str:
    """
    Return a notation of the number x in context c.

    Input:
        x (Scalar): number
        c (str): context
        w (int): width of the output string

    Output:
        s (str): notation of x
    """
    if x<10**w and x>=1 or x == 0: # diplay x with fixed comma
        s = str(round(x))
        if c == 'stm':
            s = s.zfill(w) # fill with zeros
        elif c == 'csl':
            s = format(s, '>'+str(w)) # fill with spaces
    elif c == 'ttl': # LaTeX scientific notation
        s = format(x, f'1.{max(w-6, 0)}e')
        m, e = s.split("e") # separate mantissa and exponent
        s = fr"$ {m} \times 10^{{ {int(e)} }} $"
    elif c in ('stm', 'csl'): # scientific notation with 'e'
        s = format(x, '1.0e').replace("+", "")
        s = s.replace("e0", "e").replace("e-0", "e-")
        d = w-len(s) # number of missing characters to reach the ideal width
        if d>0:
            if c=='stm' or d==1:
                if "e-" in s:
                    s = s.replace("e-", "e-0")
                else:
                    s = s.replace("e", "e0")
            else:
                s = format(x, f'1.{d-1}e')
                s = s.replace("+", "").replace("e0", "e").replace("e-0", "e-")
    else:
        raise ValueError(f"unknown context: {c}")
    return s

@beartype
def unit(
    x: str,
    c: str = 'csl',
) -> str:
    """
    Return a notation of the unit x in context c.

    Input:
        x (str): unit
        c (str): context

    Output:
        s (str): notation of x
    """
    if c == 'ttl':
        return fr"$ \mathrm{{ {x} }} $"
    else:
        return x.replace("{", "").replace("}", "").replace("^", "")

@beartype
def quantity(
    v: Scalar,
    u: str,
    c: str = 'csl',
    w: int = 4,
) -> str:
    """
    Return the notation of value v and unit u in context c.

    Input:
        v (Scalar): value
        u (str): unit
        c (str): context
        w (int): width of the value

    Output:
        s (str): notation of u and v as a physical quantity
    """
    val = number(v, c, w)
    uni = unit(u, c)
    if c == 'ttl':
        val = val.replace(r"$", "").strip()
        uni = uni.replace(r"$", "").strip()
        return fr"$ {val} {uni} $"
    else:
        return val+uni

@beartype
def equality(
    a: str,
    b: str,
    c: str = 'csl',
) -> str:
    """
    Return the notation of a=b in context c.

    Input:
        a (str): left term of equality
        b (str): right term of equality
        c (str): context

    Output:
        s (str): notation of the equality
    """
    if c == 'ttl':
        a = a.replace(r"$", "").strip()
        b = b.replace(r"$", "").strip()
        return fr"$ {a} = {b} $"
    else:
        a = a.replace("\\", "")
        if c == 'csl':
            return a+"="+b
        elif c == 'stm':
            return a+b

@beartype
def parameters(
    r: dict,
    c: str = 'csl',
    s: bool = True,
) -> str:
    """
    Return the notation describing the model parameters set r.

    Input:
        r (dict): model parameters
        c (str): context
        s (bool): display the seed

    Output:
        n (str): model parameters description
    """
    # user-defined name for the set of parameters
    if 'name' in r:
        return r['name']
    n = ""
    # version
    if 'v' in r:
        n += "-"+r['v']
    # parameter list
    p = []
    if 'd' in r:  # density
        p.append(equality('d', quantity(r['d'], "nm^{-2}", c), c))
    if 's' in r: # subarea/cell side
        p.append(equality("s", number(r['s'], c, 4)+unit("nm", c), c))
    if 'f' in r: # filling
        p.append(equality("f", str(r['f']), c))
    if 't' in r: # wall thickness
        p.append(equality("t", number(r['t'], c, 3)+unit("nm", c), c))
    if 'l' in r: # dipole length
        p.append(equality("l", number(r['l'], c, 3)+unit("nm", c), c))
    # concatenate
    if len(p) > 0:
        if c == 'stm':
            n += "_"+"_".join(p)
        elif c == 'ttl':
            p = [q.replace("$", "").strip() for q in p]
            n += r" $ \left( "+r", \ ".join(p)+r" \right) $"
        elif c == 'csl':
            n += f" ({' '.join(p)})"
        else:
            raise ValueError(f"unknown context {c}")
    return n

@beartype
def fmt(
    f: str,
    v: dict,
    c: str = 'csl',
) -> str:
    """
    Return the elements of v selected and ordered by f.

    Input:
        f (str): format
        v (dict): dictionary containing the parameters
        c (str): context

    Output:
        s (str): elements of v selected and ordered by f
    """
    if c == 'stm':
        sep = '_'
    else:
        sep = " "
    return sep.join([v[n] for n in f if not v[n] is None])
