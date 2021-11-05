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
    S = 0 if x>=0 else 1 # 0 for + / 1 for -
    m, e = f"{x:e}".split('e')
    m, e = eval(m), int(e) # mantissa and exponent
    A, B = f"{x:f}".rstrip('0').strip('-').split('.')
    a, b = len(A), len(B) # number of digits before / after the comma
    if isinstance(x, int) or x%1==0:
        if a+S<=w or a<=3:
            if c == 'stm':
                return f"{S*'-'+A.zfill(w-S)}"
            if c == 'ttl':
                return f"$ {S*'-'+A} $"
            if c == 'csl':
                return format(S*'-'+A, f'>{w}')
    else:
        if a+S<=w-2 and np.abs(x)>0.1:
            if c == 'ttl':
                return f"$ {S*'-'}{A}.{B[:w-S-1]} $"
            if c == 'csl':
                return format(f"{S*'-'}{A}.{B[:w-S-1-a]}", f'>{w}')
    u = len(str(e))
    if c == 'stm':
        z = 0
        q = len(str(m).strip('0').strip('.'))-2-S
        while z+1<=q and w>=S+1+(z+1)+1+len(str(e-(z+1))):
            z += 1
        return f"{m*10**z:1.0f}e{e-z}"
    if c == 'ttl':
        return fr"$ {format(m, f'1.{max(0, w-S-4-u)}f')} \times 10^{{ {e} }} $"
    if c == 'csl':
        return format(f"{format(m, f'1.{max(0, w-S-3-u)}f')}e{e}", f'>{w}')
    raise ValueError(f"unknown context: {c}")

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
