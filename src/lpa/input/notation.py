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
        x: number
        c: context
        w: width of the output string

    Output:
        s: notation of x
    """
    if x<10**w and x>=1 or x == 0: # diplay x with fixed comma
        s = str(round(x))
        if c == 'stm':
            s = s.zfill(w) # fill with zeros
        elif c == 'csl':
            s = format(s, '>'+str(w)) # fill with spaces
    elif c == 'ttl': # LaTeX scientific notation
        s = format(x, '1.'+str(max(w-6, 0))+'e')
        m, e = s.split("e") # separate mantissa and exponent
        s = "$ "+m+" \\times 10^{"+str(int(e))+"} $"
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
                fmt = '1.'+str(d-1)+'e'
                s = format(x, fmt)
                s = s.replace("+", "").replace("e0", "e").replace("e-0", "e-")
    else:
        raise ValueError("unknown context: "+str(c))
    return s

@beartype
def unit(
    x: str,
    c: str = 'csl',
) -> str:
    """
    Return a notation of the unit x in context c.

    Input:
        x: unit
        c: context

    Output:
        s: notation of x
    """
    if c == 'ttl':
        return "$ \mathrm{"+x+"} $"
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
        v: value
        u: unit
        c: context
        w: width of the value

    Output:
        s: notation of u and v as a physical quantity
    """
    val = number(v, c, w)
    uni = unit(u, c)
    if c == 'ttl':
        val = val.replace(r"$", "").strip()
        uni = uni.replace(r"$", "").strip()
        return "$ "+val+" "+uni+" $"
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
        a: left term of equality
        b: right term of equality
        c: context

    Output:
        s: notation of the equality
    """
    if c == 'ttl':
        a = a.replace(r"$", "").strip()
        b = b.replace(r"$", "").strip()
        return r"$ "+a+" = "+b+" $"
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
    Return the notation descibing the model parameters set r.

    Input:
        r: model parameters
        c: context
        s: display the seed

    Output:
        s: model parameters description
    """
    # user-defined name for the set of parameters
    if 'name' in r:
        return r['name']
    s = ""
    # version
    if 'v' in r:
        s += "-"+r['v']
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
    if 'r' in r and s: # random seed
        p.append(equality("r", str(r['r']), c))
    # concatenate
    if len(p) > 0:
        if c == 'stm':
            s += "_"+"_".join(p)
        elif c == 'ttl':
            p = [q.replace("$", "").strip() for q in p]
            s += r" $ \left( "+", \ ".join(p)+r" \right) $ "
        elif c == 'csl':
            s += " ("+" ".join(p)+")"
        else:
            raise ValueError("unknown context "+str(c))
    return s

@beartype
def fmt(
    f: str,
    v: dict,
    c: str = 'csl',
) -> str:
    """
    Return the elements of v selected and ordered by f.

    Input:
        f: format
        v: dictionary containing the parameters
        c: context

    Output:
        s: elements of v selected and ordered by f
    """
    if c == 'stm':
        sep = '_'
    else:
        sep = " "
    return sep.join([v[n] for n in f if not v[n] is None])
