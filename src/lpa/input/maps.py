#!/usr/bin/env python
# coding: utf-8

"""
Tools for exporting the dislocation map of a distribution.
"""

import matplotlib.pyplot as plt
from . import *
from . import sets
from . import models

@beartype
def export(
    d: sets.Distribution,
    p: str = "",
    n: Optional[str] = None,
    t: Optional[str] = None,
) -> None:
    """
    Export the dislocation map of the distribution d.

    Input:
        d: distribution to be exported
        p: path where to export the file
        n: name of the exported file
        t: title of the plot

    Complexity:
        O( len(d) )
    """
    if p!="" and p[-1]!="/":
        p += "/"
    if n == None:
        n = d.fileName(t=False)
    if t == None:
        t = d.plotTitle()
    # export
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect(1)
    b = d.s * 0.05 # borders width
    s, w = 100, 0.2 # marker size and line width
    if d.str_m in ["rrdd", "rcdd"] and (d.c is None or 'pbcr' in d.c):
        ax.grid(True, zorder=0) # subareas grid
        ticks = models.ticks(d.g, d.s, d.r['s'])
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        if len(ticks) > 10:
            labels = ["" for i in range(len(ticks))]
            labels[0], labels[-1] = round(ticks[0]), round(ticks[-1])
            ax.set_xticklabels(labels)
            ax.set_yticklabels(labels)
    if d.g == 'circle':
        g = plt.Circle((0,0), d.s, color='k', fill=False, zorder=50)
        if d.c is None:
            ax.set_xlim([-d.s-b, d.s+b])
            ax.set_ylim([-d.s-b, d.s+b])
        elif d.c == 'idbc':
            r = 1.5
            ax.set_xlim([-b-(r+0.5)*2*d.s, b+(r+0.5)*2*d.s])
            ax.set_ylim([-b-(r+0.5)*2*d.s, b+(r+0.5)*2*d.s])
            k = 2*r + 1
            w /= k
            s /= k**2
    else:
        g = plt.Rectangle((0,0), d.s, d.s, color='k', fill=False, zorder=50)
        if d.c is None:
            ax.set_xlim([-b, d.s+b])
            ax.set_ylim([-b, d.s+b])
        elif 'pbcg' in d.c:
            r = int(d.c[4:])
            ax.set_xlim([-b-r*d.s, d.s+b+r*d.s])
            ax.set_ylim([-b-r*d.s, d.s+b+r*d.s])
            k = 2*r + 1
            w /= k
            s /= k**2
    ax.add_artist(g)
    ax.set_xlabel(r"$x \ (nm)$")
    ax.set_ylabel(r"$y \ (nm)$")
    plt.title(t)
    t = [ # partition of Burgers vectors set
        (r"$\bot$", "Burgers vector sense $+$", d.p[d.b>0]),
        (r"$\top$", "Burgers vector sense $-$", d.p[d.b<0]),
    ]
    for k in t:
        x, y = k[2][:,0], k[2][:,1]
        plt.scatter(
            x,
            y,
            marker=k[0],
            label=k[1],
            zorder=100,
            linewidths=w,
            s=s,
        )
    l = plt.legend(facecolor='white', framealpha=1)
    l.set_zorder(150)
    plt.savefig(p+n+".pdf")
    plt.close('all')
