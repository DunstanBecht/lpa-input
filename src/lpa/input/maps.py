#!/usr/bin/env python
# coding: utf-8

"""
Tools for exporting the dislocation map of a distribution.
"""

import matplotlib.pyplot as plt
from . import *
from . import __version__
from . import sets
from . import models

@beartype
def export(
    d: sets.Distribution,
    **kwargs,
) -> None:
    """
    Export the dislocation map of the distribution d.

    Input:
        d (Distribution): distribution to be exported
      **expdir (str): export directory
      **expfmt (str): export format
      **expstm (str): export stem
      **supttl (str): map sup title
      **subttl (str): map sub title

    Complexity:
        O( len(d) )
    """
    # optional parameters
    expdir = getkwa('expdir', kwargs, str, '')
    expfmt = getkwa('expfmt', kwargs, str, 'pdf')
    expstm = getkwa('expstm', kwargs, str, d.name('dgsmcS', 'stm'))
    supttl = getkwa('supttl', kwargs, str, d.name('dgsc', 'ttl'))
    subttl = getkwa('subttl', kwargs, str, d.name('m', 'ttl'))
    endkwa(kwargs)
    # fig
    fig, ax = plt.subplots(figsize=(6, 6))
    # aspect
    ax.set_aspect(1)
    b = d.s * 0.05 # borders width
    s, w = 100, 0.2 # marker size and line width
    # grid
    if d.m in (models.RRDD, models.RCDD) and d.c is None:
        ax.grid(True, zorder=0) # subareas or cells grid
        ticks = models.ticks(d.g, d.s, d.r['s'])
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        if len(ticks) > 10: # the grid is too thin to display all the ticks
            labels = ["" for i in range(len(ticks))]
            labels[0], labels[-1] = round(ticks[0]), round(ticks[-1])
            ax.set_xticklabels(labels)
            ax.set_yticklabels(labels)
    # region of interest
    if d.g == 'circle':
        g = plt.Circle((0,0), d.s, color='k', fill=False, zorder=50)
        if d.c is None:
            ax.set_xlim([-d.s-b, d.s+b])
            ax.set_ylim([-d.s-b, d.s+b])
        elif d.c == 'ISD':
            r = 1.5
            ax.set_xlim([-b-(r+0.5)*2*d.s, b+(r+0.5)*2*d.s])
            ax.set_ylim([-b-(r+0.5)*2*d.s, b+(r+0.5)*2*d.s])
            k = 2*r + 1
            w /= k
            s /= k**2
    else:
        g = plt.Rectangle((0,0), d.s, d.s, color='k', fill=False, zorder=50)
        if not d.c:
            ax.set_xlim([-b, d.s+b])
            ax.set_ylim([-b, d.s+b])
        elif 'PBC' in d.c or 'GBB' in d.c:
            r = int(d.c[3:])
            ax.set_xlim([-b-r*d.s, d.s+b+r*d.s])
            ax.set_ylim([-b-r*d.s, d.s+b+r*d.s])
            k = 2*r + 1
            w /= k
            s /= k**2
    ax.add_artist(g)
    # dislocations
    partition = [ # positive and negative Burgers vector senses
        (r"$ \bot $", "Burgers vector sense $+$", d.p[d.b>0]),
        (r"$ \top $", "Burgers vector sense $-$", d.p[d.b<0]),
    ]
    for k in partition:
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
    # information
    ax.set_xlabel(r"$x \ (nm)$")
    ax.set_ylabel(r"$y \ (nm)$")
    plt.suptitle(supttl)
    plt.title(subttl)
    l = plt.legend(facecolor='white', framealpha=1)
    l.set_zorder(150)
    ax.text(
        1.05,
        0.5,
        f"lpa-input ({__version__})",
        rotation=90,
        fontfamily='monospace',
        ha='left',
        va='center',
        transform=ax.transAxes,
    )
    # export
    plt.savefig(os.path.join(expdir, expstm+"."+expfmt), format=expfmt)
    plt.close('all')
