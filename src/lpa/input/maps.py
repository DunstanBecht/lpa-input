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
    **kwargs,
) -> None:
    """
    Export the dislocation map of the distribution d.

    Input:
        d: distribution to be exported
      **expdir: export directory
      **expfmt: export format
      **expstm: export stem
      **supttl: map sup title
      **subttl: map sub title

    Complexity:
        O( len(d) )
    """
    # optional parameters
    expdir = kwargs.pop('expdir', '') # export directory
    expfmt = kwargs.pop('expfmt', 'pdf') # export format
    expstm = kwargs.pop('expstm', d.name('dgsmcS', 'stm')) # export stem
    supttl = kwargs.pop('supttl', d.name('dgsc', 'ttl')) # sup title
    subttl = kwargs.pop('subttl', d.name('m', 'ttl')) # sub title
    if len(kwargs)>0:
        raise ValueError("wrong keywords: "+str(kwargs))
    # fig
    fig, ax = plt.subplots(figsize=(6, 6))
    # aspect
    ax.set_aspect(1)
    b = d.s * 0.05 # borders width
    s, w = 100, 0.2 # marker size and line width
    # grid
    if d.m in (models.RRDD, models.RCDD) and (d.c is None or d.c[:4]=='PBCR'):
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
        elif d.c == 'IDBC':
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
        elif d.c[:4]=='PBCG':
            r = int(d.c[4:])
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
    # export
    plt.savefig(os.path.join(expdir, expstm+"."+expfmt), format=expfmt)
    plt.close('all')
