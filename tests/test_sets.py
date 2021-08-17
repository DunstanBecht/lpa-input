#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module sets.
"""

from lpa.input import sets
from lpa.input.models import RDD, RRDD, RCDD

rdd = (RDD, {'d': 5e13*1e-18})
rrdd = (RRDD, {'v': 'E', 's': 200, 'f': 2})
rcdd = (RCDD, {'v': 'R', 'd': 1e15*1e-18, 's': 200, 't': 10})

# instantiate
d_rdd = sets.Distribution('circle', 1000, *rdd, S=0)
d_rrdd = sets.Distribution('square', 2000, *rrdd, S=1)
d_rcdd = sets.Distribution('square', 2000, *rcdd, c='PBCR3', S=2)
s_rdd = sets.Sample(10, 'square', 2000, *rdd, t='edge', S=3)
s_rrdd = sets.Sample(10, 'circle', 1000, *rrdd, c='IDBC', S=4)
s_rcdd = sets.Sample(10, 'square', 2000, *rcdd, c='PBCG1', S=5)

if __name__ == "__main__":

    # get stem string
    print(d_rcdd.name(c='stm', f='dgsmtcS'))
    print(d_rcdd.name(c='stm', f='dgsmtc'))
    print(d_rcdd.name(c='stm', f='dgsmt'))
    print(d_rcdd.name(c='stm', f='dgsm'))
    print(d_rcdd.name(c='stm', f='dgs'))
    print(d_rcdd.name(c='stm', f='dg'))
    print(d_rcdd.name(c='stm', f='d'), end="\n\n")

    # get title string
    print(d_rdd.name(c='ttl', f='mtcS'))
    print(d_rdd.name(c='ttl', f='mtc'))
    print(d_rdd.name(c='ttl', f='mt'))
    print(d_rdd.name(c='ttl', f='dgs'))
    print(d_rdd.name(c='ttl', f='dg'))
    print(d_rdd.name(c='ttl', f='d'), end="\n\n")

    # get console string
    print(s_rcdd.name(c='csl', f='ndgsmtcS'))
    print(s_rcdd.name(c='csl', f='ndgsmtc'))
    print(s_rcdd.name(c='csl', f='ndgsmt'))
    print(s_rcdd.name(c='csl', f='ndgsm'))
    print(s_rcdd.name(c='csl', f='ndgs'))
    print(s_rcdd.name(c='csl', f='ndg'))
    print(s_rcdd.name(c='csl', f='nd'), end="\n\n")

    # get representation
    print(repr(d_rdd))
    print(repr(d_rrdd))
    print(repr(d_rcdd))
    print(repr(s_rdd))
    print(repr(s_rrdd))
    print(repr(s_rcdd), end="\n\n")

    # evaluate a representation
    eval('sets.'+repr(d_rdd))
    eval('sets.'+repr(d_rrdd))
    eval('sets.'+repr(d_rcdd))
    eval('sets.'+repr(s_rdd))
    eval('sets.'+repr(s_rrdd))
    eval('sets.'+repr(s_rcdd))

    # print the object
    print(d_rdd)
    print(d_rrdd)
    print(d_rcdd)
    print(s_rdd)
    print(s_rrdd)
    print(s_rcdd, end="\n\n")

    # average over a sample
    rho = lambda dist: len(dist)/dist.v
    print(str(s_rdd.average(rho)))
    print(str(s_rrdd.average(rho)))
    print(str(s_rcdd.average(rho)), end="\n\n")

    input("OK")
