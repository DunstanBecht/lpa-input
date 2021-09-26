#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module sets.
"""

from lpa.input import sets
from test_models import *

rdd = (models.RDD, prm_rdd)
rrdd = (models.RRDD, prm_rrdd)
rcdd = (models.RCDD, prm_rcdd)

if __name__ == "__main__":

    # instantiate
    d_rdd = sets.Distribution('circle', 1000, *rdd, S=0)
    d_rrdd = sets.Distribution('square', 2000, *rrdd, S=0)
    d_rcdd = sets.Distribution('square', 2000, *rcdd, c='PBC3', S=0)
    s_rdd = sets.Sample(10, 'square', 2000, *rdd, t='edge', S=0)
    s_rrdd = sets.Sample(10, 'circle', 1000, *rrdd, c='ISD', S=0)
    s_rcdd = sets.Sample(10, 'square', 2000, *rcdd, c='GBB1', S=0)

    # get stem string
    print(d_rcdd.name(c='stm', f='dgsmtcS'))
    print(d_rcdd.name(c='stm', f='dgsmtc'))
    print(d_rcdd.name(c='stm', f='dgsmt'))
    print(d_rcdd.name(c='stm', f='dgsm'))
    print(d_rcdd.name(c='stm', f='dgs'))
    print(d_rcdd.name(c='stm', f='dg'))
    print(d_rcdd.name(c='stm', f='d'))
    print()

    # get title string
    print(d_rdd.name(c='ttl', f='mtcS'))
    print(d_rdd.name(c='ttl', f='mtc'))
    print(d_rdd.name(c='ttl', f='mt'))
    print(d_rdd.name(c='ttl', f='dgs'))
    print(d_rdd.name(c='ttl', f='dg'))
    print(d_rdd.name(c='ttl', f='d'))
    print()

    # get console string
    print(s_rcdd.name(c='csl', f='ndgsmtcS'))
    print(s_rcdd.name(c='csl', f='ndgsmtc'))
    print(s_rcdd.name(c='csl', f='ndgsmt'))
    print(s_rcdd.name(c='csl', f='ndgsm'))
    print(s_rcdd.name(c='csl', f='ndgs'))
    print(s_rcdd.name(c='csl', f='ndg'))
    print(s_rcdd.name(c='csl', f='nd'))
    print()

    # evaluate a representation
    print(repr(d_rdd))
    print(repr(eval('sets.'+repr(d_rdd).replace("RDD","models.RDD"))))
    print(repr(d_rrdd))
    print(repr(eval('sets.'+repr(d_rrdd).replace("RRDD","models.RRDD"))))
    print(repr(d_rcdd))
    print(repr(eval('sets.'+repr(d_rcdd).replace("RCDD","models.RCDD"))))
    print(repr(s_rdd))
    print(repr(eval('sets.'+repr(s_rdd).replace("RDD","models.RDD"))))
    print(repr(s_rrdd))
    print(repr(eval('sets.'+repr(s_rrdd).replace("RRDD","models.RRDD"))))
    print(repr(s_rcdd))
    print(repr(eval('sets.'+repr(s_rcdd).replace("RCDD","models.RCDD"))))
    print()

    # print the object
    print(d_rdd)
    print(d_rrdd)
    print(d_rcdd)
    print(s_rdd)
    print(s_rrdd)
    print(s_rcdd)
    print()

    # average over a sample
    v = lambda dist: dist.v
    print(str(s_rdd.average(v)))
    print(str(s_rrdd.average(v)))
    print(str(s_rcdd.average(v)))
    print()

    input("OK")
