#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module maps.
"""

from lpa.input import maps
from test_sets import *

# export distribution maps
maps.export(d_rdd, exfmt='svg')
maps.export(d_rrdd, exfmt='svg')
maps.export(d_rcdd, exdir='maps/', ttlsp='suptitle', ttlsb='subtitle')

input("OK")
