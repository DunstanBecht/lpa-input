#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module notation.
"""

from lpa.input import notation

# write numbers
print(notation.number(12345, 'file', 2))
print(notation.number(12345, 'file', 3))
print(notation.number(12345, 'file', 4))
print(notation.number(12345, 'file', 5))
print(notation.number(12345, 'file', 6))
print(notation.number(12345, 'file', 7), end="\n\n")
print(notation.number(12345, 'console', 2))
print(notation.number(12345, 'console', 3))
print(notation.number(12345, 'console', 4))
print(notation.number(12345, 'console', 5))
print(notation.number(12345, 'console', 6))
print(notation.number(12345, 'console', 7), end="\n\n")
print(notation.number(12345, 'title', 2))
print(notation.number(12345, 'title', 3))
print(notation.number(12345, 'title', 4))
print(notation.number(12345, 'title', 5))
print(notation.number(12345, 'title', 6))
print(notation.number(12345, 'title', 7), end="\n\n")

# write model parameters
r_none = {}
r_rrdd = {'v': 'R', 'd': 5e-3, 's': 200}
r_rcdd = {'v': 'D', 'd': 5e-3, 's': 134, 't': 14}
print("NONE"+notation.parameters(r_none)+"...")
print("RRDD"+notation.parameters(r_rrdd)+"...")
print("RCDD"+notation.parameters(r_rcdd)+"...", end="\n\n")
print("NONE"+notation.parameters(r_none, 'console')+"...")
print("RRDD"+notation.parameters(r_rrdd, 'console')+"...")
print("RCDD"+notation.parameters(r_rcdd, 'console')+"...", end="\n\n")
print("NONE"+notation.parameters(r_none, 'title')+"...")
print("RRDD"+notation.parameters(r_rrdd, 'title')+"...")
print("RCDD"+notation.parameters(r_rcdd, 'title')+"...", end="\n\n")
