#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module notation.
"""

from lpa.input import notation

# write numbers
print(notation.number(12345, 'id', 2))
print(notation.number(12345, 'id', 3))
print(notation.number(12345, 'id', 4))
print(notation.number(12345, 'id', 5))
print(notation.number(12345, 'id', 6))
print(notation.number(12345, 'id', 7), end="\n\n")
print(notation.number(12345, 'console', 2))
print(notation.number(12345, 'console', 3))
print(notation.number(12345, 'console', 4))
print(notation.number(12345, 'console', 5))
print(notation.number(12345, 'console', 6))
print(notation.number(12345, 'console', 7), end="\n\n")
print(notation.number(12345, 'tt', 2))
print(notation.number(12345, 'tt', 3))
print(notation.number(12345, 'tt', 4))
print(notation.number(12345, 'tt', 5))
print(notation.number(12345, 'tt', 6))
print(notation.number(12345, 'tt', 7), end="\n\n")

# write model parameters
r_none = {}
r_rrdd = {'v': 'R', 'd': 5e-3, 's': 200}
r_rcdd = {'v': 'D', 'd': 5e-3, 's': 134, 't': 14}
print("NONE"+notation.parameters(r_none, 'id')+"...")
print("RRDD"+notation.parameters(r_rrdd, 'id')+"...")
print("RCDD"+notation.parameters(r_rcdd, 'id')+"...", end="\n\n")
print("NONE"+notation.parameters(r_none, 'console')+"...")
print("RRDD"+notation.parameters(r_rrdd, 'console')+"...")
print("RCDD"+notation.parameters(r_rcdd, 'console')+"...", end="\n\n")
print("NONE"+notation.parameters(r_none, 'tt')+"...")
print("RRDD"+notation.parameters(r_rrdd, 'tt')+"...")
print("RCDD"+notation.parameters(r_rcdd, 'tt')+"...", end="\n\n")
