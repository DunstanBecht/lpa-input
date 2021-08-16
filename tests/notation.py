#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module notation.
"""

from lpa.input import notation

# write numbers
print(notation.number(12345, 'stm', 2))
print(notation.number(12345, 'stm', 3))
print(notation.number(12345, 'stm', 4))
print(notation.number(12345, 'stm', 5))
print(notation.number(12345, 'stm', 6))
print(notation.number(12345, 'stm', 7), end="\n\n")
print(notation.number(12345, 'csl', 2))
print(notation.number(12345, 'csl', 3))
print(notation.number(12345, 'csl', 4))
print(notation.number(12345, 'csl', 5))
print(notation.number(12345, 'csl', 6))
print(notation.number(12345, 'csl', 7), end="\n\n")
print(notation.number(12345, 'ttl', 2))
print(notation.number(12345, 'ttl', 3))
print(notation.number(12345, 'ttl', 4))
print(notation.number(12345, 'ttl', 5))
print(notation.number(12345, 'ttl', 6))
print(notation.number(12345, 'ttl', 7), end="\n\n")

# write units
print(notation.unit("nm^{-3}", 'stm'))
print(notation.unit("nm^{-3}", 'csl'))
print(notation.unit("nm^{-3}", 'ttl'), end="\n\n")

# write quantities
print(notation.quantity(2, "nm", 'stm'))
print(notation.quantity(2, "nm", 'csl'))
print(notation.quantity(2, "nm", 'ttl'), end="\n\n")

# write equalities
print(notation.equality(r"\alpha", "10", 'stm'))
print(notation.equality(r"\alpha", "10", 'csl'))
print(notation.equality(r"\alpha", "10", 'ttl'), end="\n\n")

# write model parameters
r_none = {}
r_rrdd = {'v': 'R', 'd': 5e-3, 's': 200}
r_rcdd = {'v': 'D', 'd': 5e-3, 's': 134, 't': 14}
print("NONE"+notation.parameters(r_none, 'stm')+"...")
print("RRDD"+notation.parameters(r_rrdd, 'stm')+"...")
print("RCDD"+notation.parameters(r_rcdd, 'stm')+"...", end="\n\n")
print("NONE"+notation.parameters(r_none, 'csl')+"...")
print("RRDD"+notation.parameters(r_rrdd, 'csl')+"...")
print("RCDD"+notation.parameters(r_rcdd, 'csl')+"...", end="\n\n")
print("NONE"+notation.parameters(r_none, 'ttl')+"...")
print("RRDD"+notation.parameters(r_rrdd, 'ttl')+"...")
print("RCDD"+notation.parameters(r_rcdd, 'ttl')+"...", end="\n\n")
