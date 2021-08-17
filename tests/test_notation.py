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
r1 = {}
r2 = {'d': 0}
r3 = {'s': 0, 'f': 0, 't': 0, 'l': 0, 'x': 0}
print("NONE"+notation.parameters(r1, 'stm')+"...")
print("RRDD"+notation.parameters(r2, 'stm')+"...")
print("RCDD"+notation.parameters(r3, 'stm')+"...", end="\n\n")
print("NONE"+notation.parameters(r1, 'csl')+"...")
print("RRDD"+notation.parameters(r2, 'csl')+"...")
print("RCDD"+notation.parameters(r3, 'csl')+"...", end="\n\n")
print("NONE"+notation.parameters(r1, 'ttl')+"...")
print("RRDD"+notation.parameters(r2, 'ttl')+"...")
print("RCDD"+notation.parameters(r3, 'ttl')+"...", end="\n\n")

input("OK")
