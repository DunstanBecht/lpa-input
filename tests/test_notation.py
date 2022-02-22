#!/usr/bin/env python
# coding: utf-8

"""
Script to test the module notation.
"""

from lpa.input import notation

x = -12.5 # example of number to display

"""
The following lines display x in a representation suitable for use in
a file name.
"""
print("Numbers (stm)")
print(notation.number(x, 'stm', 2))
print(notation.number(x, 'stm', 3))
print(notation.number(x, 'stm', 4))
print(notation.number(x, 'stm', 5))
print(notation.number(x, 'stm', 6))
print(notation.number(x, 'stm', 7))
print()

"""
The following lines display x in a representation suitable for use in
a terminal or a text file.
"""
print("Numbers (csl)")
print(notation.number(x, 'csl', 2))
print(notation.number(x, 'csl', 3))
print(notation.number(x, 'csl', 4))
print(notation.number(x, 'csl', 5))
print(notation.number(x, 'csl', 6))
print(notation.number(x, 'csl', 7))
print()

"""
The following lines display x in a representation suitable for use in
a LaTeX title.
"""
print("Numbers (ttl)")
print(notation.number(x, 'ttl', 2))
print(notation.number(x, 'ttl', 3))
print(notation.number(x, 'ttl', 4))
print(notation.number(x, 'ttl', 5))
print(notation.number(x, 'ttl', 6))
print(notation.number(x, 'ttl', 7))
print()

"""
The following lines display the unit 'nm^-3' in the three available
representations.
"""
print("Units (stm/csl/ttl)")
print(notation.unit("nm^{-3}", 'stm'))
print(notation.unit("nm^{-3}", 'csl'))
print(notation.unit("nm^{-3}", 'ttl'))
print()

"""
The following lines display the physical quantity '2 nm' in the three
available representations.
"""
print("Quantities (stm/csl/ttl)")
print(notation.quantity(2, "nm", 'stm'))
print(notation.quantity(2, "nm", 'csl'))
print(notation.quantity(2, "nm", 'ttl'))
print()

"""
The following lines display the equality 'alpha = 10' in the three
available representations.
"""
print("Equalities (stm/csl/ttl)")
print(notation.equality(r"\alpha", "10", 'stm'))
print(notation.equality(r"\alpha", "10", 'csl'))
print(notation.equality(r"\alpha", "10", 'ttl'))
print()

"""
The following lines define some examples of model parameter sets.
"""
r1 = {}
r2 = {'d': 0}
r3 = {'s': 0, 'f': 0, 't': 0, 'l': 0, 'x': 0}

"""
The following lines display the model parameter sets in a
representation suitable for use in a file name.
"""
print("Model parameters (stm)")
print("NONE"+notation.parameters(r1, 'stm')+"...")
print("RRDD"+notation.parameters(r2, 'stm')+"...")
print("RCDD"+notation.parameters(r3, 'stm')+"...")
print()

"""
The following lines display the model parameter sets in a
representation suitable for use in a terminal or a text file.
"""
print("Model parameters (csl)")
print("NONE"+notation.parameters(r1, 'csl')+"...")
print("RRDD"+notation.parameters(r2, 'csl')+"...")
print("RCDD"+notation.parameters(r3, 'csl')+"...")
print()

"""
The following lines display the model parameter sets in a
representation suitable for use in a LaTeX title.
"""
print("Model parameters (ttl)")
print("NONE"+notation.parameters(r1, 'ttl')+"...")
print("RRDD"+notation.parameters(r2, 'ttl')+"...")
print("RCDD"+notation.parameters(r3, 'ttl')+"...")
print()

input("OK")
