#!/usr/bin/env python3
"""Test parsing of 'None' string in charge field."""

# Try to parse a line with 'None' in the charge field
line = 'HETATM 9429   MG CLA A0501      -9.994  -4.667  -9.824                      Mg    None     01'

tokens = line.split()
print(f"Tokens: {tokens}")
print(f"Number of tokens: {len(tokens)}")

# Token 9 should be the charge
print(f"\nToken 9 (charge): '{tokens[9]}'")

# Try to convert to float
try:
    charge = float(tokens[9])
    print(f"Parsed charge: {charge}")
except ValueError as e:
    print(f"ERROR parsing charge: {e}")

# This is what my parser expects
print(f"\nExpected format has 11 tokens, got {len(tokens)}")
