#!/usr/bin/env python3
"""Minimal test for site energy calculation without warnings."""

import sys
import os
sys.path.insert(0, '.')

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.data.parameters import create_parameter_config

# Quick test
pdb_path = "Alprotein/pdb/CP24_extended_most_occ_pH8.pdb"
config = {'ChlorophyllA': 'CLA'}

print("Loading structure with enhanced parameters...")
structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "test")
pigment_system = PigmentSystem(structure)

print(f"Pigments found: {len(pigment_system.pigments)}")
print(f"Enhanced atoms available: {pigment_system.enhanced_atoms_available}")

if len(pigment_system.pigments) > 0:
    first_pigment = list(pigment_system.pigments.values())[0]
    site_atoms = first_pigment.get_site_energy_atoms()
    print(f"Site atoms: {len(site_atoms)}")
    
    calculator = SiteEnergyCalculator()
    energy = calculator.calculate_site_energy(first_pigment, pigment_system)
    print(f"Site energy: {energy:.2f} cm⁻¹")
    print(f"Used direct access: {calculator.use_direct_access}")

print("✅ Test completed!")
