#!/usr/bin/env python3
"""
Test: Validate the corrected CDC equation implementation.

This script verifies the corrected CDC equation:
(TRESP_CC / dielectric_eff) * np.sum(Q2_pigA_tile * Q2_background_tile / R2_ab_norm)
where TRESP_CC = 1.1615E5 # Units: Angstrom*cm^-1)/e^2
"""

import sys
import os
import numpy as np
import time
sys.path.insert(0, '.')

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.core.constants import TRESP_CC

pdb_path = '/Users/mohamed/Documents/Research/Projects/2025/Alprotein-Alpha/Alprotein/scratch/most_occ_pH1_one_pig.pdb'
def test_cdc_equation():
    """Test the corrected CDC equation implementation."""
    
    print("🧮 Testing Corrected CDC Equation")
    print("=" * 50)
    print(f"TRESP_CC = {TRESP_CC:.2e} (Angstrom*cm⁻¹)/e²")
    
    # Load structure with enhanced parameters
    # pdb_path = "Alprotein/pdb/CP24_extended_most_occ_pH8.pdb"
    # pdb_path = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57D/most_occ/most_occ_pH1.pdb'

    config = {'ChlorophyllA': 'CLA'}
    
    structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "test")
    pigment_system = PigmentSystem(structure)
    
    print(f"✓ Structure loaded: {len(pigment_system.pigments)} pigments")
    
    # Get background atoms
    background_atoms = pigment_system.get_protein_background_atoms()
    print(f"✓ Background atoms: {len(background_atoms)}")
    
    if len(background_atoms) == 0:
        print("❌ No background atoms found! Cannot test CDC equation.")
        return
    
    # Initialize site energy calculator with dielectric = 1.0 for testing
    calculator = SiteEnergyCalculator(dielectric_constant=1.0)
    
    # Test on first pigment
    test_pigment = list(pigment_system.pigments.values())[0]
    pigment_name = test_pigment.get_id()
    print(f"\n🎯 Testing CDC calculation for {pigment_name}")
    
    # Get site energy atoms
    site_atoms_data = test_pigment.get_site_energy_atoms()
    print(f"Site energy atoms: {len(site_atoms_data)}")
    
    if not site_atoms_data:
        print("❌ No site energy atoms found!")
        return
    
    # Show first few site atoms
    print(f"\n📊 Site Energy Atoms (first 5):")
    for i, (atom_name, q00, q11) in enumerate(site_atoms_data[:5]):
        delta_q = q11 - q00
        print(f"  {atom_name:<8} q00={q00:+7.3f} q11={q11:+7.3f} Δq={delta_q:+7.3f}")
    
    # Manual CDC calculation for verification
    print(f"\n🔍 Manual CDC Verification:")
    
    # Extract data
    atom_names = [data[0] for data in site_atoms_data]
    q_00_charges = np.array([data[1] for data in site_atoms_data])
    q_11_charges = np.array([data[2] for data in site_atoms_data])
    delta_q = q_11_charges - q_00_charges
    
    # Get positions
    target_pos = np.array([test_pigment.get_atom_coord(name) for name in atom_names])
    background_pos = np.array([atom.coord for atom in background_atoms[:100]])  # Use first 100 for speed
    background_charges = np.array([getattr(atom, 'charge', 0.0) for atom in background_atoms[:100]])
    
    print(f"Pigment atoms: {len(delta_q)}")
    print(f"Background atoms: {len(background_charges)}")
    print(f"Non-zero background charges: {np.sum(np.abs(background_charges) > 1e-6)}")
    
    # Method 1: Old nested loop method (for comparison)
    print(f"\n📈 Method 1: Old nested loop")
    old_result = 0.0
    for i, dq_i in enumerate(delta_q):
        for j, q_j in enumerate(background_charges):
            dist = np.linalg.norm(target_pos[i] - background_pos[j])
            if dist > 0:
                old_result += (dq_i * q_j) / dist
    old_result = (TRESP_CC / 1.0) * old_result
    print(f"Old method result: {old_result:+10.2f} cm⁻¹")
    
    # Method 2: New vectorized CDC equation
    print(f"\n📈 Method 2: New vectorized CDC")
    
    # Create tiled arrays
    Q2_pigA_tile = np.tile(delta_q.reshape(-1, 1), (1, len(background_charges)))
    Q2_background_tile = np.tile(background_charges.reshape(1, -1), (len(delta_q), 1))
    
    print(f"Q2_pigA_tile shape: {Q2_pigA_tile.shape}")
    print(f"Q2_background_tile shape: {Q2_background_tile.shape}")
    
    # Calculate distance matrix
    R2_ab_norm = np.linalg.norm(
        target_pos[:, np.newaxis, :] - background_pos[np.newaxis, :, :], 
        axis=2
    )
    print(f"R2_ab_norm shape: {R2_ab_norm.shape}")
    print(f"Min distance: {np.min(R2_ab_norm):.3f} Å")
    print(f"Max distance: {np.max(R2_ab_norm):.3f} Å")
    
    # Avoid division by zero
    R2_ab_norm = np.where(R2_ab_norm > 0, R2_ab_norm, np.inf)
    
    # Apply CDC equation
    new_result = (TRESP_CC / 1.0) * np.sum(
        Q2_pigA_tile * Q2_background_tile / R2_ab_norm
    )
    print(f"New method result: {new_result:+10.2f} cm⁻¹")
    
    # Method 3: Use the actual calculator method
    print(f"\n📈 Method 3: Calculator method")
    calc_result = calculator._pigment_background_interaction(
        test_pigment, background_atoms[:100], background_charges.tolist()
    )
    print(f"Calculator result: {calc_result:+10.2f} cm⁻¹")
    
    # Compare results
    print(f"\n📊 Comparison:")
    print(f"Old vs New: {abs(old_result - new_result):.6f} cm⁻¹ difference")
    print(f"New vs Calc: {abs(new_result - calc_result):.6f} cm⁻¹ difference")
    
    if abs(old_result - new_result) < 1e-6:
        print("✅ Vectorized implementation matches old method!")
    else:
        print("❌ Vectorized implementation differs from old method!")
    
    if abs(new_result - calc_result) < 1e-6:
        print("✅ Calculator method matches vectorized implementation!")
    else:
        print("❌ Calculator method differs from vectorized implementation!")
    
    # Test full site energy calculation
    print(f"\n🎯 Full Site Energy Calculation:")
    
    vacuum_energy = test_pigment.get_vacuum_energy()
    site_energy = calculator.calculate_site_energy(test_pigment, pigment_system)
    shift = site_energy - vacuum_energy
    
    print(f"Vacuum energy: {vacuum_energy:8.1f} cm⁻¹")
    print(f"Site energy:   {site_energy:8.1f} cm⁻¹")
    print(f"Total shift:   {shift:+8.1f} cm⁻¹")
    
    # Breakdown of contributions
    protein_contribution = calculator._pigment_background_interaction(
        test_pigment, background_atoms, [getattr(atom, 'charge', 0.0) for atom in background_atoms]
    )
    
    pigment_contribution = 0.0
    for other_pigment in pigment_system.pigments.values():
        if other_pigment != test_pigment:
            pigment_contribution += calculator._pigment_pigment_interaction(test_pigment, other_pigment)
    
    print(f"\n📋 Contribution Breakdown:")
    print(f"Protein background: {protein_contribution:+8.1f} cm⁻¹")
    print(f"Other pigments:     {pigment_contribution:+8.1f} cm⁻¹")
    print(f"Total shift:        {protein_contribution + pigment_contribution:+8.1f} cm⁻¹")
    
    print(f"\n✅ CDC equation test completed!")


def compare_optimization_methods():
    """Compare optimized vs legacy calculation methods."""
    
    print(f"\n🚀 Comparing Optimization Methods")
    print("=" * 40)
    
    # Load structure
    # pdb_path = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57D/most_occ/most_occ_pH1.pdb'
    config = {'ChlorophyllA': 'CLA'}
    
    structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "test")
    pigment_system = PigmentSystem(structure)
    
    calculator = SiteEnergyCalculator(dielectric_constant=1.0)
    
    # Test first pigment
    test_pigment = list(pigment_system.pigments.values())[0]
    pigment_name = test_pigment.get_id()
    
    print(f"Testing {pigment_name} with both methods...")
    
    # Method 1: Legacy calculation
    start_time = time.time()
    legacy_energy = calculator._calculate_site_energy_legacy(test_pigment, pigment_system)
    legacy_time = time.time() - start_time
    
    # Method 2: Optimized calculation
    start_time = time.time()
    optimized_energy = calculator._calculate_site_energy_optimized(test_pigment, pigment_system)
    optimized_time = time.time() - start_time
    
    print(f"\n📊 Results:")
    print(f"Legacy method:    {legacy_energy:8.1f} cm⁻¹ ({legacy_time:.4f}s)")
    print(f"Optimized method: {optimized_energy:8.1f} cm⁻¹ ({optimized_time:.4f}s)")
    print(f"Difference:       {abs(legacy_energy - optimized_energy):8.1f} cm⁻¹")
    if optimized_time > 0:
        print(f"Speedup:          {legacy_time/optimized_time:.1f}x faster")
    
    # Test all pigments
    print(f"\n🎯 Testing All Pigments:")
    print(f"{'Pigment':<15} {'Legacy':>10} {'Optimized':>10} {'Diff':>8}")
    print("-" * 50)
    
    for pigment in pigment_system.pigments.values():
        legacy = calculator._calculate_site_energy_legacy(pigment, pigment_system)
        optimized = calculator._calculate_site_energy_optimized(pigment, pigment_system)
        diff = abs(legacy - optimized)
        
        print(f"{pigment.get_id():<15} {legacy:8.1f} {optimized:8.1f} {diff:6.1f}")


def test_residue_contributions():
    """Test detailed residue contribution analysis."""
    
    print(f"\n🔬 Testing Residue Contribution Analysis")
    print("=" * 50)
    
    # Load structure
    # pdb_path = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57D/most_occ/most_occ_pH1.pdb'
    config = {'ChlorophyllA': 'CLA'}
    
    structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "test")
    pigment_system = PigmentSystem(structure)
    calculator = SiteEnergyCalculator(dielectric_constant=1.0)
    
    # Test detailed contributions
    site_energies, contributions = calculator.calculate_detailed_site_energy_contributions(pigment_system)
    
    # Show results for first pigment
    test_pigment_name = list(pigment_system.get_pigment_names())[0]
    test_contributions = contributions[test_pigment_name]
    
    print(f"\n🎯 Residue Contributions for {test_pigment_name}:")
    
    # Sort by absolute contribution
    residue_contribs = []
    for contrib_name, contrib_value in test_contributions.items():
        if contrib_name != 'vacuum':
            residue_contribs.append((contrib_name, contrib_value))
    
    residue_contribs.sort(key=lambda x: abs(x[1]), reverse=True)
    
    print(f"{'Residue/Pigment':<25} {'Contribution':>12}")
    print("-" * 40)
    
    for contrib_name, contrib_value in residue_contribs[:15]:  # Top 15
        print(f"{contrib_name:<25} {contrib_value:+10.1f} cm⁻¹")
    
    # Summary statistics
    total_contributions = [abs(c[1]) for c in residue_contribs]
    print(f"\n📊 Summary:")
    print(f"Total contributing entities: {len(residue_contribs)}")
    print(f"Strongest contribution: {max(total_contributions):.1f} cm⁻¹")
    print(f"Average contribution: {np.mean(total_contributions):.1f} cm⁻¹")
    print(f"Contributions > 10 cm⁻¹: {sum(1 for c in total_contributions if c > 10)}")
    print(f"Contributions > 50 cm⁻¹: {sum(1 for c in total_contributions if c > 50)}")


def validate_units_and_constants():
    """Validate that units and constants are correct."""
    
    print(f"\n🔧 Validating Units and Constants")
    print("=" * 40)
    
    print(f"TRESP_CC = {TRESP_CC:.2e}")
    print(f"Expected: 1.1615E5 (Angstrom*cm⁻¹)/e²")
    
    if abs(TRESP_CC - 1.1615E5) < 1e-6:
        print("✅ TRESP_CC constant is correct!")
    else:
        print("❌ TRESP_CC constant is incorrect!")
    
    # Test unit consistency with a simple calculation
    print(f"\n🧮 Unit consistency test:")
    
    # Simple test: two unit charges 1 Angstrom apart
    charge1 = 1.0  # e
    charge2 = 1.0  # e
    distance = 1.0  # Angstrom
    dielectric = 1.0
    
    energy = (TRESP_CC / dielectric) * (charge1 * charge2) / distance
    print(f"Energy for two unit charges 1Å apart: {energy:.1f} cm⁻¹")
    print(f"Expected order of magnitude: ~10⁵ cm⁻¹")
    
    if 50000 < energy < 200000:
        print("✅ Units and magnitude are reasonable!")
    else:
        print("❌ Units or magnitude may be incorrect!")


if __name__ == "__main__":
    try:
        test_cdc_equation()
        compare_optimization_methods()
        test_residue_contributions()
        validate_units_and_constants()
        
        print(f"\n🎉 All tests completed!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()