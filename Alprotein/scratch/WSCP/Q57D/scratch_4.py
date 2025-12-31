#!/usr/bin/env python3
"""
Simple script to load a PDB file and calculate site energies using the optimized Alprotein architecture.

Usage:
    python calculate_site_energies.py [pdb_file]

If no PDB file is provided, it will use the default test file.
"""

import sys
import os
from pathlib import Path

# Add current directory to path if needed
if '.' not in sys.path:
    sys.path.insert(0, '.')

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.data.parameters import create_parameter_config

def calculate_site_energies(pdb_path, use_enhanced=True):
    """
    Load a PDB file and calculate site energies for all pigments.
    
    Args:
        pdb_path (str): Path to the PDB file
        use_enhanced (bool): Whether to use the enhanced architecture (recommended)
        
    Returns:
        dict: Site energies for each pigment
    """
    print(f"Loading PDB file: {pdb_path}")
    print(f"Using enhanced architecture: {use_enhanced}")
    print("=" * 50)
    
    try:
        if use_enhanced:
            # === NEW ENHANCED ARCHITECTURE (RECOMMENDED) ===
            print("🚀 Using optimized enhanced architecture...")
            
            # Step 1: Create parameter configuration
            # First, we need to peek at the structure to see what pigments are present
            temp_structure = ProteinStructure.from_file_extended(pdb_path, "temp")
            temp_system = PigmentSystem(temp_structure)
            pigment_resnames = [p.residue.resname for p in temp_system.pigments.values()]
            config = create_parameter_config(pigment_resnames)
            print(f"📋 Parameter configuration: {config}")
            
            # Step 2: Load structure with enhanced parameters
            structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "enhanced")
            pigment_system = PigmentSystem(structure)
            
        else:
            # === LEGACY ARCHITECTURE ===
            print("🐢 Using legacy architecture...")
            structure = ProteinStructure.from_file_extended(pdb_path, "legacy")
            pigment_system = PigmentSystem(structure)
        
        # Step 3: Display system information
        print(f"🔬 Structure loaded: {structure.name}")
        print(f"🧬 Pigments found: {len(pigment_system.pigments)}")
        print(f"⚡ Enhanced atoms available: {pigment_system.enhanced_atoms_available}")
        
        if len(pigment_system.pigments) == 0:
            print("❌ No pigments found in the structure!")
            return {}
        
        # Step 4: Create site energy calculator
        calculator = SiteEnergyCalculator(dielectric_constant=1.0)
        
        # Step 5: Calculate site energies for all pigments
        print("\n🧮 Calculating site energies...")
        print("-" * 50)
        
        site_energies = {}
        total_site_atoms = 0
        
        for i, (name, pigment) in enumerate(pigment_system.pigments.items()):
            print(f"\nPigment {i+1}: {name}")
            print(f"  Type: {type(pigment).__name__}")
            print(f"  Residue: {pigment.get_resname()}")
            
            # Get site energy atoms info
            site_atoms = pigment.get_site_energy_atoms()
            print(f"  Site atoms: {len(site_atoms)}")
            total_site_atoms += len(site_atoms)
            
            if len(site_atoms) > 0:
                # Calculate site energy
                energy = calculator.calculate_site_energy(pigment, pigment_system)
                site_energies[name] = energy
                
                # Get vacuum energy for comparison
                vacuum_energy = pigment.get_vacuum_energy()
                shift = energy - vacuum_energy
                
                print(f"  Vacuum energy: {vacuum_energy:.2f} cm⁻¹")
                print(f"  Site energy: {energy:.2f} cm⁻¹")
                print(f"  Energy shift: {shift:+.2f} cm⁻¹")
                print(f"  Used direct access: {calculator.use_direct_access}")
            else:
                print("  ⚠️  No site energy atoms available - skipping calculation")
                site_energies[name] = None
        
        # Step 6: Summary
        print("\n" + "=" * 50)
        print("📊 CALCULATION SUMMARY")
        print("=" * 50)
        
        calculated_pigments = [name for name, energy in site_energies.items() if energy is not None]
        
        if calculated_pigments:
            energies = [site_energies[name] for name in calculated_pigments]
            print(f"✅ Successfully calculated: {len(calculated_pigments)}/{len(pigment_system.pigments)} pigments")
            print(f"🔢 Total site atoms used: {total_site_atoms}")
            print(f"⚡ Architecture used: {'Enhanced (O(1) access)' if calculator.use_direct_access else 'Legacy (O(n) lookup)'}")
            print(f"🎯 Energy range: {min(energies):.2f} to {max(energies):.2f} cm⁻¹")
            print(f"📈 Energy spread: {max(energies) - min(energies):.2f} cm⁻¹")
            
            # Show top contributors to red/blue shifts
            vacuum_ref = 14900.0  # Typical chlorophyll A vacuum energy
            shifts = [(name, site_energies[name] - vacuum_ref) for name in calculated_pigments]
            shifts.sort(key=lambda x: x[1])
            
            print(f"\n🔴 Most red-shifted: {shifts[0][0]} ({shifts[0][1]:+.2f} cm⁻¹)")
            print(f"🔵 Most blue-shifted: {shifts[-1][0]} ({shifts[-1][1]:+.2f} cm⁻¹)")
        else:
            print("❌ No site energies could be calculated")
        
        return site_energies
        
    except Exception as e:
        print(f"❌ Error during calculation: {e}")
        import traceback
        traceback.print_exc()
        return {}



print("🧬 Alprotein Site Energy Calculator")
print("=" * 50)

pdb_path = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57D/most_occ/most_occ_pH1.pdb'
# pdb_path = "/Users/mohamed/Documents/Research/Projects/2025/Alprotein-Alpha/Alprotein/scratch/most_occ_pH1_one_pig.pdb"
temp_structure = ProteinStructure.from_file_extended(pdb_path, "temp")
temp_system = PigmentSystem(temp_structure)
pigment_resnames = [p.residue.resname for p in temp_system.pigments.values()]
config = create_parameter_config(pigment_resnames)
print(f"📋 Parameter configuration: {config}")

# Step 2: Load structure with enhanced parameters
structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "enhanced")
pigment_system = PigmentSystem(structure)

# Step 3: Display system information
print(f"🔬 Structure loaded: {structure.name}")
print(f"🧬 Pigments found: {len(pigment_system.pigments)}")
print(f"⚡ Enhanced atoms available: {pigment_system.enhanced_atoms_available}")


# Step 4: Create site energy calculator
calculator = SiteEnergyCalculator(dielectric_constant=2.5)

# Step 5: Calculate site energies for all pigments
print("\n🧮 Calculating site energies...")
print("-" * 50)

site_energies = {}
total_site_atoms = 0

for i, (name, pigment) in enumerate(pigment_system.pigments.items()):
    print(f"\nPigment {i+1}: {name}")
    print(f"  Type: {type(pigment).__name__}")
    print(f"  Residue: {pigment.get_resname()}")
    
    # Get site energy atoms info
    site_atoms = pigment.get_site_energy_atoms()
    print(f"  Site atoms: {len(site_atoms)}")
    total_site_atoms += len(site_atoms)
    
    if len(site_atoms) > 0:
        # Calculate site energy
        energy = calculator.calculate_site_energy(pigment, pigment_system)
        site_energies[name] = energy
        
        # Get vacuum energy for comparison
        vacuum_energy = pigment.get_vacuum_energy()
        shift = energy - vacuum_energy
        
        print(f"  Vacuum energy: {vacuum_energy:.2f} cm⁻¹")
        print(f"  Site energy: {energy:.2f} cm⁻¹")
        print(f"  Energy shift: {shift:+.2f} cm⁻¹")
        print(f"  Used direct access: {calculator.use_direct_access}")
    else:
        print("  ⚠️  No site energy atoms available - skipping calculation")
        site_energies[name] = None

# Initialize site energy calculator
calculator = SiteEnergyCalculator(dielectric_constant=2.5)

print(f"\n📊 Calculating detailed residue contributions...")

# Method 1: Basic residue contributions
site_energies, contributions = calculator.calculate_detailed_site_energy_contributions(pigment_system)

# correct results:
# -----------------
# {'wscp_A_CLA_1001': 0.21145102673236932}
# {'A_ASN_2': 0.21145102673236932}

# # =================================================================
# # METHOD 1: Access charges from pigment atoms (RECOMMENDED)
# # =================================================================
# print("\n" + "="*50)
# print("METHOD 1: Access charges from PIGMENT atoms")
# print("="*50)

# if len(pigment_system.pigments) > 0:
#     # Get first pigment
#     first_pigment = list(pigment_system.pigments.values())[0]
#     print(f"\n🧬 Pigment: {first_pigment.name}")
#     print(f"   Type: {type(first_pigment).__name__}")
#     print(f"   Residue: {first_pigment.get_resname()}")
    
#     print(f"\n📊 Charge access methods for pigment atoms:")
    
#     # Method 1a: Direct attribute access (O(1) - FASTEST)
#     print(f"\n1a. Direct attribute access (Enhanced atoms):")
#     enhanced_atoms_count = 0
#     for i, atom in enumerate(first_pigment.atoms):
#         if hasattr(atom, 'calculation_ready') and atom.calculation_ready:
#             enhanced_atoms_count += 1
#             if i < 5:  # Show first 5 examples
#                 print(f"    atom.name: {atom.name}")
#                 print(f"    atom.q00: {atom.q00:.6f}")
#                 print(f"    atom.q11: {atom.q11:.6f}")
#                 if hasattr(atom, 'q01') and atom.q01 is not None:
#                     print(f"    atom.q01: {atom.q01:.6f}")
#                 print(f"    atom.calculation_ready: {atom.calculation_ready}")
#                 print()
#             elif i == 5:
#                 print(f"    ... (showing first 5 of {enhanced_atoms_count} enhanced atoms)")
#                 break
    
#     # Method 1b: Using pigment's optimized method
#     print(f"\n1b. Using pigment's optimized method:")
#     site_atoms = first_pigment.get_site_energy_atoms_direct()
#     print(f"    Total site atoms: {len(site_atoms)}")
#     for i, (atom_name, q00, q11) in enumerate(site_atoms[:3]):
#         print(f"    {atom_name}: q00={q00:.6f}, q11={q11:.6f}")
#     if len(site_atoms) > 3:
#         print(f"    ... (showing first 3 of {len(site_atoms)} atoms)")
    
#     # Method 1c: Using standard pigment method (with automatic optimization)
#     print(f"\n1c. Using standard pigment method (auto-optimized):")
#     site_atoms_standard = first_pigment.get_site_energy_atoms()
#     print(f"    Total site atoms: {len(site_atoms_standard)}")
#     for i, (atom_name, q00, q11) in enumerate(site_atoms_standard[:3]):
#         print(f"    {atom_name}: q00={q00:.6f}, q11={q11:.6f}")
    
#     # Method 1d: Filter atoms by charge type
#     print(f"\n1d. Filter atoms by charge type:")
#     q00_atoms = first_pigment.get_atoms_with_charges('q_00')
#     q11_atoms = first_pigment.get_atoms_with_charges('q_11')
#     print(f"    Atoms with q_00: {len(q00_atoms)}")
#     print(f"    Atoms with q_11: {len(q11_atoms)}")
    
#     # Example: Get specific atom by name
#     print(f"\n1e. Get specific atom by name:")
#     try:
#         mg_atom = first_pigment.residue['MG']
#         if hasattr(mg_atom, 'q00'):
#             print(f"    MG atom q00: {mg_atom.q00:.6f}")
#             print(f"    MG atom q11: {mg_atom.q11:.6f}")
#         else:
#             print(f"    MG atom has no q00 charge (not in parameter set)")
#     except KeyError:
#         print(f"    No MG atom found in this pigment")

# # =================================================================
# # METHOD 2: Access charges from protein background atoms
# # =================================================================
# print("\n" + "="*50)
# print("METHOD 2: Access charges from PROTEIN background atoms")
# print("="*50)

# # Get protein background atoms (non-pigment atoms with charges)
# background_atoms = pigment_system.get_protein_background_atoms()
# print(f"\n🧬 Background atoms found: {len(background_atoms)}")

# if len(background_atoms) > 0:
#     print(f"\n📊 Charge access for protein atoms:")
    
#     # Method 2a: Direct charge attribute access
#     print(f"\n2a. Direct charge attribute access:")
#     charged_count = 0
#     for i, atom in enumerate(background_atoms):
#         if hasattr(atom, 'charge') and atom.charge != 0.0:
#             charged_count += 1
#             if i < 5:  # Show first 5 examples
#                 residue = atom.get_parent()
#                 print(f"    {residue.resname}_{residue.id[1]}_{atom.name}: charge={atom.charge:.6f}")
#             elif i == 5:
#                 print(f"    ... (showing first 5 of {charged_count} charged atoms)")
#                 break
    
#     # Method 2b: Check for extended charges (q00, q11 from extended PDB)
#     print(f"\n2b. Check for extended charges in protein atoms:")
#     extended_count = 0
#     for atom in background_atoms:
#         if hasattr(atom, 'q_00') or hasattr(atom, 'q00'):
#             extended_count += 1
#             if extended_count <= 3:
#                 residue = atom.get_parent()
#                 q00_val = getattr(atom, 'q_00', getattr(atom, 'q00', 'N/A'))
#                 print(f"    {residue.resname}_{residue.id[1]}_{atom.name}: q00={q00_val}")
    
#     if extended_count == 0:
#         print(f"    No extended charges (q00/q11) found in protein atoms")
#         print(f"    (This is normal - protein atoms typically only have 'charge' attribute)")
#     else:
#         print(f"    Found {extended_count} protein atoms with extended charges")
