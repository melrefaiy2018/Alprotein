#!/usr/bin/env python3
"""
Guide: How to access q00 charges directly from atoms in protein and pigment contexts.

This script demonstrates different ways to access atom charges after loading
a structure with the enhanced Alprotein architecture.
"""

import sys
import os
sys.path.insert(0, '.')

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.data.parameters import create_parameter_config

def demonstrate_charge_access():
    """Demonstrate different ways to access atom charges."""
    
    print("🧬 Direct Atom Charge Access Guide")
    print("=" * 50)
    
    # Load structure with enhanced parameters
    pdb_path = "/Users/mohamed/Documents/Research/Projects/2025/Alprotein-Alpha/Alprotein/scratch/WSCP/Q57D/most_occ_pH1.pdb"
    config = {'ChlorophyllA': 'CLA'}
    
    print("Loading structure with enhanced parameters...")
    structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "enhanced")
    pigment_system = PigmentSystem(structure)
    
    print(f"✓ Structure loaded with {len(pigment_system.pigments)} pigments")
    print(f"✓ Enhanced atoms available: {pigment_system.enhanced_atoms_available}")
    
    # =================================================================
    # METHOD 1: Access charges from pigment atoms (RECOMMENDED)
    # =================================================================
    print("\n" + "="*50)
    print("METHOD 1: Access charges from PIGMENT atoms")
    print("="*50)
    
    if len(pigment_system.pigments) > 0:
        # Get first pigment
        first_pigment = list(pigment_system.pigments.values())[0]
        print(f"\n🧬 Pigment: {first_pigment.name}")
        print(f"   Type: {type(first_pigment).__name__}")
        print(f"   Residue: {first_pigment.get_resname()}")
        
        print(f"\n📊 Charge access methods for pigment atoms:")
        
        # Method 1a: Direct attribute access (O(1) - FASTEST)
        print(f"\n1a. Direct attribute access (Enhanced atoms):")
        enhanced_atoms_count = 0
        for i, atom in enumerate(first_pigment.atoms):
            if hasattr(atom, 'calculation_ready') and atom.calculation_ready:
                enhanced_atoms_count += 1
                if i < 5:  # Show first 5 examples
                    print(f"    atom.name: {atom.name}")
                    print(f"    atom.q00: {atom.q00:.6f}")
                    print(f"    atom.q11: {atom.q11:.6f}")
                    if hasattr(atom, 'q01') and atom.q01 is not None:
                        print(f"    atom.q01: {atom.q01:.6f}")
                    print(f"    atom.calculation_ready: {atom.calculation_ready}")
                    print()
                elif i == 5:
                    print(f"    ... (showing first 5 of {enhanced_atoms_count} enhanced atoms)")
                    break
        
        # Method 1b: Using pigment's optimized method
        print(f"\n1b. Using pigment's optimized method:")
        site_atoms = first_pigment.get_site_energy_atoms_direct()
        print(f"    Total site atoms: {len(site_atoms)}")
        for i, (atom_name, q00, q11) in enumerate(site_atoms[:3]):
            print(f"    {atom_name}: q00={q00:.6f}, q11={q11:.6f}")
        if len(site_atoms) > 3:
            print(f"    ... (showing first 3 of {len(site_atoms)} atoms)")
        
        # Method 1c: Using standard pigment method (with automatic optimization)
        print(f"\n1c. Using standard pigment method (auto-optimized):")
        site_atoms_standard = first_pigment.get_site_energy_atoms()
        print(f"    Total site atoms: {len(site_atoms_standard)}")
        for i, (atom_name, q00, q11) in enumerate(site_atoms_standard[:3]):
            print(f"    {atom_name}: q00={q00:.6f}, q11={q11:.6f}")
        
        # Method 1d: Filter atoms by charge type
        print(f"\n1d. Filter atoms by charge type:")
        q00_atoms = first_pigment.get_atoms_with_charges('q_00')
        q11_atoms = first_pigment.get_atoms_with_charges('q_11')
        print(f"    Atoms with q_00: {len(q00_atoms)}")
        print(f"    Atoms with q_11: {len(q11_atoms)}")
        
        # Example: Get specific atom by name
        print(f"\n1e. Get specific atom by name:")
        try:
            mg_atom = first_pigment.residue['MG']
            if hasattr(mg_atom, 'q00'):
                print(f"    MG atom q00: {mg_atom.q00:.6f}")
                print(f"    MG atom q11: {mg_atom.q11:.6f}")
            else:
                print(f"    MG atom has no q00 charge (not in parameter set)")
        except KeyError:
            print(f"    No MG atom found in this pigment")
    
    # =================================================================
    # METHOD 2: Access charges from protein background atoms
    # =================================================================
    print("\n" + "="*50)
    print("METHOD 2: Access charges from PROTEIN background atoms")
    print("="*50)
    
    # Get protein background atoms (non-pigment atoms with charges)
    background_atoms = pigment_system.get_protein_background_atoms()
    print(f"\n🧬 Background atoms found: {len(background_atoms)}")
    
    if len(background_atoms) > 0:
        print(f"\n📊 Charge access for protein atoms:")
        
        # Method 2a: Direct charge attribute access
        print(f"\n2a. Direct charge attribute access:")
        charged_count = 0
        for i, atom in enumerate(background_atoms):
            if hasattr(atom, 'charge') and atom.charge != 0.0:
                charged_count += 1
                if i < 5:  # Show first 5 examples
                    residue = atom.get_parent()
                    print(f"    {residue.resname}_{residue.id[1]}_{atom.name}: charge={atom.charge:.6f}")
                elif i == 5:
                    print(f"    ... (showing first 5 of {charged_count} charged atoms)")
                    break
        
        # Method 2b: Check for extended charges (q00, q11 from extended PDB)
        print(f"\n2b. Check for extended charges in protein atoms:")
        extended_count = 0
        for atom in background_atoms:
            if hasattr(atom, 'q_00') or hasattr(atom, 'q00'):
                extended_count += 1
                if extended_count <= 3:
                    residue = atom.get_parent()
                    q00_val = getattr(atom, 'q_00', getattr(atom, 'q00', 'N/A'))
                    print(f"    {residue.resname}_{residue.id[1]}_{atom.name}: q00={q00_val}")
        
        if extended_count == 0:
            print(f"    No extended charges (q00/q11) found in protein atoms")
            print(f"    (This is normal - protein atoms typically only have 'charge' attribute)")
        else:
            print(f"    Found {extended_count} protein atoms with extended charges")
    
    # =================================================================
    # METHOD 3: Access charges from all atoms in structure
    # =================================================================
    print("\n" + "="*50)
    print("METHOD 3: Access charges from ALL atoms in structure")
    print("="*50)
    
    print(f"\n🧬 Scanning all atoms in structure:")
    
    # Method 3a: Iterate through all atoms
    total_atoms = 0
    atoms_with_charge = 0
    atoms_with_q00 = 0
    atoms_with_q11 = 0
    
    for model in structure.pdb:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    total_atoms += 1
                    
                    # Check for standard charge
                    if hasattr(atom, 'charge') and atom.charge != 0.0:
                        atoms_with_charge += 1
                    
                    # Check for q00 charges
                    if hasattr(atom, 'q00') and atom.q00 != 0.0:
                        atoms_with_q00 += 1
                    
                    # Check for q11 charges
                    if hasattr(atom, 'q11') and atom.q11 != 0.0:
                        atoms_with_q11 += 1
    
    print(f"    Total atoms in structure: {total_atoms}")
    print(f"    Atoms with 'charge': {atoms_with_charge}")
    print(f"    Atoms with 'q00': {atoms_with_q00}")
    print(f"    Atoms with 'q11': {atoms_with_q11}")
    
    # =================================================================
    # METHOD 4: Helper functions for easy access
    # =================================================================
    print("\n" + "="*50)
    print("METHOD 4: Helper functions for easy access")
    print("="*50)
    
    def get_atom_q00(atom):
        """Helper function to get q00 charge from any atom."""
        # Try different attribute names
        for attr in ['q00', 'q_00']:
            if hasattr(atom, attr):
                return getattr(atom, attr)
        return None
    
    def get_atom_charge(atom):
        """Helper function to get any charge from atom."""
        # Try in order of preference
        for attr in ['charge', 'q00', 'q_00']:
            if hasattr(atom, attr) and getattr(atom, attr) != 0.0:
                return getattr(atom, attr)
        return 0.0
    
    def find_atoms_with_charges(structure, charge_type='q00'):
        """Find all atoms with specific charge type."""
        atoms_found = []
        attr_names = {
            'q00': ['q00', 'q_00'],
            'q11': ['q11', 'q_11'],
            'charge': ['charge']
        }
        
        for model in structure.pdb:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        for attr in attr_names.get(charge_type, [charge_type]):
                            if hasattr(atom, attr) and getattr(atom, attr) != 0.0:
                                atoms_found.append((atom, getattr(atom, attr)))
                                break
        return atoms_found
    
    print(f"\n4. Using helper functions:")
    
    # Example usage of helper functions
    if len(pigment_system.pigments) > 0:
        first_pigment = list(pigment_system.pigments.values())[0]
        first_atom = first_pigment.atoms[0]
        
        q00_value = get_atom_q00(first_atom)
        charge_value = get_atom_charge(first_atom)
        
        print(f"    First pigment atom: {first_atom.name}")
        print(f"    q00 value: {q00_value}")
        print(f"    charge value: {charge_value}")
    
    # Find all q00 atoms
    q00_atoms_all = find_atoms_with_charges(structure, 'q00')
    print(f"    Total atoms with q00 charges: {len(q00_atoms_all)}")
    
    # =================================================================
    # SUMMARY AND BEST PRACTICES
    # =================================================================
    print("\n" + "="*50)
    print("🎯 SUMMARY AND BEST PRACTICES")
    print("="*50)
    
    print(f"""
✅ RECOMMENDED APPROACHES:

1. For PIGMENT atoms (fastest, O(1) access):
   - Use: atom.q00, atom.q11, atom.q01 directly
   - Check: hasattr(atom, 'calculation_ready') and atom.calculation_ready
   - Method: pigment.get_site_energy_atoms_direct()

2. For PROTEIN background atoms:
   - Use: atom.charge (standard charge)
   - Check: hasattr(atom, 'charge') and atom.charge != 0.0
   - Method: pigment_system.get_protein_background_atoms()

3. For ANY atom (generic):
   - Use helper functions to try multiple attribute names
   - Check for both 'q00' and 'q_00' variants

⚡ PERFORMANCE NOTES:
   - Enhanced atoms (pigments): O(1) direct access
   - Protein atoms: Standard charge attribute access
   - Avoid dictionary lookups in tight loops

🔍 DEBUGGING TIPS:
   - Use hasattr() to check if charge attributes exist
   - Check atom.calculation_ready for enhanced atoms
   - Use dir(atom) to see all available attributes
""")

if __name__ == "__main__":
    demonstrate_charge_access()