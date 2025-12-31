#!/usr/bin/env python3
"""
Simple example: How to get q00 charges from atoms in different contexts.
"""

import sys
from pathlib import Path

sys.path.insert(0, '.')

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem

def simple_charge_access_example():
    """Simple example of accessing q00 charges."""
    
    print("Simple Charge Access Example")
    print("=" * 40)
    
    # Load structure with enhanced parameters
    script_dir = Path(__file__).resolve().parent
    pdb_path = script_dir / "Alprotein" / "pdb" / "CP24_extended_most_occ_pH8.pdb"
    config = {'ChlorophyllA': 'CLA'}
    
    structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "demo")
    pigment_system = PigmentSystem(structure)
    
    # =================================================================
    # EXAMPLE 1: Get q00 from pigment atoms (ENHANCED)
    # =================================================================
    print("\n1. Getting q00 from PIGMENT atoms:")
    
    if len(pigment_system.pigments) > 0:
        pigment = list(pigment_system.pigments.values())[0]
        print(f"   Pigment: {pigment.name}")
        
        # Method A: Direct access (FASTEST)
        print(f"\n   Method A - Direct access:")
        for atom in pigment.atoms:
            if hasattr(atom, 'calculation_ready') and atom.calculation_ready:
                print(f"     {atom.name}: q00 = {atom.q00:.6f}")
                break  # Just show first one
        
        # Method B: Using pigment method
        print(f"\n   Method B - Using pigment method:")
        site_atoms = pigment.get_site_energy_atoms()
        if site_atoms:
            atom_name, q00, q11 = site_atoms[0]
            print(f"     {atom_name}: q00 = {q00:.6f}, q11 = {q11:.6f}")
        
        # Method C: Get specific atom
        print(f"\n   Method C - Get specific atom:")
        try:
            mg_atom = pigment.residue['MG']
            if hasattr(mg_atom, 'q00'):
                print(f"     MG: q00 = {mg_atom.q00:.6f}")
            else:
                print(f"     MG: no q00 charge")
        except KeyError:
            print(f"     No MG atom found")
    
    # =================================================================
    # EXAMPLE 2: Get charges from protein atoms
    # =================================================================
    print(f"\n2. Getting charges from PROTEIN atoms:")
    
    background_atoms = pigment_system.get_protein_background_atoms()
    print(f"   Background atoms: {len(background_atoms)}")
    
    if background_atoms:
        atom = background_atoms[0]
        residue = atom.get_parent()
        print(f"   First background atom: {residue.resname}_{residue.id[1]}_{atom.name}")
        
        # Check for different charge types
        if hasattr(atom, 'charge'):
            print(f"     charge = {atom.charge:.6f}")
        
        if hasattr(atom, 'q00'):
            print(f"     q00 = {atom.q00:.6f}")
        else:
            print(f"     no q00 charge (normal for protein atoms)")
    
    # =================================================================
    # EXAMPLE 3: Helper functions for any atom
    # =================================================================
    print(f"\n3. Helper functions for ANY atom:")
    
    def get_q00_safe(atom):
        """Safely get q00 from any atom."""
        for attr in ['q00', 'q_00']:
            if hasattr(atom, attr):
                return getattr(atom, attr)
        return None
    
    def get_any_charge(atom):
        """Get any available charge from atom."""
        for attr in ['q00', 'q_00', 'charge']:
            if hasattr(atom, attr):
                value = getattr(atom, attr)
                if value != 0.0:
                    return attr, value
        return None, 0.0
    
    # Test helper functions
    if len(pigment_system.pigments) > 0:
        pigment = list(pigment_system.pigments.values())[0]
        test_atom = pigment.atoms[0]
        
        q00_val = get_q00_safe(test_atom)
        charge_type, charge_val = get_any_charge(test_atom)
        
        print(f"   Test atom: {test_atom.name}")
        print(f"     q00 (safe): {q00_val}")
        print(f"     any charge: {charge_type} = {charge_val}")
    
    # =================================================================
    # EXAMPLE 4: Practical usage patterns
    # =================================================================
    print(f"\n4. Practical usage patterns:")
    
    # Pattern A: Calculate electrostatic interaction
    print(f"\n   Pattern A - Calculate interaction energy:")
    if len(pigment_system.pigments) > 0:
        pigment = list(pigment_system.pigments.values())[0]
        background_atoms = pigment_system.get_protein_background_atoms()
        
        total_interaction = 0.0
        interactions_count = 0
        
        # For pigment atoms with q00 charges
        for pig_atom in pigment.atoms:
            if hasattr(pig_atom, 'calculation_ready') and pig_atom.calculation_ready:
                # Get charge difference (excitation)
                delta_q = pig_atom.q11 - pig_atom.q00
                
                # Interact with first few background atoms (example)
                for bg_atom in background_atoms[:5]:
                    if hasattr(bg_atom, 'charge') and bg_atom.charge != 0.0:
                        # Calculate distance
                        import numpy as np
                        distance = np.linalg.norm(pig_atom.coord - bg_atom.coord)
                        
                        if distance > 0.1:  # Avoid division by zero
                            interaction = (delta_q * bg_atom.charge) / distance
                            total_interaction += interaction
                            interactions_count += 1
                
                break  # Just one pigment atom for demo
        
        print(f"     Sample interaction energy: {total_interaction:.6f}")
        print(f"     Interactions calculated: {interactions_count}")
    
    # Pattern B: Find all atoms with specific charge
    print(f"\n   Pattern B - Find atoms with q00 > 0.5:")
    high_charge_atoms = []
    
    for model in structure.pdb:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if hasattr(atom, 'q00') and atom.q00 > 0.5:
                        high_charge_atoms.append((atom, atom.q00))
    
    print(f"     Found {len(high_charge_atoms)} atoms with q00 > 0.5")
    if high_charge_atoms:
        atom, charge = high_charge_atoms[0]
        residue = atom.get_parent()
        print(f"     Example: {residue.resname}_{residue.id[1]}_{atom.name} = {charge:.6f}")
    
    print(f"\n✅ Examples completed!")

if __name__ == "__main__":
    simple_charge_access_example()
