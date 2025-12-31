#!/usr/bin/env python3
"""
Test: Verify that protein atoms have proper q00 charges.
"""

import sys
import os
sys.path.insert(0, '.')

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem

def test_protein_q00_charges():
    """Test that protein atoms have q00 charges."""
    
    print("🧪 Testing Protein Q00 Charges")
    print("=" * 40)
    
    # Load structure with enhanced parameters
    pdb_path = "Alprotein/pdb/CP24_extended_most_occ_pH8.pdb"
    config = {'ChlorophyllA': 'CLA'}
    
    structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "test")
    pigment_system = PigmentSystem(structure)
    
    print(f"✓ Structure loaded")
    print(f"✓ Pigments found: {len(pigment_system.pigments)}")
    
    # Get protein background atoms
    background_atoms = pigment_system.get_protein_background_atoms()
    print(f"✓ Background atoms: {len(background_atoms)}")
    
    if background_atoms:
        # Test first few atoms
        print(f"\n📊 Testing protein atom charges:")
        
        atoms_with_charge = 0
        atoms_with_q00 = 0
        atoms_with_both = 0
        charge_mismatch = 0
        
        for i, atom in enumerate(background_atoms[:10]):  # Test first 10
            residue = atom.get_parent()
            
            has_charge = hasattr(atom, 'charge') and atom.charge is not None and atom.charge != 0.0
            has_q00 = hasattr(atom, 'q00') and atom.q00 is not None and atom.q00 != 0.0
            
            if has_charge:
                atoms_with_charge += 1
            if has_q00:
                atoms_with_q00 += 1
            if has_charge and has_q00:
                atoms_with_both += 1
                # Check if they match (both are not None at this point)
                if abs(atom.charge - atom.q00) > 1e-6:
                    charge_mismatch += 1
            
            if i < 5:  # Show details for first 5
                print(f"  {residue.resname}_{residue.id[1]}_{atom.name}:")
                print(f"    charge: {getattr(atom, 'charge', 'None')}")
                print(f"    q00: {getattr(atom, 'q00', 'None')}")
                if has_charge and has_q00:
                    print(f"    match: {abs(atom.charge - atom.q00) < 1e-6}")
        
        print(f"\n📈 Summary (first 10 atoms):")
        print(f"  Atoms with 'charge': {atoms_with_charge}")
        print(f"  Atoms with 'q00': {atoms_with_q00}")
        print(f"  Atoms with both: {atoms_with_both}")
        print(f"  Charge mismatches: {charge_mismatch}")
        
        if atoms_with_both > 0 and charge_mismatch == 0:
            print(f"✅ All protein atoms have matching charge and q00!")
        elif atoms_with_q00 == 0:
            print(f"❌ No protein atoms have q00 charges")
        else:
            print(f"⚠️  Some protein atoms have mismatched charges")
    
    # Also test pigment atoms
    print(f"\n🧬 Testing pigment atom charges:")
    if len(pigment_system.pigments) > 0:
        pigment = list(pigment_system.pigments.values())[0]
        
        pigment_charge_count = 0
        pigment_q00_count = 0
        pigment_q11_count = 0
        
        for atom in pigment.atoms:
            if hasattr(atom, 'charge') and atom.charge is not None and atom.charge != 0.0:
                pigment_charge_count += 1
            if hasattr(atom, 'q00') and atom.q00 is not None and atom.q00 != 0.0:
                pigment_q00_count += 1
            if hasattr(atom, 'q11') and atom.q11 is not None and atom.q11 != 0.0:
                pigment_q11_count += 1
        
        print(f"  Pigment atoms with 'charge': {pigment_charge_count}")
        print(f"  Pigment atoms with 'q00': {pigment_q00_count}")
        print(f"  Pigment atoms with 'q11': {pigment_q11_count}")
        
        if pigment_q00_count > 0 and pigment_q11_count > 0:
            print(f"✅ Pigment atoms have proper q00/q11 charges!")
        else:
            print(f"❌ Pigment atoms missing q00/q11 charges")
    
    print(f"\n✅ Test completed!")

if __name__ == "__main__":
    test_protein_q00_charges()
