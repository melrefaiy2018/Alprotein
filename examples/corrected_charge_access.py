#!/usr/bin/env python3
"""
CORRECTED: How to access q00 charges from protein and pigment atoms.
"""

import sys
import os
sys.path.insert(0, '.')

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem

def explain_protein_charges():
    """Explain the corrected charge access for protein atoms."""
    
    print("🧬 CORRECTED: Protein and Pigment Charge Access")
    print("=" * 50)
    
    # Load structure with enhanced parameters
    pdb_path = "Alprotein/pdb/CP24_extended_most_occ_pH8.pdb"
    config = {'ChlorophyllA': 'CLA'}
    
    structure = ProteinStructure.from_file_with_parameters(pdb_path, config, "demo")
    pigment_system = PigmentSystem(structure)
    
    print(f"✓ Structure loaded with enhanced parameters")
    
    # =================================================================
    # EXPLANATION: What charges exist where
    # =================================================================
    print(f"\n📚 EXPLANATION:")
    print(f"   🧬 PROTEIN atoms (non-pigment residues like ASN, GLU, etc.):")
    print(f"      - Have 'charge' attribute (from PDB file)")
    print(f"      - NOW ALSO have 'q00' attribute (same value as charge)")
    print(f"      - Do NOT have 'q11' (no excited state for protein)")
    print(f"")
    print(f"   🌿 PIGMENT atoms (CLA, CHL residues with enhanced parameters):")
    print(f"      - Have 'charge' attribute (ground state, same as q00)")
    print(f"      - Have 'q00' attribute (ground state charge)")
    print(f"      - Have 'q11' attribute (excited state charge)")
    print(f"      - May have 'q01' attribute (transition charge)")
    
    # =================================================================
    # EXAMPLE 1: Protein atoms
    # =================================================================
    print(f"\n" + "="*50)
    print(f"EXAMPLE 1: PROTEIN ATOMS")
    print(f"="*50)
    
    background_atoms = pigment_system.get_protein_background_atoms()
    
    if background_atoms:
        atom = background_atoms[0]
        residue = atom.get_parent()
        
        print(f"\nProtein atom example: {residue.resname}_{residue.id[1]}_{atom.name}")
        print(f"  atom.charge = {atom.charge:.6f}")
        print(f"  atom.q00 = {atom.q00:.6f}")
        print(f"  hasattr(atom, 'q11') = {hasattr(atom, 'q11')}")
        print(f"  Values match: {abs(atom.charge - atom.q00) < 1e-6}")
        
        print(f"\n✅ CORRECT WAY to get protein charges:")
        print(f"   # Both work and give the same result:")
        print(f"   protein_charge = atom.charge    # {atom.charge:.6f}")
        print(f"   protein_q00 = atom.q00          # {atom.q00:.6f}")
    
    # =================================================================
    # EXAMPLE 2: Pigment atoms
    # =================================================================
    print(f"\n" + "="*50)
    print(f"EXAMPLE 2: PIGMENT ATOMS")
    print(f"="*50)
    
    if len(pigment_system.pigments) > 0:
        pigment = list(pigment_system.pigments.values())[0]
        
        # Find first atom with calculation parameters
        enhanced_atom = None
        for atom in pigment.atoms:
            if hasattr(atom, 'calculation_ready') and atom.calculation_ready:
                enhanced_atom = atom
                break
        
        if enhanced_atom:
            residue = enhanced_atom.get_parent()
            print(f"\nPigment atom example: {residue.resname}_{residue.id[1]}_{enhanced_atom.name}")
            print(f"  atom.charge = {enhanced_atom.charge:.6f}")
            print(f"  atom.q00 = {enhanced_atom.q00:.6f}")
            print(f"  atom.q11 = {enhanced_atom.q11:.6f}")
            if hasattr(enhanced_atom, 'q01'):
                print(f"  atom.q01 = {enhanced_atom.q01:.6f}")
            print(f"  atom.calculation_ready = {enhanced_atom.calculation_ready}")
            
            print(f"\n✅ CORRECT WAY to get pigment charges:")
            print(f"   # All these work:")
            print(f"   ground_state = atom.q00         # {enhanced_atom.q00:.6f}")
            print(f"   excited_state = atom.q11        # {enhanced_atom.q11:.6f}")
            print(f"   current_charge = atom.charge    # {enhanced_atom.charge:.6f} (usually = q00)")
            print(f"   charge_diff = atom.q11 - atom.q00  # {enhanced_atom.q11 - enhanced_atom.q00:.6f}")
    
    # =================================================================
    # PRACTICAL EXAMPLES
    # =================================================================
    print(f"\n" + "="*50)
    print(f"PRACTICAL USAGE EXAMPLES")
    print(f"="*50)
    
    print(f"\n1. Calculate electrostatic interaction (pigment-protein):")
    print(f"""
def calculate_interaction(pigment_atom, protein_atom):
    # For pigment atom: use q11 - q00 (excitation charge change)
    delta_q = pigment_atom.q11 - pigment_atom.q00
    
    # For protein atom: use q00 (or charge, they're the same)
    protein_charge = protein_atom.q00  # or protein_atom.charge
    
    # Calculate interaction
    distance = np.linalg.norm(pigment_atom.coord - protein_atom.coord)
    interaction = (delta_q * protein_charge) / distance
    return interaction
""")
    
    print(f"\n2. Get all atoms with q00 charges:")
    print(f"""
def get_all_q00_atoms(structure):
    q00_atoms = []
    for model in structure.pdb:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if hasattr(atom, 'q00') and atom.q00 != 0.0:
                        q00_atoms.append((atom, atom.q00))
    return q00_atoms

# This will now include BOTH protein and pigment atoms!
""")
    
    print(f"\n3. Safe charge access function:")
    print(f"""
def get_ground_state_charge(atom):
    '''Get ground state charge from any atom.'''
    # Try q00 first (works for both protein and pigment)
    if hasattr(atom, 'q00'):
        return atom.q00
    # Fallback to charge (for older code)
    elif hasattr(atom, 'charge'):
        return atom.charge
    else:
        return 0.0
""")
    
    print(f"\n✅ Summary:")
    print(f"   🧬 Protein atoms: Use atom.q00 or atom.charge (same value)")
    print(f"   🌿 Pigment atoms: Use atom.q00, atom.q11, atom.q01 as needed")
    print(f"   🔧 Both types now have q00 for consistency!")

if __name__ == "__main__":
    explain_protein_charges()
