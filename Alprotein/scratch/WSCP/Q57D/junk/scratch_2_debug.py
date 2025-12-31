#!/usr/bin/env python3
"""
Debug version of scratch_2.py to investigate missing protein contributions
"""

import os
import sys
from pathlib import Path

# Add proper imports
from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.utils.atom_level_cdc import (
    analyze_site_energy_contributions_with_atoms,
    find_strongest_atomic_contribution,
    analyze_backbone_vs_sidechain
)

def debug_protein_structure(protein):
    """Debug function to examine protein structure"""
    print("\n🔍 DEBUGGING PROTEIN STRUCTURE")
    print("=" * 50)
    
    # Get all residues
    residues = list(protein.pdb.get_residues())
    print(f"Total residues in structure: {len(residues)}")
    
    # Categorize residues
    pigment_residues = []
    protein_residues = []
    water_residues = []
    other_residues = []
    
    for residue in residues:
        resname = residue.get_resname()
        if resname in ['CLA', 'CHL']:
            pigment_residues.append(resname)
        elif resname == 'HOH':
            water_residues.append(resname)
        elif resname in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
                         'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
            protein_residues.append(resname)
        else:
            other_residues.append(resname)
    
    print(f"Pigment residues: {len(pigment_residues)} ({set(pigment_residues)})")
    print(f"Protein residues: {len(protein_residues)} (first 10: {protein_residues[:10]})")
    print(f"Water residues: {len(water_residues)}")
    print(f"Other residues: {len(other_residues)} ({set(other_residues)})")
    
    return len(protein_residues) > 0

def debug_pigment_system(system):
    """Debug function to examine pigment system"""
    print("\n🔍 DEBUGGING PIGMENT SYSTEM")
    print("=" * 50)
    
    print(f"Number of pigments: {len(system.pigments)}")
    for i, pigment in enumerate(system.pigments):
        print(f"Pigment {i+1}: {pigment.residue.get_resname()}_{pigment.residue.get_id()[1]}")
    
    # Check if pigments have CDC dictionaries loaded
    for pigment in system.pigments:
        if hasattr(pigment, 'cdc_dict') and pigment.cdc_dict:
            print(f"✅ Pigment {pigment.residue.get_resname()}_{pigment.residue.get_id()[1]} has CDC dictionary")
        else:
            print(f"❌ Pigment {pigment.residue.get_resname()}_{pigment.residue.get_id()[1]} missing CDC dictionary")

def debug_calculator_setup(system, calculator):
    """Debug the calculator setup"""
    print("\n🔍 DEBUGGING CALCULATOR SETUP")
    print("=" * 50)
    
    print(f"Calculator type: {type(calculator)}")
    print(f"Dielectric constant: {calculator.dielectric}")
    
    # Try to see what the calculator can access
    print(f"System has {len(system.pigments)} pigments")
    print(f"Protein structure has {len(list(system.protein.pdb.get_residues()))} residues")

# Main execution
def main():
    # File paths
    PDB_file = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57E/outside_rotomers/most_occ/most_occ_pH1.pdb'
    
    print(f"🔄 Loading PDB file: {PDB_file}")
    
    # Load structure and calculate
    protein = ProteinStructure.from_file(PDB_file, 'Complex')
    has_protein_residues = debug_protein_structure(protein)
    
    if not has_protein_residues:
        print("❌ ERROR: No protein residues found!")
        return
    
    # Setup pigment system
    system = PigmentSystem(protein)
    pigment_count = system.add_pigments_by_residue(
        resname="CLA",
        tresp_dict_name="CLA_IPPC",
        cdc_dict_name="CLA"
    )
    
    print(f"✅ Added {pigment_count} CLA pigments")
    debug_pigment_system(system)
    
    # Setup calculator
    calculator = SiteEnergyCalculator(dielectric=2.0)
    debug_calculator_setup(system, calculator)
    
    print(f"\n🔄 Analyzing site energy contributions...")
    print(f"Working directory: {os.getcwd()}")
    print(f"Cutoff distance: 20 Å")
    
    # This is the key function call that should include protein contributions
    try:
        results = analyze_site_energy_contributions_with_atoms(
            system, 
            calculator, 
            str(os.getcwd()), 
            cutoff=20
        )
        
        print("✅ Analysis completed successfully")
        
        # Let's examine the results structure
        print(f"\nResults keys: {list(results.keys())}")
        
        for pigment_id, contributions in results.items():
            print(f"\n📊 Contributions to {pigment_id}:")
            for source, data in contributions.items():
                if source != 'vacuum':
                    print(f"  {source}: {data['total_contribution']:.2f} cm⁻¹")
            
            # Count non-pigment, non-vacuum contributions
            protein_contributions = {k: v for k, v in contributions.items() 
                                   if k not in ['vacuum'] and not k.endswith('_CLA_1001')}
            
            if protein_contributions:
                print(f"  ✅ Found {len(protein_contributions)} protein contributions")
            else:
                print(f"  ❌ No protein contributions found!")
        
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
