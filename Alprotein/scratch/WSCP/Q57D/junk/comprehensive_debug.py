#!/usr/bin/env python3
"""
Fixed version with comprehensive debugging and multiple approaches
"""

import os
import sys
from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.utils.atom_level_cdc import (
    analyze_site_energy_contributions_with_atoms,
    find_strongest_atomic_contribution,
    analyze_backbone_vs_sidechain
)

def try_manual_cdc_calculation(system, calculator, cutoff=20):
    """Manual CDC calculation to debug missing protein contributions"""
    print(f"\n🔄 Manual CDC Analysis (cutoff={cutoff} Å)")
    print("=" * 50)
    
    for i, pigment in enumerate(system.pigments):
        # pigment_id = f"{pigment.residue.get_id()[0]}_{pigment.residue.get_resname()}_{pigment.residue.get_id()[1]}"
        # print(f"\nAnalyzing contributions to {pigment_id}:")
        
        # Get pigment center
        pigment_atoms = list(pigment.residue.get_atoms())
        pigment_center = sum([atom.get_coord() for atom in pigment_atoms]) / len(pigment_atoms)
        
        protein_contributions = []
        
        # Check all other residues
        for residue in system.protein.pdb.get_residues():
            if residue == pigment.residue:
                continue
                
            # Skip other pigments for now
            if residue.get_resname() in ['CLA', 'CHL']:
                continue
                
            # Calculate distance
            residue_atoms = list(residue.get_atoms())
            if not residue_atoms:
                continue
                
            residue_center = sum([atom.get_coord() for atom in residue_atoms]) / len(residue_atoms)
            distance = ((pigment_center - residue_center) ** 2).sum() ** 0.5
            
            if distance <= cutoff:
                res_id = f"{residue.get_id()[0]}_{residue.get_resname()}_{residue.get_id()[1]}"
                protein_contributions.append((res_id, distance))
        
        print(f"  Found {len(protein_contributions)} protein residues within {cutoff} Å")
        
        # Show closest 5
        protein_contributions.sort(key=lambda x: x[1])
        for res_id, dist in protein_contributions[:5]:
            print(f"    {res_id}: {dist:.2f} Å")
        
        if len(protein_contributions) == 0:
            print(f"  ❌ No protein residues found within {cutoff} Å!")
            print(f"  Try increasing cutoff distance...")

# Main execution with multiple approaches
PDB_file = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57E/outside_rotomers/most_occ/most_occ_pH1.pdb'

print("🔧 APPROACH 1: Fixed Original Script")
print("=" * 60)

# Load structure
protein = ProteinStructure.from_file(PDB_file, 'Complex')
system = PigmentSystem(protein)
system.add_pigments_by_residue(
    resname="CLA",
    tresp_dict_name="CLA_IPPC",
    cdc_dict_name="CLA"
)

calculator = SiteEnergyCalculator(dielectric=2.0)

# Try with different cutoff distances
for cutoff in [1]:
    print(f"\n🔄 Trying cutoff = {cutoff} Å")
    try:
        results = analyze_site_energy_contributions_with_atoms(
            system, calculator, str(os.getcwd()), cutoff=cutoff
        )
        
        # Count non-pigment contributions
        total_protein_contributions = 0
        for pigment_id, contributions in results.items():
            protein_contribs = {k: v for k, v in contributions.items() 
                              if k not in ['vacuum'] and not k.endswith('_CLA_1001')}
            total_protein_contributions += len(protein_contribs)
            
            if protein_contribs:
                print(f"  ✅ {pigment_id}: {len(protein_contribs)} protein contributions")
                # Show top 3 contributions
                sorted_contribs = sorted(protein_contribs.items(), 
                                       key=lambda x: abs(x[1]['total_contribution']), 
                                       reverse=True)
                for contrib_id, data in sorted_contribs[:3]:
                    print(f"    {contrib_id}: {data['total_contribution']:.2f} cm⁻¹")
        
        if total_protein_contributions == 0:
            print(f"  ❌ Still no protein contributions with cutoff {cutoff} Å")
        else:
            print(f"  ✅ Found {total_protein_contributions} total protein contributions!")
            break
            
    except Exception as e:
        print(f"  ❌ Error with cutoff {cutoff}: {e}")

# Manual analysis to understand the issue
try_manual_cdc_calculation(system, calculator, cutoff=30)

print(f"\n🔧 APPROACH 2: Alternative CDC Calculation")
print("=" * 60)

# Try direct site energy calculation
try:
    for i, pigment in enumerate(system.pigments):
        pigment_id = f"{pigment.residue.get_id()[0]}_{pigment.residue.get_resname()}_{pigment.residue.get_id()[1]}"
        
        # Try to calculate site energy directly
        site_energy = calculator.calculate_site_energy(pigment, system.protein, system.pigments)
        print(f"{pigment_id} site energy: {site_energy:.2f} cm⁻¹")
        
except Exception as e:
    print(f"Direct calculation failed: {e}")

print(f"\n📋 RECOMMENDATIONS:")
print("=" * 60)
print("1. ✅ Added missing 'import os' - this was causing the original error")
print("2. 🔍 Try running the debug script to see detailed analysis")
print("3. 📈 If still no protein contributions, the issue might be:")
print("   • Missing CDC parameters for amino acid residues")
print("   • Bug in the analyze_site_energy_contributions_with_atoms function")
print("   • Need to explicitly load amino acid CDC dictionaries")
print("4. 🧪 Try the manual CDC calculation approach above")
print("5. 📞 Contact the Alprotein developers if the issue persists")
