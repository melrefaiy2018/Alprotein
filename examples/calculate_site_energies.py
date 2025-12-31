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

def main():
    """Main function to run the site energy calculation."""
    print("🧬 Alprotein Site Energy Calculator")
    print("=" * 50)
    
    # Determine PDB file to use
    if len(sys.argv) > 1:
        pdb_path = sys.argv[1]
    else:
        # Use default test file
        pdb_path = "Alprotein/pdb/CP24_extended_most_occ_pH8.pdb"
        print(f"💡 No PDB file specified, using default: {pdb_path}")
    
    # Check if file exists
    if not os.path.exists(pdb_path):
        print(f"❌ Error: PDB file not found: {pdb_path}")
        print(f"💡 Usage: python {sys.argv[0]} [pdb_file]")
        sys.exit(1)
    
    # Calculate site energies using enhanced architecture
    site_energies = calculate_site_energies(pdb_path, use_enhanced=True)
    
    # Optional: Save results to file
    if site_energies:
        output_file = f"{Path(pdb_path).stem}_site_energies.txt"
        with open(output_file, 'w') as f:
            f.write("# Site Energy Calculation Results\n")
            f.write(f"# PDB file: {pdb_path}\n")
            f.write("# Pigment_Name\tSite_Energy_cm-1\n")
            
            for name, energy in site_energies.items():
                if energy is not None:
                    f.write(f"{name}\t{energy:.6f}\n")
                else:
                    f.write(f"{name}\tN/A\n")
        
        print(f"\n💾 Results saved to: {output_file}")
    
    print("\n✨ Calculation complete!")

if __name__ == "__main__":
    main()
