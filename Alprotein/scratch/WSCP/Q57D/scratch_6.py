
import sys
import os
from pathlib import Path

# # Add current directory to path if needed
# if '.' not in sys.path:
#     sys.path.insert(0, '.')

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.data.parameters import create_parameter_config



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
calculator = SiteEnergyCalculator(dielectric_constant=2.5, e0a=14950, e0b=0)
site_energies, contributions = calculator.calculate_detailed_site_energy_contributions(pigment_system)
vacuum_energies = {'CLA': 14950}
calculator.calculate_all_site_energies(pigment_system, vacuum_energies)