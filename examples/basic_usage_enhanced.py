"""Example showing enhanced architecture usage."""

from pathlib import Path

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator

script_dir = Path(__file__).resolve().parent
pdb_file = script_dir.parent / "Alprotein" / "pdb" / "CP24_extended_most_occ_pH8.pdb"

# # Option 1: Auto-detect and load with default parameters
# structure = ProteinStructure.from_file_with_parameters(
#     pdb_file,
#     parameter_config={
#         'ChlorophyllA': 'CLA',
#         'ChlorophyllB': 'CHL'
#     },
#     name="MyProtein"
# )

# Option 2: Load structure first, add parameters later
structure = ProteinStructure.from_file_with_parameters(pdb_file, name="MyProtein")
config = {'ChlorophyllA': 'CLA'}
structure.enrich_with_parameters(config)

# Create system and run calculations
pigment_system = PigmentSystem(structure)
calculator = SiteEnergyCalculator()

print(f"Found {len(pigment_system.pigments)} pigments")
print(f"Using enhanced atoms: {pigment_system.enhanced_atoms_available}")

for i, pigment in enumerate(pigment_system.pigments.values()):
    energy = calculator.calculate_site_energy(pigment, pigment_system)
    print(f"Pigment {i} ({pigment.residue.resname}): {energy:.2f} cm⁻¹")
