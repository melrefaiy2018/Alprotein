from Alprotein import ProteinStructure, PigmentSystem, HamiltonianCalculator
# from Alprotein import ChlorophyllA


# PDB_file = 'most_occ_pH1.pdb'
# PDB_file = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57E/outside_rotomers/most_occ/most_occ_pH1.pdb'
# PDB_file = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57D/most_occ/most_occ_pH3.pdb'
# PDB_file = "/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57D/most_occ/test_pdb/most_occ_pH1.pdb"
PDB_file = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/paper_data_2025/With_water/Q57D/most_occ/most_occ_pH1.pdb'
# Load structure and calculate
protein = ProteinStructure.from_file_extended(PDB_file, 'Complex')
system = PigmentSystem(protein)
system.add_pigments_by_residue(
                                resname="CLA",
                                tresp_dict_name="CLA_IPPC",
                                cdc_dict_name="CLA"
                                )

calculator = HamiltonianCalculator(system, 
                                   dielectric_cdc=2.5, 
                                   dielectric_tresp=1, 
                                   E_0a=14950, 
                                   E_0b=0, 
                                   f_val=0.72)
hamiltonian = calculator.construct_hamiltonian()

# Get the detailed CDC results
total_shifts, contributions = calculator.site_energy_calc.calculate_detailed_site_energy_contributions(system)

print(total_shifts)        # overall site energy shift for each pigment
print(contributions)       # breakdown by each contributing residue/pigment