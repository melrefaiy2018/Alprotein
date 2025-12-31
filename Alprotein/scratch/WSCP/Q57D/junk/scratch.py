
from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.utils.atom_level_cdc import (
    analyze_site_energy_contributions_with_atoms,
    find_strongest_atomic_contribution,
    analyze_backbone_vs_sidechain
)


# # Here, is the pymembrane version
# dir_path = f'{os.getcwd()}/'
# # load the pdb file and the hamiltonian:
# # =====================================
# pdb_file = f'{dir_path}extended_most_occ_pH8.pdb'

# # build atomic protein:
# # ======================
# pdb_atomic = PDB.PDBParser().get_structure('CP24', pdb_file)
# cp24_atomic = ElectrostaticProteinAtomic(pdb_atomic, path_extended_pdb=pdb_file, name='CP24')

# # define the protein parameters:
# # ==============================
# cp24_atomic.fval_tresp = 0.72
# cp24_atomic.dielectric_tresp = 1
# cp24_atomic.dielectric_dipole = 1.4
# cp24_atomic.dielectric_cdc = 2

# cp24_atomic.prepare_pigments('CLA', ChlorophyllAtomic,
#                         list_chain=None,
#                         lineshape=kubo_renger_cla,
#                         dict_coupling_data=coupling_renger_by_type_mcce['CLA_IPPC'],
#                         dict_q0011_charges=q0011_renger_by_type_mcce['CLA'],
#                         disorder=40)

# cp24_atomic.prepare_pigments('CHL', ChlorophyllAtomic,
#                         list_chain=None,
#                         lineshape=kubo_renger_cla,
#                         dict_coupling_data=coupling_renger_by_type_mcce['CHL_IPPC'],
#                         dict_q0011_charges=q0011_renger_by_type_mcce['CHL'],
#                         disorder=30)

# # calculate site energy shift:
# # ============================
# e0_a = 14900
# e0_b = 15350 
# N_ens = 100
# temp = 300

# dict_total_site_energy_shift, dict_total_contribution_site_energy = cp24_atomic.calculate_cdc_site_shift(verbose=False)
# dict_site_energy_raw = cp24_atomic.calculate_total_site_energy(dict_total_site_energy_shift, E_0a=e0_a, E_0b=e0_b)
# print(dict_site_energy_raw)




# here is the Alprotein version:
# ==============================
pdb_file = script_dir.parent / "Alprotein" / "pdb" / "CP24_most_occ_pH7.pdb"
cutoff = 20.0
if not pdb_file.exists():
    logger.error(f"Missing PDB file: {pdb_file}")
    return

protein = ProteinStructure.from_file_extended(str(pdb_file), 'atom_analysis')
pigment_system = PigmentSystem(protein)
pigment_system.add_pigments_by_residue(
                resname="CLA",
                tresp_dict_name="CLA_IPPC",
                cdc_dict_name="CLA"
            )
# pigment_system.add_pigments_by_residue(
#     resname="CHL",
#     tresp_dict_name="CHL_IPPC",
#     cdc_dict_name="CHL"
# )
calculator = SiteEnergyCalculator(dielectric=2.0)
print(calculator)
results = analyze_site_energy_contributions_with_atoms(
    pigment_system, calculator, str(calc_dir), cutoff=cutoff)

logger.info("\n✅ Analysis completed. Generating outputs...")

# Save structured output
json_file = calc_dir / "cdc_contributions_with_atoms.json"
txt_file = calc_dir / "cdc_summary_with_atoms.txt"
pd.DataFrame(results).to_json(json_file, indent=2)
with open(txt_file, 'w') as f:
    f.write("CDC Summary Table\n")
    f.write("="*60 + "\n")
    for pid, pdata in results.items():
        for rid, rdata in pdata.items():
            f.write(f"{pid}-{rid}: {rdata.get('total_contribution', 0):+.2f} cm^-1\n")

# Generate bar plot of dominant atoms
plot_dominant_atom_bar(results, out_dir, cutoff)
logger.info(f"Figures saved to: {out_dir}")