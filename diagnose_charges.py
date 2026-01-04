#!/usr/bin/env python3
"""Diagnose why some atoms don't have charges in extended file."""

from Alprotein import ProteinStructure

# Load the extended file
print("Loading extended file...")
protein = ProteinStructure.from_file_extended('examples/data/extended_most_occ_pH7_IsiA.pdb', 'test')

atoms = list(protein.pdb.get_atoms())
print(f"Total atoms loaded: {len(atoms)}")

# Check charges
atoms_with_charge_attr = [a for a in atoms if hasattr(a, 'charge')]
atoms_with_nonzero_charge = [a for a in atoms if hasattr(a, 'charge') and a.charge != 0]

print(f"Atoms with 'charge' attribute: {len(atoms_with_charge_attr)}")
print(f"Atoms with non-zero charge: {len(atoms_with_nonzero_charge)}")

# Find atoms WITHOUT charge attribute
atoms_without_charge = [a for a in atoms if not hasattr(a, 'charge')]
print(f"\nAtoms WITHOUT 'charge' attribute: {len(atoms_without_charge)}")

if atoms_without_charge:
    print("\nFirst 10 atoms without charge:")
    for i, atom in enumerate(atoms_without_charge[:10]):
        res = atom.get_parent()
        chain = res.get_parent()

        # Check if atom is in extended_atoms dict
        atom_key = (chain.id, res.id[1], atom.name)
        in_dict = atom_key in protein.extended_atoms

        print(f"  {i+1}. {atom.name} in {res.resname} {res.id[1]} chain {chain.id}")
        print(f"     Serial: {atom.serial_number}")
        print(f"     In extended_atoms: {in_dict}")
        if in_dict:
            ext_data = protein.extended_atoms[atom_key]
            print(f"     Extended data: {ext_data}")

# Check extended_atoms dictionary
print(f"\nExtended atoms dictionary size: {len(protein.extended_atoms)}")

# Sample some extended_atoms entries
print("\nSample extended_atoms entries:")
for i, (key, value) in enumerate(list(protein.extended_atoms.items())[:5]):
    print(f"  Key {key}: {value}")

# Check if there are atoms with charge=0.0 vs no charge
atoms_with_zero_charge = [a for a in atoms if hasattr(a, 'charge') and a.charge == 0.0]
print(f"\nAtoms with charge == 0.0: {len(atoms_with_zero_charge)}")

# Compare atom serial numbers
serials_with_charge = set(a.serial_number for a in atoms_with_charge_attr)
serials_without_charge = set(a.serial_number for a in atoms_without_charge)

print(f"\nSerial number ranges:")
if serials_with_charge:
    print(f"  With charge: {min(serials_with_charge)} - {max(serials_with_charge)}")
if serials_without_charge:
    print(f"  Without charge: {min(serials_without_charge)} - {max(serials_without_charge)}")
