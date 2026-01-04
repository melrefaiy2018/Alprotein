#!/usr/bin/env python3
"""Check for duplicate atom keys due to conformers."""

from Alprotein import ProteinStructure

# Load the extended file
print("Loading extended file...")
protein = ProteinStructure.from_file_extended('examples/data/extended_most_occ_pH7_IsiA.pdb', 'test')

atoms = list(protein.pdb.get_atoms())
print(f"Total atoms loaded by BioPython: {len(atoms)}")
print(f"Extended atoms dictionary size: {len(protein.extended_atoms)}")

# Build key count from atoms
atom_key_counts = {}
for atom in atoms:
    res = atom.get_parent()
    chain = res.get_parent()
    key = (chain.id, res.id[1], atom.name)
    atom_key_counts[key] = atom_key_counts.get(key, 0) + 1

# Find duplicates
duplicate_keys = {k: v for k, v in atom_key_counts.items() if v > 1}
print(f"\nDuplicate atom keys: {len(duplicate_keys)}")

if duplicate_keys:
    print("\nFirst 10 duplicate keys:")
    for i, (key, count) in enumerate(list(duplicate_keys.items())[:10]):
        print(f"  {key}: {count} atoms")

# Check if extended_atoms has fewer entries due to overwrites
expected_entries = len(atoms)
actual_entries = len(protein.extended_atoms)
print(f"\nExpected extended_atoms entries: {expected_entries}")
print(f"Actual extended_atoms entries: {actual_entries}")
print(f"Difference (overwrites): {expected_entries - actual_entries}")

# Check atoms with zero charge
atoms_with_zero = [a for a in atoms if hasattr(a, 'charge') and a.charge == 0.0]
print(f"\nAtoms with charge=0.0: {len(atoms_with_zero)}")

if atoms_with_zero:
    print("\nFirst 10 atoms with zero charge:")
    for i, atom in enumerate(atoms_with_zero[:10]):
        res = atom.get_parent()
        chain = res.get_parent()
        key = (chain.id, res.id[1], atom.name)

        print(f"  {i+1}. {atom.name} in {res.resname} {res.id[1]} chain {chain.id}")
        print(f"     Serial: {atom.serial_number}")

        # Check extended_atoms
        if key in protein.extended_atoms:
            ext_data = protein.extended_atoms[key]
            print(f"     Extended charge: {ext_data.get('charge', 'N/A')}")
            print(f"     Identifier: {ext_data.get('identifier', 'N/A')}")

# Check q00 attribute
atoms_with_q00 = [a for a in atoms if hasattr(a, 'q00')]
atoms_with_nonzero_q00 = [a for a in atoms if hasattr(a, 'q00') and a.q00 != 0.0]
print(f"\nAtoms with q00 attribute: {len(atoms_with_q00)}")
print(f"Atoms with non-zero q00: {len(atoms_with_nonzero_q00)}")
