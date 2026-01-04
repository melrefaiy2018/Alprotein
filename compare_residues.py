#!/usr/bin/env python3
"""Compare residue structures between the two files."""

from Alprotein import ProteinStructure

print("=" * 80)
print("Comparing HIS 31 residue structure")
print("=" * 80)

# Load both files
print("\nLoading extended file...")
protein_ext = ProteinStructure.from_file_extended('examples/data/extended_most_occ_pH7_IsiA.pdb', 'ext')

print("Loading most_occ file...")
protein_most = ProteinStructure.from_file_extended('examples/data/most_occ_pH7_IsiA.pdb', 'most')

# Find HIS 31 in both structures
print("\n" + "=" * 80)
print("Finding HIS 31 in both structures")
print("=" * 80)

# Extended file
residues_ext = list(protein_ext.pdb.get_residues())
his31_ext = None
for res in residues_ext:
    if res.resname == 'HIS' and res.id[1] == 31:
        his31_ext = res
        break

# Most_occ file
residues_most = list(protein_most.pdb.get_residues())
his31_most_list = []
for res in residues_most:
    if res.resname == 'HIS' and res.id[1] == 31:
        his31_most_list.append(res)

print(f"\nExtended file - HIS 31 residues found: {1 if his31_ext else 0}")
if his31_ext:
    atoms_ext = list(his31_ext.get_atoms())
    print(f"  Atoms in residue: {len(atoms_ext)}")
    print(f"  Residue ID: {his31_ext.id}")
    print(f"  Atom names: {[a.name for a in atoms_ext]}")
    print(f"  Identifiers: {set(getattr(a, 'identifier', 'N/A') for a in atoms_ext)}")

print(f"\nMost_occ file - HIS 31 residues found: {len(his31_most_list)}")
for i, res in enumerate(his31_most_list):
    atoms = list(res.get_atoms())
    print(f"  Residue {i+1}:")
    print(f"    Atoms in residue: {len(atoms)}")
    print(f"    Residue ID: {res.id}")
    print(f"    Atom names: {[a.name for a in atoms]}")
    if atoms:
        print(f"    Conformer IDs: {set(getattr(a, 'conformer_id', 'N/A') for a in atoms)}")
        print(f"    Identifiers: {set(getattr(a, 'identifier', 'N/A') for a in atoms)}")

# Count total residues
print("\n" + "=" * 80)
print("Total residue counts")
print("=" * 80)
print(f"Extended file: {len(residues_ext)} residues")
print(f"Most_occ file: {len(residues_most)} residues")

# Check atoms per residue distribution
print("\n" + "=" * 80)
print("Residues with >20 atoms (likely mixed conformers)")
print("=" * 80)

large_res_ext = [r for r in residues_ext if len(list(r.get_atoms())) > 20]
print(f"\nExtended file: {len(large_res_ext)} residues with >20 atoms")
for res in large_res_ext[:5]:
    atoms = list(res.get_atoms())
    identifiers = set(getattr(a, 'identifier', 'N/A') for a in atoms)
    print(f"  {res.resname} {res.id[1]}: {len(atoms)} atoms, identifiers: {identifiers}")

large_res_most = [r for r in residues_most if len(list(r.get_atoms())) > 20]
print(f"\nMost_occ file: {len(large_res_most)} residues with >20 atoms")
for res in large_res_most[:5]:
    atoms = list(res.get_atoms())
    conformers = set(getattr(a, 'conformer_id', 'N/A') for a in atoms)
    print(f"  {res.resname} {res.id[1]}: {len(atoms)} atoms, conformers: {conformers}")
