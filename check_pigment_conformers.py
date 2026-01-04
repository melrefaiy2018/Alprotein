#!/usr/bin/env python3
"""Check pigment conformer loading between files."""

from Alprotein import ProteinStructure

print("=" * 80)
print("Analyzing Pigment Conformer Loading")
print("=" * 80)

# Load both files
print("\nLoading extended file...")
protein_ext = ProteinStructure.from_file_extended('examples/data/extended_most_occ_pH7_IsiA.pdb', 'ext')

print("Loading most_occ file...")
protein_most = ProteinStructure.from_file_extended('examples/data/most_occ_pH7_IsiA.pdb', 'most')

# Check CLA 501 specifically
print("\n" + "=" * 80)
print("CLA 501 Residue Analysis")
print("=" * 80)

def get_cla_501(protein):
    for residue in protein.pdb.get_residues():
        if residue.resname == 'CLA' and residue.id[1] == 501:
            return residue
    return None

cla_ext = get_cla_501(protein_ext)
cla_most = get_cla_501(protein_most)

if cla_ext:
    atoms_ext = list(cla_ext.get_atoms())
    print(f"\nExtended file - CLA 501:")
    print(f"  Total atoms: {len(atoms_ext)}")
    print(f"  Residue ID: {cla_ext.id}")

    # Group by identifier
    by_identifier = {}
    for atom in atoms_ext:
        ident = getattr(atom, 'identifier', 'unknown')
        if ident not in by_identifier:
            by_identifier[ident] = []
        by_identifier[ident].append(atom.name)

    print(f"  Atoms grouped by identifier:")
    for ident, atom_names in sorted(by_identifier.items()):
        print(f"    '{ident}': {len(atom_names)} atoms - {atom_names[:5]}...")

if cla_most:
    atoms_most = list(cla_most.get_atoms())
    print(f"\nMost_occ file - CLA 501:")
    print(f"  Total atoms: {len(atoms_most)}")
    print(f"  Residue ID: {cla_most.id}")

    # Group by conformer_id
    by_conformer = {}
    for atom in atoms_most:
        conf = getattr(atom, 'conformer_id', 'unknown')
        if conf not in by_conformer:
            by_conformer[conf] = []
        by_conformer[conf].append(atom.name)

    print(f"  Atoms grouped by conformer_id:")
    for conf, atom_names in sorted(by_conformer.items()):
        print(f"    '{conf}': {len(atom_names)} atoms - {atom_names[:5]}...")

# Check if accessing MG atom gives same result
print("\n" + "=" * 80)
print("Accessing MG atom from both residues")
print("=" * 80)

if cla_ext:
    mg_ext = cla_ext['MG']
    print(f"\nExtended file - MG atom:")
    print(f"  Coordinates: {mg_ext.coord}")
    print(f"  Identifier: {getattr(mg_ext, 'identifier', 'N/A')}")
    print(f"  Serial: {mg_ext.serial_number}")

if cla_most:
    mg_most = cla_most['MG']
    print(f"\nMost_occ file - MG atom:")
    print(f"  Coordinates: {mg_most.coord}")
    print(f"  Conformer ID: {getattr(mg_most, 'conformer_id', 'N/A')}")
    print(f"  Identifier: {getattr(mg_most, 'identifier', 'N/A')}")
    print(f"  Serial: {mg_most.serial_number}")

# Check coordinate difference
if cla_ext and cla_most:
    import numpy as np
    coord_diff = np.linalg.norm(mg_ext.coord - mg_most.coord)
    print(f"\nCoordinate difference for MG: {coord_diff:.6f} Å")

# Count all pigment residues
print("\n" + "=" * 80)
print("All Pigment Residues")
print("=" * 80)

def count_pigment_residues(protein):
    pigment_residues = [r for r in protein.pdb.get_residues() if r.resname == 'CLA']

    total_atoms = sum(len(list(r.get_atoms())) for r in pigment_residues)

    # Check for residues with many atoms (likely mixed conformers)
    large_residues = [r for r in pigment_residues if len(list(r.get_atoms())) > 100]

    return len(pigment_residues), total_atoms, len(large_residues)

num_res_ext, atoms_ext, large_ext = count_pigment_residues(protein_ext)
num_res_most, atoms_most, large_most = count_pigment_residues(protein_most)

print(f"\nExtended file:")
print(f"  Pigment residues: {num_res_ext}")
print(f"  Total pigment atoms: {atoms_ext}")
print(f"  Residues with >100 atoms (mixed conformers): {large_ext}")

print(f"\nMost_occ file:")
print(f"  Pigment residues: {num_res_most}")
print(f"  Total pigment atoms: {atoms_most}")
print(f"  Residues with >100 atoms (mixed conformers): {large_most}")
