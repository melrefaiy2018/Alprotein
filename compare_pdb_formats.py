#!/usr/bin/env python3
"""
Compare PDB Format 1 and Format 2 parsing to verify correctness.

This script loads both formats and compares all parsed fields to identify
any discrepancies in atom properties, charges, coordinates, etc.
"""

import numpy as np
from Alprotein import ProteinStructure


def compare_atoms(atom1, atom2, atom_idx):
    """Compare two atoms and return differences."""
    differences = []

    # Basic attributes
    if atom1.name != atom2.name:
        differences.append(f"  name: {atom1.name!r} vs {atom2.name!r}")

    if atom1.element != atom2.element:
        differences.append(f"  element: {atom1.element!r} vs {atom2.element!r}")

    if atom1.serial_number != atom2.serial_number:
        differences.append(f"  serial_number: {atom1.serial_number} vs {atom2.serial_number}")

    # Coordinates (with small tolerance for floating point)
    coord_diff = np.linalg.norm(atom1.coord - atom2.coord)
    if coord_diff > 1e-6:
        differences.append(f"  coord: {atom1.coord} vs {atom2.coord} (diff: {coord_diff:.6f})")

    # Charges
    charge1 = getattr(atom1, 'charge', None)
    charge2 = getattr(atom2, 'charge', None)
    if charge1 != charge2:
        # Allow small floating point differences
        if charge1 is None or charge2 is None or abs(charge1 - charge2) > 1e-6:
            differences.append(f"  charge: {charge1} vs {charge2}")

    # Identifier
    id1 = getattr(atom1, 'identifier', None)
    id2 = getattr(atom2, 'identifier', None)
    if id1 != id2:
        differences.append(f"  identifier: {id1!r} vs {id2!r}")

    # Conformer ID
    conf1 = getattr(atom1, 'conformer_id', None)
    conf2 = getattr(atom2, 'conformer_id', None)
    if conf1 != conf2:
        differences.append(f"  conformer_id: {conf1!r} vs {conf2!r}")

    # Occupancy
    if atom1.occupancy != atom2.occupancy:
        # Allow small floating point differences
        if atom1.occupancy is None or atom2.occupancy is None or abs(atom1.occupancy - atom2.occupancy) > 1e-6:
            differences.append(f"  occupancy: {atom1.occupancy} vs {atom2.occupancy}")

    # B-factor
    if atom1.bfactor != atom2.bfactor:
        if atom1.bfactor is None or atom2.bfactor is None or abs(atom1.bfactor - atom2.bfactor) > 1e-6:
            differences.append(f"  bfactor: {atom1.bfactor} vs {atom2.bfactor}")

    # Residue information
    res1 = atom1.get_parent()
    res2 = atom2.get_parent()

    if res1.resname != res2.resname:
        differences.append(f"  resname: {res1.resname!r} vs {res2.resname!r}")

    if res1.id[1] != res2.id[1]:  # residue number
        differences.append(f"  resseq: {res1.id[1]} vs {res2.id[1]}")

    chain1 = res1.get_parent().id
    chain2 = res2.get_parent().id
    if chain1 != chain2:
        differences.append(f"  chain_id: {chain1!r} vs {chain2!r}")

    return differences


def main():
    print("=" * 80)
    print("PDB Format Comparison: Format 1 (MCCE) vs Format 2 (Standard Variable)")
    print("=" * 80)

    # File paths
    format1_file = '/Users/mohamed/Documents/Research/Projects/GitHub_repo/AlProtein/examples/data/extended_most_occ_pH7_IsiA.pdb'
    format2_file = '/Users/mohamed/Documents/Research/Projects/GitHub_repo/AlProtein/examples/data/most_occ_pH7_IsiA.pdb'

    print(f"\nFormat 1 file: {format1_file}")
    print(f"Format 2 file: {format2_file}\n")

    # Load both files
    print("Loading Format 1...")
    protein1 = ProteinStructure.from_file_extended(format1_file, 'format1')
    atoms1 = list(protein1.pdb.get_atoms())

    print("Loading Format 2...")
    protein2 = ProteinStructure.from_file_extended(format2_file, 'format2')
    atoms2 = list(protein2.pdb.get_atoms())

    # Normalize metadata fields that are format-specific
    protein1.normalize_metadata()
    protein2.normalize_metadata()

    print("\n" + "=" * 80)
    print("OVERALL STATISTICS")
    print("=" * 80)

    print(f"\nTotal atoms:")
    print(f"  Format 1: {len(atoms1)}")
    print(f"  Format 2: {len(atoms2)}")
    print(f"  Difference: {abs(len(atoms1) - len(atoms2))}")

    # Charged atoms
    charged1 = protein1.get_atoms_with_charges()
    charged2 = protein2.get_atoms_with_charges()
    print(f"\nAtoms with charges:")
    print(f"  Format 1: {len(charged1)} ({100*len(charged1)/len(atoms1):.1f}%)")
    print(f"  Format 2: {len(charged2)} ({100*len(charged2)/len(atoms2):.1f}%)")

    # Extended atoms
    print(f"\nExtended atoms dictionary:")
    print(f"  Format 1: {len(protein1.extended_atoms)} entries")
    print(f"  Format 2: {len(protein2.extended_atoms)} entries")

    # Compare conformer information
    conformers1 = protein1.get_conformer_info()
    conformers2 = protein2.get_conformer_info()
    print(f"\nConformers:")
    print(f"  Format 1: {len(conformers1)} conformers")
    if conformers1:
        for conf_id, atoms in list(conformers1.items())[:3]:
            print(f"    {conf_id}: {len(atoms)} atoms")
    print(f"  Format 2: {len(conformers2)} conformers")
    if conformers2:
        for conf_id, atoms in list(conformers2.items())[:3]:
            print(f"    {conf_id}: {len(atoms)} atoms")

    print("\n" + "=" * 80)
    print("DETAILED ATOM COMPARISON (First 20 matching atoms)")
    print("=" * 80)

    # Compare first N atoms that can be matched by name
    num_to_compare = min(20, len(atoms1), len(atoms2))

    print(f"\nComparing first {num_to_compare} atoms...\n")

    total_differences = 0
    atoms_with_differences = 0

    for i in range(num_to_compare):
        atom1 = atoms1[i]
        atom2 = atoms2[i]

        differences = compare_atoms(atom1, atom2, i)

        if differences:
            atoms_with_differences += 1
            total_differences += len(differences)

            res1 = atom1.get_parent()
            print(f"Atom {i}: {atom1.name} in {res1.resname} {res1.id[1]}")
            for diff in differences:
                print(diff)
            print()

    if atoms_with_differences == 0:
        print("✅ All compared atoms are IDENTICAL!\n")
    else:
        print(f"⚠️  Found differences in {atoms_with_differences}/{num_to_compare} atoms")
        print(f"   Total field differences: {total_differences}\n")

    print("=" * 80)
    print("SAMPLE ATOM DETAILS")
    print("=" * 80)

    # Show detailed info for first atom from each format
    print("\n--- Format 1 Sample (Atom 0) ---")
    atom1 = atoms1[0]
    res1 = atom1.get_parent()
    print(f"Name: {atom1.name}")
    print(f"Element: {atom1.element}")
    print(f"Serial: {atom1.serial_number}")
    print(f"Coordinates: {atom1.coord}")
    print(f"Charge: {getattr(atom1, 'charge', 'N/A')}")
    print(f"Identifier: {getattr(atom1, 'identifier', 'N/A')}")
    print(f"Conformer ID: {getattr(atom1, 'conformer_id', 'N/A')}")
    print(f"Occupancy: {atom1.occupancy}")
    print(f"B-factor: {atom1.bfactor}")
    print(f"Residue: {res1.resname} {res1.id}")
    print(f"Chain: {res1.get_parent().id}")

    print("\n--- Format 2 Sample (Atom 0) ---")
    atom2 = atoms2[0]
    res2 = atom2.get_parent()
    print(f"Name: {atom2.name}")
    print(f"Element: {atom2.element}")
    print(f"Serial: {atom2.serial_number}")
    print(f"Coordinates: {atom2.coord}")
    print(f"Charge: {getattr(atom2, 'charge', 'N/A')}")
    print(f"Identifier: {getattr(atom2, 'identifier', 'N/A')}")
    print(f"Conformer ID: {getattr(atom2, 'conformer_id', 'N/A')}")
    print(f"Occupancy: {atom2.occupancy}")
    print(f"B-factor: {atom2.bfactor}")
    print(f"Residue: {res2.resname} {res2.id}")
    print(f"Chain: {res2.get_parent().id}")

    print("\n" + "=" * 80)
    print("CHARGE DISTRIBUTION COMPARISON")
    print("=" * 80)

    # Compare charge distributions
    charges1 = [getattr(a, 'charge', 0.0) for a in atoms1 if hasattr(a, 'charge')]
    charges2 = [getattr(a, 'charge', 0.0) for a in atoms2 if hasattr(a, 'charge')]

    print(f"\nFormat 1 charges:")
    print(f"  Min: {min(charges1):.4f}")
    print(f"  Max: {max(charges1):.4f}")
    print(f"  Mean: {np.mean(charges1):.4f}")
    print(f"  Sum: {sum(charges1):.4f}")

    print(f"\nFormat 2 charges:")
    print(f"  Min: {min(charges2):.4f}")
    print(f"  Max: {max(charges2):.4f}")
    print(f"  Mean: {np.mean(charges2):.4f}")
    print(f"  Sum: {sum(charges2):.4f}")

    print("\n" + "=" * 80)
    print("COORDINATE COMPARISON")
    print("=" * 80)

    # Compare coordinate ranges
    coords1 = np.array([a.coord for a in atoms1])
    coords2 = np.array([a.coord for a in atoms2])

    print(f"\nFormat 1 coordinates:")
    print(f"  X range: [{coords1[:,0].min():.3f}, {coords1[:,0].max():.3f}]")
    print(f"  Y range: [{coords1[:,1].min():.3f}, {coords1[:,1].max():.3f}]")
    print(f"  Z range: [{coords1[:,2].min():.3f}, {coords1[:,2].max():.3f}]")
    print(f"  Center: [{coords1[:,0].mean():.3f}, {coords1[:,1].mean():.3f}, {coords1[:,2].mean():.3f}]")

    print(f"\nFormat 2 coordinates:")
    print(f"  X range: [{coords2[:,0].min():.3f}, {coords2[:,0].max():.3f}]")
    print(f"  Y range: [{coords2[:,1].min():.3f}, {coords2[:,1].max():.3f}]")
    print(f"  Z range: [{coords2[:,2].min():.3f}, {coords2[:,2].max():.3f}]")
    print(f"  Center: [{coords2[:,0].mean():.3f}, {coords2[:,1].mean():.3f}, {coords2[:,2].mean():.3f}]")

    # Calculate RMSD between first N matching atoms
    n_match = min(len(atoms1), len(atoms2))
    rmsd = np.sqrt(np.mean((coords1[:n_match] - coords2[:n_match])**2))
    print(f"\nRMSD between first {n_match} atoms: {rmsd:.6f} Å")

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    if len(atoms1) == len(atoms2):
        print("\n✅ Atom counts match")
    else:
        print(f"\n⚠️  Atom counts differ by {abs(len(atoms1) - len(atoms2))}")

    if rmsd < 0.001:
        print("✅ Coordinates are identical (RMSD < 0.001 Å)")
    elif rmsd < 0.1:
        print(f"✅ Coordinates are very similar (RMSD = {rmsd:.6f} Å)")
    else:
        print(f"⚠️  Coordinates differ significantly (RMSD = {rmsd:.6f} Å)")

    charge_diff = abs(sum(charges1) - sum(charges2))
    if charge_diff < 0.001:
        print("✅ Total charges are identical")
    else:
        print(f"⚠️  Total charges differ by {charge_diff:.4f}")

    if atoms_with_differences == 0:
        print("✅ All field comparisons passed")
    else:
        print(f"⚠️  {atoms_with_differences}/{num_to_compare} atoms have differences")

    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
