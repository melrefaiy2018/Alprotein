#!/usr/bin/env python3
"""
Inspect a single residue from a PDB file to see all atom fields.

This script loads a PDB file and prints detailed information about
all atoms in a specified residue, showing all parsed attributes.
"""

from Alprotein import ProteinStructure
import numpy as np


def inspect_residue(pdb_file, residue_index=0):
    """
    Load PDB and inspect all fields for one residue.

    Args:
        pdb_file: Path to PDB file
        residue_index: Index of residue to inspect (default: 0 = first residue)
    """
    print("=" * 80)
    print(f"RESIDUE INSPECTION: {pdb_file}")
    print("=" * 80)

    # Load the structure
    print(f"\nLoading {pdb_file}...")
    protein = ProteinStructure.from_file_extended(pdb_file, 'test')

    # Get all residues
    residues = list(protein.pdb.get_residues())
    print(f"Total residues in structure: {len(residues)}")

    if residue_index >= len(residues):
        print(f"ERROR: Residue index {residue_index} out of range (max: {len(residues)-1})")
        return

    # Get the target residue
    residue = residues[residue_index]
    atoms = list(residue.get_atoms())

    print("\n" + "=" * 80)
    print(f"RESIDUE {residue_index}: {residue.resname} {residue.id}")
    print("=" * 80)

    # Residue info
    print(f"\nResidue Name: {residue.resname}")
    print(f"Residue ID: {residue.id}")
    print(f"  - hetfield: {residue.id[0]!r}")
    print(f"  - sequence number: {residue.id[1]}")
    print(f"  - insertion code: {residue.id[2]!r}")

    # Parent chain info
    chain = residue.get_parent()
    print(f"\nChain ID: {chain.id!r}")

    print(f"\nNumber of atoms: {len(atoms)}")

    # Detailed atom information
    print("\n" + "=" * 80)
    print("DETAILED ATOM INFORMATION")
    print("=" * 80)

    for i, atom in enumerate(atoms):
        print(f"\n--- Atom {i}: {atom.name} ---")

        # Basic BioPython attributes
        print("\nStandard BioPython Attributes:")
        print(f"  name: {atom.name!r}")
        print(f"  element: {atom.element!r}")
        print(f"  serial_number: {atom.serial_number}")
        print(f"  fullname: {atom.fullname!r}")
        print(f"  altloc: {atom.altloc!r}")

        # Coordinates
        print(f"\nCoordinates:")
        print(f"  coord: {atom.coord}")
        print(f"    x: {atom.coord[0]:.6f}")
        print(f"    y: {atom.coord[1]:.6f}")
        print(f"    z: {atom.coord[2]:.6f}")

        # Occupancy and B-factor
        print(f"\nOccupancy and B-factor:")
        print(f"  occupancy: {atom.occupancy}")
        print(f"  bfactor: {atom.bfactor}")

        # Extended attributes (from custom parser)
        print(f"\nExtended Attributes (from custom parser):")

        # Charge
        if hasattr(atom, 'charge'):
            print(f"  charge: {atom.charge}")
        else:
            print(f"  charge: NOT SET")

        # Identifier
        if hasattr(atom, 'identifier'):
            print(f"  identifier: {atom.identifier!r}")
        else:
            print(f"  identifier: NOT SET")

        # Conformer ID
        if hasattr(atom, 'conformer_id'):
            print(f"  conformer_id: {atom.conformer_id!r}")
        else:
            print(f"  conformer_id: NOT SET")

        # Occupancy extended
        if hasattr(atom, 'occupancy_extended'):
            print(f"  occupancy_extended: {atom.occupancy_extended}")
        else:
            print(f"  occupancy_extended: NOT SET")

        # Calculation attributes
        print(f"\nCalculation Attributes:")

        if hasattr(atom, 'q00'):
            print(f"  q00 (ground state charge): {atom.q00}")
        else:
            print(f"  q00: NOT SET")

        if hasattr(atom, 'q11'):
            print(f"  q11 (excited state charge): {atom.q11}")
        else:
            print(f"  q11: NOT SET")

        if hasattr(atom, 'q01'):
            print(f"  q01 (transition charge): {atom.q01}")
        else:
            print(f"  q01: NOT SET")

        # All attributes (catch anything we missed)
        print(f"\nAll atom attributes:")
        atom_attrs = [attr for attr in dir(atom) if not attr.startswith('_')]
        for attr in sorted(atom_attrs):
            if attr not in ['name', 'element', 'serial_number', 'fullname', 'altloc',
                           'coord', 'occupancy', 'bfactor', 'charge', 'identifier',
                           'conformer_id', 'occupancy_extended', 'q00', 'q11', 'q01']:
                try:
                    value = getattr(atom, attr)
                    if not callable(value):
                        print(f"  {attr}: {value!r}")
                except:
                    pass

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Check which atoms have charges
    atoms_with_charge = [a for a in atoms if hasattr(a, 'charge') and a.charge != 0]
    print(f"\nAtoms with charge attribute: {len([a for a in atoms if hasattr(a, 'charge')])}/{len(atoms)}")
    print(f"Atoms with non-zero charge: {len(atoms_with_charge)}/{len(atoms)}")

    if atoms_with_charge:
        total_charge = sum(a.charge for a in atoms_with_charge)
        print(f"Total charge for this residue: {total_charge:.6f}")

        # Charge distribution
        charges = [a.charge for a in atoms_with_charge]
        print(f"Charge range: [{min(charges):.4f}, {max(charges):.4f}]")

    # Check for conformer IDs
    conformer_ids = set(getattr(a, 'conformer_id', '') for a in atoms)
    if conformer_ids - {''}:
        print(f"\nConformer IDs present: {sorted(conformer_ids - {''})}")
    else:
        print(f"\nNo conformer IDs (Format 2 style)")

    # Check identifiers
    identifiers = set(getattr(a, 'identifier', '') for a in atoms)
    print(f"\nIdentifiers present: {sorted(identifiers - {''})}")

    print("\n" + "=" * 80)


if __name__ == "__main__":
    # Inspect first residue from the extended file
    pdb_file = 'examples/data/extended_most_occ_pH7_IsiA.pdb'

    print("\nInspecting FIRST residue (index 0):")
    inspect_residue(pdb_file, residue_index=0)

    # Optionally inspect another residue
    # print("\n\n")
    # print("Inspecting SECOND residue (index 1):")
    # inspect_residue(pdb_file, residue_index=1)
