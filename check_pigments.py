#!/usr/bin/env python3
"""Check distribution of pigment vs protein atoms."""

from Alprotein import ProteinStructure

print("=" * 80)
print("Checking pigment vs protein atoms")
print("=" * 80)

# Load both files
print("\nLoading extended file...")
protein_ext = ProteinStructure.from_file_extended('examples/data/extended_most_occ_pH7_IsiA.pdb', 'ext')

print("\nLoading most_occ file...")
protein_most = ProteinStructure.from_file_extended('examples/data/most_occ_pH7_IsiA.pdb', 'most')

# Pigment resnames
pigment_resnames = {'CLA', 'CHL', 'BCL', 'CHD', 'PEO', 'CAR'}

def analyze_atoms(protein, name):
    print(f"\n{name}:")
    print("=" * 80)

    atoms = list(protein.pdb.get_atoms())
    print(f"Total atoms: {len(atoms)}")

    # Categorize by residue type
    pigment_atoms = []
    protein_atoms = []
    other_atoms = []

    for atom in atoms:
        residue = atom.get_parent()
        resname = residue.resname.strip()

        if resname in pigment_resnames:
            pigment_atoms.append(atom)
        elif resname in ['HOH', 'WAT']:  # Water
            other_atoms.append(atom)
        else:
            protein_atoms.append(atom)

    print(f"\nAtom distribution:")
    print(f"  Protein atoms: {len(protein_atoms)}")
    print(f"  Pigment atoms: {len(pigment_atoms)}")
    print(f"  Other atoms (water, etc.): {len(other_atoms)}")

    # Check charges
    protein_with_charge = [a for a in protein_atoms if hasattr(a, 'charge') and a.charge is not None]
    protein_nonzero_charge = [a for a in protein_atoms if hasattr(a, 'charge') and a.charge != 0.0]

    pigment_with_charge = [a for a in pigment_atoms if hasattr(a, 'charge') and a.charge is not None]
    pigment_nonzero_charge = [a for a in pigment_atoms if hasattr(a, 'charge') and a.charge != 0.0]

    print(f"\nProtein atom charges:")
    print(f"  With charge attribute: {len(protein_with_charge)}/{len(protein_atoms)} ({100*len(protein_with_charge)/max(1,len(protein_atoms)):.1f}%)")
    print(f"  With non-zero charge: {len(protein_nonzero_charge)}/{len(protein_atoms)} ({100*len(protein_nonzero_charge)/max(1,len(protein_atoms)):.1f}%)")
    print(f"  With zero charge: {len(protein_atoms) - len(protein_nonzero_charge)}")

    print(f"\nPigment atom charges:")
    print(f"  With charge attribute: {len(pigment_with_charge)}/{len(pigment_atoms)} ({100*len(pigment_with_charge)/max(1,len(pigment_atoms)):.1f}%)")
    print(f"  With non-zero charge: {len(pigment_nonzero_charge)}/{len(pigment_atoms)} ({100*len(pigment_nonzero_charge)/max(1,len(pigment_atoms)):.1f}%)")

    # Check pigment residues
    pigment_residues = set()
    for atom in pigment_atoms:
        residue = atom.get_parent()
        pigment_residues.add((residue.resname, residue.id[1]))

    print(f"\nPigment residues: {len(pigment_residues)}")
    for resname, resnum in sorted(pigment_residues)[:10]:
        print(f"  {resname} {resnum}")

# Analyze both
analyze_atoms(protein_ext, "Extended file")
analyze_atoms(protein_most, "Most_occ file")

# Compare specific residues
print("\n" + "=" * 80)
print("Checking specific atom with zero charge")
print("=" * 80)

# Check HIS 31 CG (which we know has zero charge)
atoms_ext = list(protein_ext.pdb.get_atoms())
atoms_most = list(protein_most.pdb.get_atoms())

def find_atom(atoms, resnum, atomname):
    for atom in atoms:
        residue = atom.get_parent()
        if residue.id[1] == resnum and atom.name == atomname:
            return atom
    return None

cg_ext = find_atom(atoms_ext, 31, 'CG')
cg_most = find_atom(atoms_most, 31, 'CG')

print("\nHIS 31 CG atom:")
if cg_ext:
    print(f"  Extended file: charge={getattr(cg_ext, 'charge', 'N/A')}, q00={getattr(cg_ext, 'q00', 'N/A')}")
if cg_most:
    print(f"  Most_occ file: charge={getattr(cg_most, 'charge', 'N/A')}, q00={getattr(cg_most, 'q00', 'N/A')}")
