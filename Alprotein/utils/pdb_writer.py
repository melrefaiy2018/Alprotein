"""Utilities for writing structures in extended PDB format."""

from __future__ import annotations

from typing import Any

from Bio.PDB.Structure import Structure


class ExtendedPDBWriter:
    """Write a BioPython structure to an extended PDB file with charges."""

    def __init__(self, filename: str) -> None:
        """Create writer for ``filename``."""
        self.filename = filename

    def __str__(self) -> str:
        """Return a readable representation of the writer."""
        return f"ExtendedPDBWriter(filename={self.filename})"

    __repr__ = __str__

    def __eq__(self, other: Any) -> bool:
        """Return ``True`` if ``other`` writes to the same file."""
        if not isinstance(other, ExtendedPDBWriter):
            return False
        return self.filename == other.filename

    def write_pdb_line(self, atom_record: str, **atom_info: Any) -> None:
        """Write a single PDB line to ``filename``."""
        if atom_info["residue_name"] in (
            "CLA",
            "CLB",
            "CHL",
            "BCR",
            "SQD",
            "LHG",
            "LMG",
            "LMU",
            "PQN",
            "DGD",
        ):
            atom_record = "HETATM"

        pdb_line = self._format_atom_record(atom_record)
        pdb_line += self._format_atom_index(atom_info["atom_index"])
        pdb_line += self._format_atom_name(atom_info["atom_name"])
        pdb_line += self._format_residue_name(atom_info["residue_name"])
        pdb_line += self._format_chain_id(atom_info["chain_id"])
        pdb_line += self._format_residue_index(atom_info["residue_index"])
        pdb_line += self._format_xyz(atom_info["x"])
        pdb_line += self._format_xyz(atom_info["y"])
        pdb_line += " "
        pdb_line += self._format_xyz(atom_info["z"])
        pdb_line += self._format_element_name(atom_info["element_name"])
        pdb_line += self._format_charge(atom_info["residue_name"], atom_info["charge"])
        pdb_line += "\n"

        with open(self.filename, "a", encoding="utf-8") as pdb_file:
            pdb_file.write(pdb_line)

    @staticmethod
    def _format_atom_record(atom_record: str) -> str:
        return f"{atom_record:<6s}"

    @staticmethod
    def _format_atom_index(atom_index: int) -> str:
        return f"{atom_index:>5d} "

    @staticmethod
    def _format_atom_name(atom_name: str) -> str:
        return f"{ExtendedPDBWriter._construct_atomtype_string(atom_name)} "

    @staticmethod
    def _format_residue_name(residue_name: str) -> str:
        return f"{residue_name:>3s} "

    @staticmethod
    def _format_chain_id(chain_id: str) -> str:
        return f"{chain_id:>1s}"

    @staticmethod
    def _format_residue_index(residue_index: str) -> str:
        return f"{residue_index:>4s}    "

    @staticmethod
    def _format_xyz(value: float) -> str:
        return f"{float(value):>8.3f}"

    @staticmethod
    def _format_element_name(element_name: str) -> str:
        return " " * 21 + f"{element_name:>2s}  "

    @staticmethod
    def _format_charge(residue_name: str, charge: float) -> str:
        if residue_name in ("CLA", "CLB", "CHL"):
            return "  None"
        return f"{float(charge):7.4f}"

    @staticmethod
    def _construct_atomtype_string(atomtype: str) -> str:
        if len(atomtype) > 4:
            raise ValueError("atom type can't be more than 5 characters")
        if len(atomtype) == 4:
            return f"{atomtype:<4s}"
        if len(atomtype) == 3:
            return f"{atomtype:>4s}"
        return f"{atomtype:^4s}"

    def write_structure(self, structure: Structure) -> None:
        """Write all atoms from ``structure`` to ``filename``."""
        for chain in structure.pdb.get_chains():
            for residue in chain.get_residues():
                for atom in residue.get_atoms():
                    element = "Mg" if atom.get_name() == "MG" else atom.element
                    atom_info = {
                        "atom_record": "ATOM",
                        "atom_index": atom.get_serial_number(),
                        "atom_name": atom.get_name(),
                        "residue_name": residue.resname,
                        "chain_id": chain.id,
                        "residue_index": str(residue.id[1]),
                        "x": float(atom.coord[0]),
                        "y": float(atom.coord[1]),
                        "z": float(atom.coord[2]),
                        "element_name": element,
                        "charge": getattr(atom, "charge", 0),
                    }
                    self.write_pdb_line(**atom_info)

