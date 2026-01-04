"""
Unit tests for PDB parser dual format support.

Tests format detection, parsing logic for both Format 1 (Extended MCCE)
and Format 2 (Standard with variable spacing), and backward compatibility.
"""

import pytest
import numpy as np
from Bio import PDB

from Alprotein.core.protein_structure import ProteinStructure


class TestFormatDetection:
    """Test automatic format detection."""

    def test_detect_extended_mcce_format_with_conformer(self):
        """Verify Format 1 detection with conformer ID."""
        line = 'ATOM      1  N   PRO 40005_000 -13.450 -11.506  21.344   1.500      -0.055      BK____M000'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        assert protein._detect_format(line) == 'extended_mcce'

    def test_detect_extended_mcce_format_with_letter_prefix(self):
        """Verify Format 1 detection with letter-prefixed resseq."""
        line = 'ATOM      1  N   VAL s0019    -15.382 -8.237  22.728                     N   -0.443     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        assert protein._detect_format(line) == 'standard_variable'

    def test_detect_standard_variable_format(self):
        """Verify Format 2 detection with 11 tokens."""
        line = 'ATOM      1    N PRO 40005     -13.450 -11.506  21.344                      N   -0.055     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        assert protein._detect_format(line) == 'standard_variable'

    def test_detect_standard_variable_format_various_tokens(self):
        """Verify Format 2 detection works with different values."""
        line = 'ATOM      5    O PRO 40005     -15.629  -9.576  20.228                      O   -0.402     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        assert protein._detect_format(line) == 'standard_variable'


class TestFormat1Parsing:
    """Test Extended MCCE format (Format 1) parsing."""

    def test_parse_format1_with_conformer(self):
        """Test Format 1 parsing extracts all fields correctly."""
        line = 'ATOM      1  N   PRO 40005_000 -13.450 -11.506  21.344   1.500      -0.055      BK____M000'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_extended_mcce_format(line, 0)

        assert data is not None
        assert data['serial'] == 1
        assert data['name'] == 'N'
        assert data['resname'] == 'PRO'
        assert data['resseq'] == 40005
        assert data['conformer_id'] == '000'
        assert data['charge'] == -0.055
        assert data['identifier'] == 'BK'
        assert data['occupancy'] == 1.500
        assert data['chain_id'] == '4'  # Col 21 contains '4' which is parsed as chain ID
        assert data['element'] == 'N'

    def test_parse_format1_coordinates(self):
        """Verify Format 1 parses coordinates correctly."""
        line = 'ATOM      1  N   PRO 40005_000 -13.450 -11.506  21.344   1.500      -0.055      BK____M000'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_extended_mcce_format(line, 0)

        assert np.allclose(data['coord'], np.array([-13.450, -11.506, 21.344]))

    def test_parse_format1_different_conformer(self):
        """Test Format 1 with different conformer ID."""
        line = 'ATOM      6  CB  PRO 40005_001 -15.684 -11.816  22.011   2.000      -0.022      01O000M000'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_extended_mcce_format(line, 0)

        assert data['conformer_id'] == '001'
        assert data['identifier'] == '01'
        assert data['charge'] == -0.022
        assert data['occupancy'] == 2.000

    def test_parse_format1_element_inference(self):
        """Verify Format 1 infers element from atom name."""
        line = 'ATOM      2  CA  PRO 40005_000 -14.779 -11.759  20.776   2.000      -0.002      BK____M000'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_extended_mcce_format(line, 0)

        assert data['element'] == 'C'  # Inferred from 'CA'


class TestFormat2Parsing:
    """Test Standard variable spacing format (Format 2) parsing."""

    def test_parse_format2_basic(self):
        """Test Format 2 parsing extracts all fields correctly."""
        line = 'ATOM      1    N PRO 40005     -13.450 -11.506  21.344                      N   -0.055     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_standard_variable_format(line, 0)

        assert data is not None
        assert data['serial'] == 1
        assert data['name'] == 'N'
        assert data['resname'] == 'PRO'
        assert data['resseq'] == 40005
        assert data['element'] == 'N'
        assert data['charge'] == -0.055
        assert data['identifier'] == 'BK'
        assert data['conformer_id'] == '000'
        assert data['chain_id'] == 'A'
        assert data['occupancy'] == 1.0   # Default occupancy

    def test_parse_format2_coordinates(self):
        """Verify Format 2 parses coordinates correctly."""
        line = 'ATOM      1    N PRO 40005     -13.450 -11.506  21.344                      N   -0.055     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_standard_variable_format(line, 0)

        assert np.allclose(data['coord'], np.array([-13.450, -11.506, 21.344]))

    def test_parse_format2_negative_coords(self):
        """Test Format 2 handles negative coordinates."""
        line = 'ATOM      5    O PRO 40005     -15.629  -9.576  20.228                      O   -0.402     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_standard_variable_format(line, 0)

        assert data['coord'][0] == -15.629
        assert data['coord'][1] == -9.576
        assert data['coord'][2] == 20.228

    def test_parse_format2_element_inference(self):
        """Verify Format 2 element inference works."""
        line = 'ATOM      2   CA PRO 40005     -14.779 -11.759  20.776                      C   -0.002     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_standard_variable_format(line, 0)

        assert data['element'] == 'C'

    def test_parse_format2_different_identifier(self):
        """Test Format 2 with different identifier."""
        line = 'ATOM      6   CB PRO 40005     -15.684 -11.816  22.011                      C   -0.022     01'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_standard_variable_format(line, 0)

        assert data['identifier'] == '01'
        assert data['charge'] == -0.022

    def test_parse_format2_positive_charge(self):
        """Test Format 2 with positive charge."""
        line = 'ATOM      3   HA PRO 40005     -14.816 -12.655  20.189                      H    0.093     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_standard_variable_format(line, 0)

        assert data['charge'] == 0.093

    def test_parse_format2_wrong_token_count(self):
        """Test Format 2 raises error for wrong token count."""
        # Only 10 tokens instead of 11
        line = 'ATOM      1    N PRO 40005     -13.450 -11.506  21.344                      N   -0.055'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')

        with pytest.raises(ValueError, match="Expected 10 or 11 tokens"):
            protein._parse_standard_variable_format(line, 0)


class TestUnifiedParsingInterface:
    """Test the unified _parse_extended_atom_line interface."""

    def test_routing_to_format1(self):
        """Verify routing to Format 1 parser."""
        line = 'ATOM      1  N   PRO 40005_000 -13.450 -11.506  21.344   1.500      -0.055      BK____M000'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_extended_atom_line(line, 0)

        assert data is not None
        assert data['conformer_id'] == '000'  # Format 1 has conformers

    def test_routing_to_format2(self):
        """Verify routing to Format 2 parser."""
        line = 'ATOM      1    N PRO 40005     -13.450 -11.506  21.344                      N   -0.055     BK'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_extended_atom_line(line, 0)

        assert data is not None
        assert data['conformer_id'] == ''  # Format 2 has no conformers

    def test_short_line_returns_none(self):
        """Verify short lines return None."""
        line = 'ATOM      1  N'
        protein = ProteinStructure(PDB.Structure.Structure('test'), 'test')
        data = protein._parse_extended_atom_line(line, 0)

        assert data is None


class TestBackwardCompatibility:
    """Test backward compatibility with existing files."""

    def test_load_format1_file(self):
        """Test loading existing Format 1 file."""
        pdb_path = 'Alprotein/pdb/CP24_most_occ_pH7.pdb'
        protein = ProteinStructure.from_file_extended(pdb_path, 'test')

        atoms = list(protein.pdb.get_atoms())
        assert len(atoms) > 0, "Should load atoms successfully"

        # Verify charges loaded
        atoms_with_charges = protein.get_atoms_with_charges()
        assert len(atoms_with_charges) > 0, "Should have atoms with charges"

        # Check conformer IDs parsed
        atoms_with_conformers = [a for a in atoms if hasattr(a, 'conformer_id') and a.conformer_id != '']
        assert len(atoms_with_conformers) > 0, "Should have atoms with conformer IDs"

    def test_validation_format1(self):
        """Test validation works for Format 1."""
        pdb_path = 'Alprotein/pdb/CP24_most_occ_pH7.pdb'
        protein = ProteinStructure.from_file_extended(pdb_path, 'test')

        validation = protein.validate_extended_format()
        assert validation['atoms_with_charges'] > 0, "Should have charged atoms"


class TestIntegrationBothFormats:
    """Integration tests for both formats."""

    def test_load_format2_file(self):
        """Test loading Format 2 file."""
        pdb_path = 'Alprotein/pdb/CP24_extended_most_occ_pH8.pdb'
        protein = ProteinStructure.from_file_extended(pdb_path, 'test')

        atoms = list(protein.pdb.get_atoms())
        assert len(atoms) > 0, "Should load atoms successfully"

        # Verify charges loaded
        atoms_with_charges = protein.get_atoms_with_charges()
        assert len(atoms_with_charges) > 0, "Should have atoms with charges"

        # Verify elements assigned
        atoms_with_elements = [a for a in atoms if hasattr(a, 'element') and a.element]
        assert len(atoms_with_elements) == len(atoms), "All atoms should have elements"

    def test_both_formats_produce_compatible_structures(self):
        """Verify both formats produce compatible atom structures."""
        pdb_path_fmt1 = 'Alprotein/pdb/CP24_most_occ_pH7.pdb'
        pdb_path_fmt2 = 'Alprotein/pdb/CP24_extended_most_occ_pH8.pdb'

        protein1 = ProteinStructure.from_file_extended(pdb_path_fmt1, 'test1')
        protein2 = ProteinStructure.from_file_extended(pdb_path_fmt2, 'test2')

        atoms1 = list(protein1.pdb.get_atoms())
        atoms2 = list(protein2.pdb.get_atoms())

        # Both should have atoms
        assert len(atoms1) > 0
        assert len(atoms2) > 0

        # Both should have charged atoms via get_atoms_with_charges()
        charged_atoms1 = protein1.get_atoms_with_charges()
        charged_atoms2 = protein2.get_atoms_with_charges()

        assert len(charged_atoms1) > 0, "Format 1 should have charged atoms"
        assert len(charged_atoms2) > 0, "Format 2 should have charged atoms"

        # Both should have extended_atoms dictionary populated
        assert len(protein1.extended_atoms) > 0, "Format 1 should have extended atoms"
        assert len(protein2.extended_atoms) > 0, "Format 2 should have extended atoms"
