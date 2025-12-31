"""
Enhanced ProteinStructure class for handling extended PDB formats with additional metadata.
"""

import copy
import numpy as np
from Bio import PDB
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from functools import cached_property
from typing import Any, Dict, List, Optional, Set, Tuple
import re
import warnings

from ..utils.pdb_writer import ExtendedPDBWriter
from .enhanced_atom import EnhancedAtom


class ExtendedAtom:
    """
    Enhanced atom class that stores additional metadata from extended PDB formats.
    """
    
    def __init__(self, biopython_atom: Atom, **kwargs):
        """
        Initialize extended atom with BioPython atom and additional metadata.
        
        Args:
            biopython_atom: BioPython Atom object
            **kwargs: Additional metadata (charge, conformer_id, etc.)
        """
        self.atom = biopython_atom
        self.metadata = kwargs
        
        # Store commonly used attributes directly
        self.charge = kwargs.get('charge', 0.0)
        self.conformer_id = kwargs.get('conformer_id', '')
        self.identifier = kwargs.get('identifier', '')
        self.occupancy_extended = kwargs.get('occupancy_extended', 1.0)
        
    def __getattr__(self, name):
        """Delegate to BioPython atom for standard attributes."""
        return getattr(self.atom, name)
    
    def __repr__(self):
        """Return a string representation including atom name and charge."""
        return f"ExtendedAtom({self.atom.name}, charge={self.charge:.3f})"


class ProteinStructure:
    """
    Enhanced protein structure class for handling extended PDB formats.
    
    This class wraps a BioPython Structure object and provides additional
    functionality for charge assignment, conformer handling, and extended metadata.
    """
    
    def __init__(self, pdb_structure: Structure, name: str):
        """
        Initialize a ProteinStructure.
        
        Args:
            pdb_structure: BioPython Structure object
            name: Name identifier for the structure
            
        Raises:
            TypeError: If pdb_structure is not a BioPython Structure
        """
        if not isinstance(pdb_structure, Structure):
            raise TypeError("pdb_structure must be a Biopython Structure object.")
        
        self.pdb = pdb_structure
        self.name = str(name)
        self.extended_atoms = {}  # Store extended atom information
        self.conformer_info = {}  # Store conformer-specific information
        self.parsing_info = {}    # Store information about the parsing process
        
    @classmethod
    def from_file(cls, pdb_path: str, name: str) -> 'ProteinStructure':
        """
        Create a ProteinStructure from a PDB file using standard BioPython parser.
        
        Args:
            pdb_path: Path to the PDB file
            name: Name identifier for the structure
            
        Returns:
            ProteinStructure instance
        """
        try:
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure(name, pdb_path)
            protein_structure = cls(structure, name)

            # Fix element types for non-standard residues
            corrected_count = protein_structure.fix_element_types()
            if corrected_count > 0:
                print(f"Corrected element types for {corrected_count} atoms")

            # Check a few lines of the file for charge information. If it looks
            # like charges are present but the standard parser ignored them,
            # warn the user that they may want to use ``from_file_extended``.
            try:
                with open(pdb_path, "r", encoding="utf-8") as fh:
                    for _ in range(20):
                        line = fh.readline()
                        if not line:
                            break
                        if line.startswith(("ATOM", "HETATM")) and len(line) > 66:
                            extra = line[66:].strip()
                            if extra and re.search(r"[+-]?\d+\.\d+", extra):
                                warnings.warn(
                                    "Charge fields detected in PDB file but not parsed. "
                                    "Use from_file_extended to load charges.",
                                    UserWarning,
                                )
                                break
            except OSError:
                pass

            return protein_structure
        except FileNotFoundError:
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        except Exception as e:
            raise Exception(f"Failed to parse PDB file {pdb_path}: {str(e)}")
    
    @classmethod
    def from_file_extended(cls, pdb_path: str, name: str) -> 'ProteinStructure': #TODO: this is now loading the mcce file without conversion. so we need to update the docstring and the method name.
        """
        Create a ProteinStructure from an extended PDB file with enhanced parsing.
        
        This method creates a custom BioPython structure while preserving all
        extended information from the PDB file.
        
        Args:
            pdb_path: Path to the extended PDB file
            name: Name identifier for the structure
            
        Returns:
            ProteinStructure instance with extended metadata
        """
        try:
            # First create using standard parser to get basic structure
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure(name, pdb_path)
            
            # Create our enhanced structure
            protein_structure = cls(structure, name)
            
            # Parse extended information
            protein_structure._parse_extended_pdb_file(pdb_path)
            
            return protein_structure
            
        except FileNotFoundError:
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        except Exception as e:
            raise Exception(f"Failed to parse extended PDB file {pdb_path}: {str(e)}")

    @classmethod
    def from_file_with_parameters(cls, pdb_path, parameter_config=None, name=None):
        """Load PDB file with optional immediate parameter enrichment.
        
        Args:
            pdb_path (str): Path to PDB file
            parameter_config (dict): Mapping of pigment types to parameter names
            name (str): Optional name for the structure
            
        Returns:
            ProteinStructure: Enhanced structure with optional parameters
        """
        # First load using existing method
        instance = cls.from_file_extended(pdb_path, name)
        
        # Convert to enhanced atoms
        instance._convert_to_enhanced_atoms()
        
        # Ensure all protein atoms have q00 = charge after loading
        instance._assign_q00_to_protein_atoms()
        
        # Optionally enrich with parameters
        if parameter_config:
            instance.enrich_with_parameters(parameter_config)
        
        return instance

    def _convert_to_enhanced_atoms(self):
        """Convert all Biopython atoms to EnhancedAtom objects."""
        for model in self.pdb:
            for chain in model:
                for residue in chain:
                    # Create enhanced atoms list
                    enhanced_atoms = []
                    for atom in residue:
                        enhanced_atom = EnhancedAtom(atom)
                        enhanced_atoms.append(enhanced_atom)
                    
                    # Replace the residue's atom list
                    residue.child_list = enhanced_atoms
                    
                    # Update the residue's internal dictionary
                    residue.child_dict = {atom.id: atom for atom in enhanced_atoms}

    def enrich_with_parameters(self, parameter_config):
        """Enrich pigment atoms with calculation parameters.
        
        Args:
            parameter_config (dict): Mapping like {'ChlorophyllA': 'chl_a_default'}
        """
        from ..data.parameters import get_site_energy_parameters
        
        pigment_residues = self._identify_pigment_residues()
        
        for residue in pigment_residues:
            pigment_type = self._get_pigment_type(residue.resname)
            if pigment_type and pigment_type in parameter_config:
                param_name = parameter_config[pigment_type]
                params = get_site_energy_parameters(param_name)
                self._attach_parameters_to_residue(residue, params, pigment_type)

    def _identify_pigment_residues(self):
        """Identify all pigment residues in the structure."""
        pigment_resnames = ['CLA', 'CHL', 'BCL', 'CHD', 'PEO', 'CAR']  # Add more as needed
        pigment_residues = []
        
        for model in self.pdb:
            for chain in model:
                for residue in chain:
                    if residue.resname in pigment_resnames:
                        pigment_residues.append(residue)
        
        return pigment_residues

    def _get_pigment_type(self, resname):
        """Map residue name to pigment type."""
        mapping = {
            'CLA': 'ChlorophyllA',
            'CHL': 'ChlorophyllA', 
            'BCL': 'ChlorophyllB',
            'CHD': 'ChlorophyllB',
            # Add more mappings as needed
        }
        return mapping.get(resname)

    def _attach_parameters_to_residue(self, residue, params, pigment_type):
        """Attach calculation parameters to all atoms in a residue."""
        # FIXED: Use the correct field names from SITE_ENERGY_DATA
        atom_names = params.get('atom', [])  # Note: field is 'atom', not 'atom_names'
        q00_charges = params.get('q_00', [])
        q11_charges = params.get('q_11', [])
        q01_charges = params.get('q_01', [])
        
        # Create lookup dictionary for O(1) access
        charge_lookup = {}
        for i, name in enumerate(atom_names):
            charge_lookup[name] = {
                'q00': q00_charges[i] if i < len(q00_charges) else 0.0,
                'q11': q11_charges[i] if i < len(q11_charges) else 0.0,
                'q01': q01_charges[i] if i < len(q01_charges) else 0.0,
            }
        
        # Attach parameters to matching atoms
        for atom in residue:
            if atom.name in charge_lookup:
                charges = charge_lookup[atom.name]
                atom.q00 = charges['q00']
                atom.q11 = charges['q11']
                atom.q01 = charges['q01']
                atom.is_pigment_atom = True
                atom.pigment_type = pigment_type
                atom.calculation_ready = True
    
    @classmethod
    def from_extended_lines(cls, pdb_lines: List[str], name: str) -> 'ProteinStructure':
        """
        Create a ProteinStructure from extended PDB format lines.
        
        This method creates a complete BioPython structure from scratch,
        parsing the extended format properly.
        
        Args:
            pdb_lines: List of PDB format lines
            name: Name identifier for the structure
            
        Returns:
            ProteinStructure instance
        """
        # Create empty structure
        structure = Structure(name)
        current_model = Model(0)
        structure.add(current_model)
        
        protein_structure = cls(structure, name)
        chains = {}
        
        # Parse each line
        for line_num, line in enumerate(pdb_lines):
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    atom_data = protein_structure._parse_extended_atom_line(line, line_num)
                    if atom_data:
                        # Get or create chain
                        chain_id = atom_data['chain_id']
                        if chain_id not in chains:
                            chains[chain_id] = Chain(chain_id)
                            current_model.add(chains[chain_id])
                        
                        chain = chains[chain_id]
                        
                        # Get or create residue
                        res_id = (atom_data['hetfield'], atom_data['resseq'], atom_data['icode'])
                        try:
                            residue = chain[res_id]
                        except KeyError:
                            residue = Residue(res_id, atom_data['resname'], atom_data['segid'])
                            chain.add(residue)
                        
                        # Create atom
                        atom = Atom(
                            name=atom_data['name'],
                            coord=atom_data['coord'],
                            bfactor=atom_data['bfactor'],
                            occupancy=atom_data['occupancy'],
                            altloc=atom_data['altloc'],
                            fullname=atom_data['fullname'],
                            serial_number=atom_data['serial'],
                            element=atom_data['element']
                        )
                        
                        # Add extended attributes
                        atom.charge = atom_data['charge']
                        atom.conformer_id = atom_data.get('conformer_id', '')
                        atom.identifier = atom_data.get('identifier', '')
                        atom.occupancy_extended = atom_data.get('occupancy_extended', 1.0)
                        
                        # Add to residue
                        residue.add(atom)
                        
                        # Store extended information
                        atom_key = f"{chain_id}_{res_id[1]}_{atom_data['name']}"
                        protein_structure.extended_atoms[atom_key] = {
                            'charge': atom_data['charge'],
                            'conformer_id': atom_data.get('conformer_id', ''),
                            'identifier': atom_data.get('identifier', ''),
                            'original_line': line.strip()
                        }
                        
                except Exception as e:
                    print(f"Warning: Could not parse line {line_num + 1}: {e}")
                    continue
        
        # Fix element types for all atoms after parsing
        corrected_elements = protein_structure.fix_element_types()
        
        protein_structure.parsing_info = {
            'total_lines': len(pdb_lines),
            'atoms_parsed': len(protein_structure.extended_atoms),
            'chains_found': list(chains.keys()),
            'element_types_corrected': corrected_elements
        }
        
        return protein_structure
    
    def _fix_element_type(self, atom_name: str, resname: str, current_element: str = None) -> str:
        """
        Fix element type for atoms, especially in non-standard residues like CLA, CHL.
        
        Args:
            atom_name: Atom name from PDB (e.g., 'CHA', 'N1A', 'MG')
            resname: Residue name (e.g., 'CLA', 'CHL', 'ASP')
            current_element: Current element assignment (might be 'X' for unknown)
            
        Returns:
            Corrected element symbol
        """
        # If we already have a valid element (not 'X' or empty), keep it
        if current_element and current_element.strip() not in ['X', '']:
            return current_element.strip()
        
        # Clean atom name
        atom_name = atom_name.strip()
        
        # Special cases for common non-standard residues
        if resname in ['CLA', 'CHL']:  # Chlorophyll
            # Common chlorophyll atoms
            if atom_name == 'MG':
                return 'Mg'
            elif atom_name.startswith(('C', 'CHA', 'CHB', 'CHC', 'CHD')):
                return 'C'
            elif atom_name.startswith('N'):
                return 'N'
            elif atom_name.startswith('O'):
                return 'O'
        
        # For other residues or general case, use first letter heuristic
        # Remove digits and common modifiers
        first_char = atom_name[0].upper()
        
        # Handle special cases where first letter might not be element
        if first_char == 'H':
            return 'H'  # Hydrogen
        elif first_char == 'C':
            return 'C'  # Carbon
        elif first_char == 'N':
            return 'N'  # Nitrogen
        elif first_char == 'O':
            return 'O'  # Oxygen
        elif first_char == 'P':
            return 'P'  # Phosphorus
        elif first_char == 'S':
            return 'S'  # Sulfur
        elif first_char == 'F':
            return 'F'  # Fluorine
        elif first_char == 'I':
            return 'I'  # Iodine
        elif first_char == 'B':
            return 'B'  # Boron
        elif first_char == 'K':
            return 'K'  # Potassium
        elif first_char == 'M':
            # Could be Mg, Mn, Mo, etc. Check second character
            if len(atom_name) > 1:
                second_char = atom_name[1].upper()
                if second_char == 'G':
                    return 'Mg'  # Magnesium
                elif second_char == 'N':
                    return 'Mn'  # Manganese
                elif second_char == 'O':
                    return 'Mo'  # Molybdenum
            return 'M'  # Fallback to first letter
        elif first_char == 'Z':
            # Zinc
            if len(atom_name) > 1 and atom_name[1].upper() == 'N':
                return 'Zn'
            return 'Z'
        elif first_char == 'F' and len(atom_name) > 1 and atom_name[1].upper() == 'E':
            return 'Fe'  # Iron
        elif first_char == 'C' and len(atom_name) > 1:
            second_char = atom_name[1].upper()
            if second_char == 'A':
                return 'Ca'  # Calcium
            elif second_char == 'L':
                return 'Cl'  # Chlorine
            elif second_char == 'U':
                return 'Cu'  # Copper
            return 'C'  # Default to Carbon
        else:
            # For other cases, just use the first character
            return first_char
    
    def _parse_extended_atom_line(self, line: str, line_num: int) -> Dict[str, Any]:
            """
            Parse an extended PDB atom line with enhanced format.
            
            Expected format:
            ATOM      1  N   ASN A0002_000 -41.421  26.792   2.542   1.500      -0.434      BK____M000
            Positions:                                           54-60   66-78      84+
                                                            occupancy charge   identifier
            
            Args:
                line: PDB line to parse
                line_num: Line number for error reporting
                
            Returns:
                Dictionary with parsed atom data or None if parsing fails
            """
            if len(line) < 54:  # Minimum length for coordinates
                return None
            
            try:
                # Standard PDB fields
                record_type = line[0:6].strip()
                serial = int(line[6:11].strip())
                name = line[12:16].strip()
                altloc = line[16:17].strip()
                resname = line[17:20].strip()
                chain_id = line[21:22].strip() or 'A'  # Default to 'A' if empty
                
                # Parse coordinates
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coord = np.array([x, y, z], dtype=np.float64)
                
                # Parse occupancy and B-factor
                occupancy = float(line[54:60].strip()) if len(line) > 60 and line[54:60].strip() else 1.0
                bfactor = float(line[60:66].strip()) if len(line) > 66 and line[60:66].strip() else 0.0
                
                # Parse element - check if provided in standard PDB element field (positions 76-78)
                element = ''
                if len(line) > 78:
                    element = line[76:78].strip()
                
                # Fix element type using our correction method
                element = self._fix_element_type(name, resname, element)
                
                # Handle extended residue ID format (e.g., A0002_000)
                resseq_field = line[22:30].strip()  # Extended field for residue sequence
                if '_' in resseq_field:
                    # Extended format: A0002_000
                    resseq_str, conformer_id = resseq_field.split('_', 1)
                    resseq = int(re.sub(r'[A-Za-z]', '', resseq_str))  # Extract numeric part
                else:
                    # Standard format
                    resseq = int(line[22:26].strip())
                    conformer_id = ''
                
                # Parse icode
                icode = line[26:27].strip() if len(line) > 27 else ' '
                
                # Parse charge and identifier based on your specific format
                charge = 0.0
                identifier = ''
                
                # For your format, charge appears to be around positions 66-78
                # and identifier starts around position 84
                if len(line) > 66:
                    # Extract the charge field - look for it in the region after B-factor
                    # Your example: "      -0.434      BK____M000"
                    charge_and_id_region = line[66:]
                    
                    # Split into potential fields
                    parts = charge_and_id_region.split()
                    
                    for i, part in enumerate(parts):
                        part = part.strip()
                        if not part:
                            continue
                        
                        # Try to identify charge (should be a float, typically with . and +/-)
                        if ('.' in part or part.startswith(('+', '-'))) and not part.startswith(('BK', 'CLA', 'CHL')):
                            try:
                                potential_charge = float(part)
                                # Reasonable charge range for atoms
                                if -5.0 <= potential_charge <= 5.0:
                                    charge = potential_charge
                                    continue
                            except ValueError:
                                pass
                        
                        # If it's not a charge, it might be the identifier
                        # Identifiers typically contain letters
                        if any(c.isalpha() for c in part) and part not in ['1.00', '0.00']:
                            identifier = part
                
                # Fallback: if no identifier found in parts, check end of line
                if not identifier and len(line) > 80:
                    # Look for identifier at the end
                    end_part = line[80:].strip()
                    if end_part and any(c.isalpha() for c in end_part):
                        identifier = end_part
                
                return {
                    'name': name,
                    'fullname': f" {name:<3}",  # Formatted for BioPython
                    'resname': resname,
                    'chain_id': chain_id,
                    'resseq': resseq,
                    'icode': icode,
                    'coord': coord,
                    'occupancy': occupancy,
                    'bfactor': bfactor,
                    'altloc': altloc,
                    'element': element,
                    'serial': serial,
                    'hetfield': ' ' if record_type == 'ATOM' else 'H',
                    'segid': '',
                    'conformer_id': conformer_id,
                    'charge': charge,
                    'identifier': identifier,
                    'occupancy_extended': occupancy,
                    'original_line': line.strip()
                }
                
            except Exception as e:
                print(f"Warning: Could not parse line {line_num + 1}: {e}")
                return None
            

    def fix_element_types(self) -> int:
        """
        Fix element types for all atoms in the structure, particularly for non-standard residues.
        
        Returns:
            Number of atoms with corrected element types
        """
        corrected_count = 0
        
        for atom in self.pdb.get_atoms():
            residue = atom.get_parent()
            resname = residue.get_resname().strip()
            atom_name = atom.name.strip()
            current_element = getattr(atom, 'element', '').strip()
            
            # Fix element type
            corrected_element = self._fix_element_type(atom_name, resname, current_element)
            
            # Update if different
            if corrected_element != current_element:
                atom.element = corrected_element
                corrected_count += 1
        
        return corrected_count
    
    def _parse_extended_pdb_file(self, pdb_path: str) -> int:
        """
        Parse extended information from PDB file and associate with existing atoms.
        
        Args:
            pdb_path: Path to the PDB file
            
        Returns:
            Number of atoms with extended information assigned
        """
        assigned_count = 0
        
        try:
            with open(pdb_path, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            return 0
        
        # Create atom mapping for quick lookup
        atom_map = {}
        for atom in self.pdb.get_atoms():
            residue = atom.get_parent()
            chain_id = residue.get_parent().id
            res_id = residue.id[1]
            atom_name = atom.name.strip()
            key = (chain_id, res_id, atom_name)
            atom_map[key] = atom
        
        # Parse extended information
        for line_num, line in enumerate(lines):
            if line.startswith(('ATOM', 'HETATM')):
                atom_data = self._parse_extended_atom_line(line, line_num)
                if atom_data:
                    # Find matching atom
                    key = (atom_data['chain_id'], atom_data['resseq'], atom_data['name'])
                    
                    if key in atom_map:
                        atom = atom_map[key]
                        
                        # Assign extended attributes
                        atom.charge = atom_data['charge']
                        atom.q00 = atom_data['charge']  # For protein atoms, q00 = charge
                        atom.conformer_id = atom_data.get('conformer_id', '')
                        atom.identifier = atom_data.get('identifier', '')
                        atom.occupancy_extended = atom_data.get('occupancy_extended', 1.0)
                        
                        # Fix element type if needed
                        if 'element' in atom_data and atom_data['element']:
                            corrected_element = self._fix_element_type(
                                atom_data['name'], 
                                atom_data['resname'], 
                                atom_data['element']
                            )
                            atom.element = corrected_element
                        
                        # Store in extended atoms dictionary
                        atom_key = f"{atom_data['chain_id']}_{atom_data['resseq']}_{atom_data['name']}"
                        self.extended_atoms[atom_key] = {
                            'charge': atom_data['charge'],
                            'conformer_id': atom_data.get('conformer_id', ''),
                            'identifier': atom_data.get('identifier', ''),
                            'original_line': line.strip()
                        }
                        
                        assigned_count += 1
                    else:
                        # Try fuzzy matching for different chain naming
                        for (stored_chain, stored_res, stored_atom), atom_obj in atom_map.items():
                            if stored_res == atom_data['resseq'] and stored_atom == atom_data['name']:
                                atom_obj.charge = atom_data['charge']
                                atom_obj.q00 = atom_data['charge']  # For protein atoms, q00 = charge
                                atom_obj.conformer_id = atom_data.get('conformer_id', '')
                                atom_obj.identifier = atom_data.get('identifier', '')
                                
                                # Fix element type
                                if 'element' in atom_data and atom_data['element']:
                                    corrected_element = self._fix_element_type(
                                        atom_data['name'], 
                                        atom_data['resname'], 
                                        atom_data['element']
                                    )
                                    atom_obj.element = corrected_element
                                
                                assigned_count += 1
                                break
        
        # After parsing, ensure all protein atoms (non-pigment) have q00 = charge
        self._assign_q00_to_protein_atoms()
        
        # Always fix element types for all atoms after parsing
        corrected_elements = self.fix_element_types()
        
        self.parsing_info = {
            'extended_atoms_found': assigned_count,
            'total_lines_processed': len(lines),
            'element_types_corrected': corrected_elements
        }
        
        return assigned_count
    
    def _assign_q00_to_protein_atoms(self):
        """
        Ensure all protein atoms (non-pigment atoms) have q00 = charge.
        This is called after parsing to make sure the assignment is complete.
        """
        pigment_resnames = {'CLA', 'CHL', 'BCL', 'CHD', 'PEO', 'CAR'}
        
        assigned_count = 0
        for atom in self.pdb.get_atoms():
            residue = atom.get_parent()
            resname = residue.get_resname().strip()
            
            # Only assign q00 to non-pigment atoms
            if resname not in pigment_resnames:
                if hasattr(atom, 'charge') and atom.charge is not None:
                    atom.q00 = atom.charge
                    assigned_count += 1
                elif not hasattr(atom, 'q00') or atom.q00 is None:
                    # If no charge found, set q00 to 0.0 for protein atoms
                    atom.q00 = 0.0
        
        if assigned_count > 0:
            print(f"✓ Assigned q00 charges to {assigned_count} protein atoms")
    
    def set_atomic_charges(self, charge_file_path: str = None) -> None:
        """
        Add atomic charges from a file or extract from extended PDB format.
        
        Args:
            charge_file_path: Optional path to separate charge file
        """
        if charge_file_path is None:
            # Try to use already parsed extended information
            charges_assigned = len([atom for atom in self.pdb.get_atoms() 
                                 if hasattr(atom, 'charge') and atom.charge != 0.0])
            if charges_assigned > 0:
                print(f"Using {charges_assigned} charges from extended PDB format")
                return
            else:
                raise ValueError("No charges found in structure and no charge file provided")
        
        # Load from separate file (existing functionality)
        charges_assigned = 0
        try:
            with open(charge_file_path, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"Charge file not found: {charge_file_path}")
        
        for atom in self.pdb.get_atoms():
            for line in lines:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and \
                   int(line[22:26].strip()) == atom.get_parent().id[1] and \
                   line[12:16].strip() == atom.name:
                    try:
                        charge = float(line[54:62].strip())
                        atom.charge = charge
                        charges_assigned += 1
                    except (ValueError, IndexError):
                        pass
                    break
        
        print(f"Assigned charges to {charges_assigned} atoms from {charge_file_path}")
    
    def get_atoms_with_charges(self) -> List:
        """Get all atoms that have assigned charges."""
        return [atom for atom in self.pdb.get_atoms() 
                if hasattr(atom, 'charge') and atom.charge != 0.0]
    
    def get_conformer_info(self) -> Dict[str, List]:
        """
        Get information about different conformers in the structure.
        
        Returns:
            Dictionary mapping conformer IDs to lists of atoms
        """
        conformers = {}
        for atom in self.pdb.get_atoms():
            if hasattr(atom, 'conformer_id') and atom.conformer_id:
                if atom.conformer_id not in conformers:
                    conformers[atom.conformer_id] = []
                conformers[atom.conformer_id].append(atom)
        return conformers
    
    def get_atoms_by_conformer(self, conformer_id: str) -> List:
        """Get all atoms belonging to a specific conformer."""
        return [atom for atom in self.pdb.get_atoms() 
                if hasattr(atom, 'conformer_id') and atom.conformer_id == conformer_id]
    
    def validate_extended_format(self) -> Dict[str, Any]:
        """
        Validate the extended PDB format parsing.
        
        Returns:
            Dictionary with validation results
        """
        total_atoms = len(list(self.pdb.get_atoms()))
        atoms_with_charges = len(self.get_atoms_with_charges())
        atoms_with_conformers = len([atom for atom in self.pdb.get_atoms() 
                                   if hasattr(atom, 'conformer_id') and atom.conformer_id])
        
        validation = {
            'total_atoms': total_atoms,
            'atoms_with_charges': atoms_with_charges,
            'atoms_with_conformers': atoms_with_conformers,
            'charge_coverage': atoms_with_charges / total_atoms if total_atoms > 0 else 0,
            'conformer_coverage': atoms_with_conformers / total_atoms if total_atoms > 0 else 0,
            'extended_atoms_stored': len(self.extended_atoms),
            'conformer_ids': list(self.get_conformer_info().keys()),
            'parsing_info': self.parsing_info
        }
        
        return validation
    
    def get_residues_by_name(self, resname: str, conformer_id: str = None) -> List:
        """
        Get residues by name, optionally filtered by conformer.
        
        Args:
            resname: Residue name to search for
            conformer_id: Optional conformer ID filter
            
        Returns:
            List of matching residues
        """
        residues = [res for res in self.pdb.get_residues() 
                   if res.get_resname().strip() == resname]
        
        if conformer_id is not None:
            # Filter by conformer
            filtered_residues = []
            for res in residues:
                atoms_in_conformer = [atom for atom in res.get_atoms() 
                                    if hasattr(atom, 'conformer_id') and 
                                    atom.conformer_id == conformer_id]
                if atoms_in_conformer:
                    filtered_residues.append(res)
            residues = filtered_residues
        
        return residues
    
    def export_extended_pdb(self, output_path: str) -> None:
        """
        Export structure back to extended PDB format.
        
        Args:
            output_path: Output file path
        """
        open(output_path, 'w', encoding='utf-8').close()
        writer = ExtendedPDBWriter(output_path)
        writer.write_structure(self.pdb)
        with open(output_path, 'a', encoding='utf-8') as f:
            f.write('END\n')
    
    @cached_property
    def mass(self) -> float:
        """Calculate total mass of the structure."""
        return sum(getattr(atom, 'mass', 0.0) for atom in self.pdb.get_atoms())
    
    @cached_property
    def center_of_mass(self) -> np.ndarray:
        """Calculate center of mass of the structure."""
        com = np.zeros(3, dtype=np.float64)
        total_mass = 0.0
        
        for atom in self.pdb.get_atoms():
            atom_mass = getattr(atom, 'mass', 0.0)
            if atom_mass > 0:
                com += atom.coord * atom_mass
                total_mass += atom_mass
        
        return com / total_mass if total_mass > 0 else com
    
    def get_background_atoms(self, exclude_atoms: Optional[Set] = None) -> List:
        """
        Get all protein atoms that are not in the exclude set and have charges.
        
        Args:
            exclude_atoms: Set of atoms to exclude (e.g., pigment atoms)
            
        Returns:
            List of background atoms with charges
        """
        background_atoms = []
        exclude_atoms = exclude_atoms or set()
        for atom in self.pdb.get_atoms():
            if atom not in exclude_atoms and hasattr(atom, 'charge') and atom.charge != 0.0:
                background_atoms.append(atom)
        return background_atoms
    
    def __repr__(self) -> str:
        """Return a summary representation of the structure."""
        validation = self.validate_extended_format()
        return (
            f"<ProteinStructure: {self.name} "
            f"({validation['total_atoms']} atoms, "
            f"{validation['atoms_with_charges']} with charges, "
            f"{len(validation['conformer_ids'])} conformers)>"
        )
