"""
Abstract pigment classes and pigment implementations.

This module provides a minimal abstract interface for pigment calculations
while maintaining BioPython familiarity and backward compatibility.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Optional, Any, Union
import numpy as np
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from ..data.parameters import get_coupling_parameters, get_site_energy_parameters


class AbstractPigment(ABC):
    """
    Abstract base class defining the minimal interface required for pigment calculations.
    
    This class enforces consistency across different pigment types while allowing
    flexibility in implementation details. Only calculation-critical methods are abstract.
    """
    
    def __init__(self, residue_object: Residue, **kwargs):
        """
        Initialize the pigment with a BioPython residue.
        
        Args:
            residue_object: BioPython residue object
            **kwargs: Pigment-specific parameters
        """
        self.residue = residue_object
        self.name = self._create_name()
        self.params = {}
        self._load_calculation_params(**kwargs)
    
    def _create_name(self) -> str:
        """Create a unique name for the pigment - concrete implementation."""
        chain_id = self.residue.get_parent().id
        res_name = self.residue.resname.strip()
        res_id = self.residue.id[1]
        return f'{chain_id}_{res_name}_{res_id}'
    
    # ==========================================
    # ABSTRACT METHODS - Must be implemented by subclasses
    # ==========================================
    
    @abstractmethod
    def get_site_energy_atoms(self) -> List[Tuple[str, float, float]]:
        """
        Get atoms and charges for site energy calculations.
        
        Returns:
            List of (atom_name, q_00_charge, q_11_charge) tuples
        """
        pass
    
    @abstractmethod
    def get_coupling_atoms(self) -> List[Tuple[str, float]]:
        """
        Get atoms and transition charges for coupling calculations.
        
        Returns:
            List of (atom_name, transition_charge) tuples
        """
        pass
    
    @abstractmethod
    def get_transition_dipole_vector(self, transition: str = 'qy') -> np.ndarray:
        """
        Get the transition dipole moment vector for a specific transition.
        
        Args:
            transition: Type of transition ('qy', 'qx', etc.)
            
        Returns:
            Normalized 3D vector
        """
        pass
    
    @abstractmethod
    def get_vacuum_energy(self) -> float:
        """
        Get the vacuum transition energy for this pigment type.
        
        Returns:
            Vacuum energy in cm⁻¹
        """
        pass
    
    @abstractmethod
    def validate_for_calculations(self) -> Dict[str, bool]:
        """
        Validate that this pigment has all required atoms and parameters for calculations.
        
        Returns:
            Dictionary with validation results for different calculation types
        """
        pass
    
    @abstractmethod
    def _load_calculation_params(self, **kwargs) -> None:
        """
        Load pigment-specific calculation parameters.
        
        Args:
            **kwargs: Parameter specifications (tresp_dict_name, cdc_dict_name, etc.)
        """
        pass
    
    # ==========================================
    # CONCRETE METHODS - Shared implementations
    # ==========================================
    
    def get_atom_coord(self, atom_name: str) -> np.ndarray:
        """
        Get coordinates of a specific atom by name.
        BioPython-style method that all pigments can use.
        
        Args:
            atom_name: Name of the atom
            
        Returns:
            3D coordinates as numpy array
            
        Raises:
            ValueError: If atom is not found in the residue
        """
        try:
            return self.residue[atom_name.strip()].coord
        except KeyError:
            available_atoms = [atom.name for atom in self.residue.get_atoms()]
            raise ValueError(
                f"Atom '{atom_name.strip()}' not found in pigment {self.name}. "
                f"Available atoms: {available_atoms}"
            )
    
    def get_atoms(self) -> List[Atom]:
        """
        Get all atoms in this pigment (BioPython-style method).
        
        Returns:
            List of BioPython Atom objects
        """
        return list(self.residue.get_atoms())
    
    def get_id(self) -> str:
        """Get the unique identifier for this pigment (BioPython-style)."""
        return self.name
    
    def get_resname(self) -> str:
        """Get the residue name (BioPython-style)."""
        return self.residue.resname.strip()
    
    def get_coord(self) -> np.ndarray:
        """
        Get the center coordinates of the pigment (BioPython-style).
        Default implementation - subclasses can override for pigment-specific centers.
        """
        # Default: geometric center of all atoms
        coords = np.array([atom.coord for atom in self.residue.get_atoms()])
        return np.mean(coords, axis=0)
    
    def has_atom(self, atom_name: str) -> bool:
        """Check if atom exists in this pigment."""
        try:
            self.get_atom_coord(atom_name)
            return True
        except ValueError:
            return False
    
    def get_atoms_with_charges(self, charge_type: str = 'q_00') -> List[Atom]:
        """
        Get all atoms in this pigment that have the specified charge type.
        
        Args:
            charge_type: Type of charge to check for ('q_00', 'q_11', 'q_01', 'charge')
            
        Returns:
            List of atoms with the specified charge type
        """
        atoms_with_charges = []
        for atom in self.residue.get_atoms():
            if hasattr(atom, charge_type):
                atoms_with_charges.append(atom)
        return atoms_with_charges
    
    def get_atoms_by_charge_type(self, charge_type: str = 'q_00') -> List[Tuple[Atom, float]]:
        """
        Get atoms that have a specific charge type assigned.
        
        Args:
            charge_type: 'q_00', 'q_11', 'q_01', or 'charge'
            
        Returns:
            List of (atom, charge_value) tuples
        """
        atoms_charges = []
        for atom in self.residue.get_atoms():
            if hasattr(atom, charge_type):
                charge_value = getattr(atom, charge_type)
                atoms_charges.append((atom, charge_value))
        return atoms_charges
    
    def set_charge_state(self, state: str = 'ground') -> None:
        """
        Set the active charge state for all atoms in this pigment.
        
        Args:
            state: 'ground' (q_00), 'excited' (q_11), or 'transition' (q_01)
        """
        charge_mapping = {
            'ground': 'q_00',
            'excited': 'q_11', 
            'transition': 'q_01'
        }
        
        if state not in charge_mapping:
            raise ValueError(f"Invalid state '{state}'. Must be one of: {list(charge_mapping.keys())}")
        
        charge_attr = charge_mapping[state]
        atoms_updated = 0
        
        for atom in self.residue.get_atoms():
            if hasattr(atom, charge_attr):
                atom.charge = getattr(atom, charge_attr)
                atoms_updated += 1
        
        print(f"Set {atoms_updated} atoms to {state} state ({charge_attr}) for {self.name}")
    
    def get_charge_summary(self) -> Dict[str, int]:
        """
        Get a summary of charge assignments for this pigment.
        
        Returns:
            Dictionary with counts of atoms having each charge type
        """
        summary = {
            'q_00': 0,
            'q_11': 0, 
            'q_01': 0,
            'charge': 0
        }
        
        for atom in self.residue.get_atoms():
            for charge_type in summary.keys():
                if hasattr(atom, charge_type):
                    summary[charge_type] += 1
        
        return summary
    
    @property
    def xyz(self) -> np.ndarray:
        """
        Get the center position of the pigment (backward compatibility).
        
        Returns:
            3D coordinates as numpy array
        """
        return self.get_coord()
    
    def get_dipole_dir(self, tdm: str = 'qy') -> np.ndarray:
        """
        Get the transition dipole moment direction (backward compatibility).
        
        Args:
            tdm: Transition dipole type ('qy' or 'qx')
            
        Returns:
            Normalized 3D direction vector
        """
        return self.get_transition_dipole_vector(tdm)
    
    def __repr__(self) -> str:
        """Return a readable identifier for the pigment."""
        return f"{self.__class__.__name__}({self.name})"
    
    def __iter__(self):
        """Allow iteration over atoms (BioPython-style)."""
        return iter(self.residue.get_atoms())


class ChlorophyllA(AbstractPigment):
    """
    Chlorophyll A implementation of the abstract pigment interface.
    """
    def __init__(self, residue, params_attached=False, cdc_dict_name=None):
        """Initialize pigment with optional pre-attached parameters.
        
        Args:
            residue: Biopython residue object
            params_attached (bool): If True, parameters are already attached to atoms
            cdc_dict_name (str): Parameter set name (used only if params_attached=False)
        """
        # Call parent constructor but skip parameter loading
        self.residue = residue
        self.name = self._create_name()
        self.params = {}
        self.atoms = list(residue)
        
        print(f"DEBUG: Initializing ChlorophyllA {self.name}")
        print(f"  params_attached: {params_attached}")
        print(f"  cdc_dict_name: {cdc_dict_name}")
        
        # CRITICAL FIX: Only load parameters if they are NOT already attached
        if params_attached:
            print(f"  Using pre-attached parameters for {self.name}")
            self._extract_attached_parameters()
            # Mark that we have enhanced atom access
            self._has_enhanced_atoms = True
        else:
            print(f"  Loading parameters from dictionary for {self.name}")
            # Default to 'CLA' if no parameter name provided
            if cdc_dict_name is None:
                cdc_dict_name = 'CLA'  # Use valid default
                print(f"  Using default CDC parameter set: {cdc_dict_name}")
            else:
                print(f"  Using provided CDC parameter set: {cdc_dict_name}")
            self._load_calculation_params(cdc_dict_name=cdc_dict_name)
            self._has_enhanced_atoms = False
        
        # Validate what we loaded
        site_atoms = self.get_site_energy_atoms()
        coupling_atoms = self.get_coupling_atoms()
        print(f"  Result: {len(site_atoms)} site energy atoms, {len(coupling_atoms)} coupling atoms")

    def _extract_attached_parameters(self):
        """Extract parameters from enhanced atoms for backward compatibility."""
        atom_names = []
        q00_charges = []
        q11_charges = []
        q01_charges = []
        
        for atom in self.atoms:
            if hasattr(atom, 'calculation_ready') and atom.calculation_ready:
                atom_names.append(atom.name)
                q00_charges.append(atom.q00)
                q11_charges.append(atom.q11)
                if hasattr(atom, 'q01') and atom.q01 is not None:
                    q01_charges.append(atom.q01)
        
        # Store in params dict for backward compatibility
        self.params['atom_names'] = atom_names
        self.params['q_00'] = q00_charges
        self.params['q_11'] = q11_charges
        if q01_charges:
            self.params['q_01'] = q01_charges

    def get_site_energy_atoms_direct(self):
        """Get atoms with charges using direct access (optimized version).
        
        Returns:
            list: Tuples of (atom_name, q00, q11) for atoms with attached parameters
        """
        return [
            (atom.name, atom.q00, atom.q11)
            for atom in self.atoms 
            if hasattr(atom, 'calculation_ready') and atom.calculation_ready
        ]

    def get_transition_charges_direct(self):
        """Get transition charges using direct access (optimized version).
        
        Returns:
            list: Tuples of (atom_name, q01) for atoms with transition charges
        """
        return [
            (atom.name, atom.q01)
            for atom in self.atoms
            if (hasattr(atom, 'calculation_ready') and 
                atom.calculation_ready and 
                hasattr(atom, 'q01') and 
                atom.q01 is not None)
        ]
    
    def get_site_energy_atoms(self) -> List[Tuple[str, float, float]]:
        """Get atoms and charges for CDC site energy calculations."""
        # Use optimized direct access if enhanced atoms are available
        if hasattr(self, '_has_enhanced_atoms') and self._has_enhanced_atoms:
            return self.get_site_energy_atoms_direct()
        
        # Fallback to legacy parameter dictionary method
        if not self.params.get('cdc_loaded'):
            return []
        
        atoms_charges = []
        atom_names = self.params['q_atoms']
        q_00_charges = self.params['q_00']
        q_11_charges = self.params['q_11']
        
        for i, atom_name in enumerate(atom_names):
            if i < len(q_00_charges) and i < len(q_11_charges):
                if self.has_atom(atom_name):
                    atoms_charges.append((
                        atom_name,
                        q_00_charges[i],
                        q_11_charges[i]
                    ))
        
        return atoms_charges
    
    def get_coupling_atoms(self) -> List[Tuple[str, float]]:
        """Get atoms and transition charges for TrEsp coupling calculations."""
        # Use optimized direct access if enhanced atoms are available
        if hasattr(self, '_has_enhanced_atoms') and self._has_enhanced_atoms:
            return self.get_transition_charges_direct()
        
        # Fallback to legacy parameter dictionary method
        if not self.params.get('tresp_loaded'):
            return []
        
        atoms_charges = []
        atom_names = self.params['tresp_atoms']
        transition_charges = self.params['tresp_pc']
        
        for i, atom_name in enumerate(atom_names):
            if i < len(transition_charges):
                if self.has_atom(atom_name):
                    atoms_charges.append((atom_name, transition_charges[i]))
        
        return atoms_charges
    
    def get_transition_dipole_vector(self, transition: str = 'qy') -> np.ndarray:
        """
        Get transition dipole direction for chlorophyll A.
        
        Uses N1B-N1D vector for Qy, N1A-N1C vector for Qx.
        """
        # Define atom patterns for both naming conventions
        atom_patterns = {
            'qy': [('N1B', 'N1D'), ('NB', 'ND')],
            'qx': [('N1A', 'N1C'), ('NA', 'NC')]
        }
        
        # Try each pattern for the given transition type
        for atoms in atom_patterns.get(transition, []):
            try:
                coord1 = self.get_atom_coord(atoms[0])
                coord2 = self.get_atom_coord(atoms[1])
                
                v = coord1 - coord2
                norm = np.linalg.norm(v)
                return v / norm if norm > 0 else np.zeros(3)
            except ValueError:
                continue
        
        print(f"Warning: Could not calculate {transition} dipole direction for {self.name}")
        return np.zeros(3)
    
    def get_vacuum_energy(self) -> float:
        """Vacuum energy for chlorophyll A."""
        return 14900.0  # cm⁻¹
    
    def get_coord(self) -> np.ndarray:
        """
        Get center position - for chlorophylls, use Mg position if available.
        """
        try:
            return self.get_atom_coord('MG')
        except ValueError:
            # Fallback to nitrogen center
            try:
                nitrogen_patterns = [
                    ['N1A', 'N1B', 'N1C', 'N1D'],
                    ['NA', 'NB', 'NC', 'ND']
                ]
                
                for pattern in nitrogen_patterns:
                    try:
                        coords = [self.get_atom_coord(n) for n in pattern]
                        return np.mean(coords, axis=0)
                    except ValueError:
                        continue
            except ValueError:
                pass
            
            # Final fallback to geometric center
            return super().get_coord()
    
    def validate_for_calculations(self) -> Dict[str, bool]:
        """
        Validate chlorophyll A for different calculation types.
        """
        validation = {
            'site_energy_ready': False,
            'coupling_ready': False,
            'dipole_ready': False,
            'structure_complete': False
        }
        
        # Check site energy readiness
        site_atoms = self.get_site_energy_atoms()
        validation['site_energy_ready'] = len(site_atoms) > 0
        
        # Check coupling readiness
        coupling_atoms = self.get_coupling_atoms()
        validation['coupling_ready'] = len(coupling_atoms) > 0
        
        # Check dipole calculation readiness
        qy_dipole = self.get_transition_dipole_vector('qy')
        qx_dipole = self.get_transition_dipole_vector('qx')
        validation['dipole_ready'] = (np.linalg.norm(qy_dipole) > 0 and 
                                    np.linalg.norm(qx_dipole) > 0)
        
        # Check key structural atoms
        key_atoms = ['MG', 'N1A', 'N1B', 'N1C', 'N1D', 'CHA', 'CHB', 'CHC', 'CHD']
        missing_key_atoms = [atom for atom in key_atoms if not self.has_atom(atom)]
        validation['structure_complete'] = len(missing_key_atoms) == 0
        
        return validation
    
    def _load_calculation_params(self, tresp_dict_name: Optional[str] = None, 
                               cdc_dict_name: Optional[str] = None, 
                               **kwargs) -> None:
        """
        Load calculation parameters for chlorophyll A.
        """
        print(f"    DEBUG: _load_calculation_params called for {self.name}")
        print(f"      tresp_dict_name: {tresp_dict_name}")
        print(f"      cdc_dict_name: {cdc_dict_name}")
        print(f"      kwargs: {kwargs}")
        
        # Log any unexpected parameters for debugging
        expected_params = {'tresp_dict_name', 'cdc_dict_name'}
        unexpected_params = set(kwargs.keys()) - expected_params
        if unexpected_params:
            print(f"      Note: Ignoring unexpected parameters for {self.name}: {unexpected_params}")
        
        # Load TrEsp parameters
        if tresp_dict_name:
            try:
                tresp_params = get_coupling_parameters(tresp_dict_name)
                self.params.update(tresp_params)
                self.params['tresp_loaded'] = True
                print(f"      Loaded TrEsp parameters '{tresp_dict_name}' for {self.name}")
            except KeyError:
                print(f"      Warning: TrEsp parameters '{tresp_dict_name}' not found")
                self.params['tresp_loaded'] = False
        else:
            print(f"      No TrEsp parameters requested")
            self.params['tresp_loaded'] = False
        
        # Load CDC parameters
        if cdc_dict_name:
            try:
                print(f"      Attempting to load CDC parameters '{cdc_dict_name}'")
                cdc_params = get_site_energy_parameters(cdc_dict_name)
                print(f"      Successfully loaded CDC parameters: {list(cdc_params.keys())}")
                self.params.update(cdc_params)
                self.params['q_atoms'] = [name.strip() for name in self.params['atom']]
                self.params['cdc_loaded'] = True
                print(f"      Loaded CDC parameters '{cdc_dict_name}' for {self.name}: {len(self.params['q_atoms'])} atoms")
            except KeyError as e:
                print(f"      Error: CDC parameters '{cdc_dict_name}' not found: {e}")
                self.params['cdc_loaded'] = False
            except Exception as e:
                print(f"      Error loading CDC parameters: {e}")
                self.params['cdc_loaded'] = False
        else:
            print(f"      No CDC parameters requested")
            self.params['cdc_loaded'] = False
        
        # Assign charges to atoms
        print(f"      Assigning charges to atoms...")
        self._assign_pigment_charges()
    
    def _assign_pigment_charges(self) -> None:
        """
        Assign charges to pigment atoms from parameter files.
        This is the same implementation as the original Chlorophyll class.
        """
        charges_assigned = {'q_00': 0, 'q_11': 0, 'q_01': 0}
        missing_atoms = {'cdc': [], 'tresp': []}
        
        # Assign q_00 and q_11 charges from CDC parameters
        if self.params.get('cdc_loaded'):
            atom_names = self.params['q_atoms']
            q_00_charges = self.params['q_00']
            q_11_charges = self.params['q_11']
            
            # Validate charge array lengths
            if len(atom_names) != len(q_00_charges) or len(atom_names) != len(q_11_charges):
                print(f"Warning: Charge array length mismatch for {self.name}")
                print(f"  Atoms: {len(atom_names)}, q_00: {len(q_00_charges)}, q_11: {len(q_11_charges)}")
            
            for i, atom_name in enumerate(atom_names):
                try:
                    atom = self.residue[atom_name.strip()]
                    
                    # Assign ground state charge (q_00)
                    if i < len(q_00_charges):
                        atom.q_00 = q_00_charges[i]
                        charges_assigned['q_00'] += 1
                    
                    # Assign excited state charge (q_11) 
                    if i < len(q_11_charges):
                        atom.q_11 = q_11_charges[i]
                        charges_assigned['q_11'] += 1
                    
                    # Set the default charge attribute to q_00 for compatibility
                    if i < len(q_00_charges):
                        atom.charge = q_00_charges[i]
                    
                except KeyError:
                    missing_atoms['cdc'].append(atom_name.strip())
                    continue
        
        # Assign q_01 charges from TrEsp parameters
        if self.params.get('tresp_loaded'):
            tresp_atom_names = self.params['tresp_atoms']
            tresp_charges = self.params['tresp_pc']
            
            # Validate charge array length
            if len(tresp_atom_names) != len(tresp_charges):
                print(f"Warning: TrEsp charge array length mismatch for {self.name}")
                print(f"  Atoms: {len(tresp_atom_names)}, charges: {len(tresp_charges)}")
            
            for i, atom_name in enumerate(tresp_atom_names):
                try:
                    atom = self.residue[atom_name.strip()]
                    
                    # Assign transition charge (q_01)
                    if i < len(tresp_charges):
                        atom.q_01 = tresp_charges[i]
                        charges_assigned['q_01'] += 1
                    
                except KeyError:
                    missing_atoms['tresp'].append(atom_name.strip())
                    continue
        
        # Report results
        print(f"Assigned charges to {self.name}: q_00={charges_assigned['q_00']}, q_11={charges_assigned['q_11']}, q_01={charges_assigned['q_01']}")
        
        # Report missing atoms if any
        if missing_atoms['cdc']:
            print(f"  Missing CDC atoms: {missing_atoms['cdc']}")
        if missing_atoms['tresp']:
            print(f"  Missing TrEsp atoms: {missing_atoms['tresp']}")
        
        # Validate total charges (optional)
        self._validate_charge_totals()
    
    def _validate_charge_totals(self) -> None:
        """
        Validate that the total charges are reasonable.
        Reports the sum of each charge type for debugging.
        """
        charge_sums = {'q_00': 0.0, 'q_11': 0.0, 'q_01': 0.0}
        
        for atom in self.residue.get_atoms():
            for charge_type in charge_sums.keys():
                if hasattr(atom, charge_type):
                    charge_sums[charge_type] += getattr(atom, charge_type)
        
        # Report charge sums (helpful for debugging)
        print(f"  Charge sums for {self.name}: " + 
              ", ".join([f"{ct}={charge_sums[ct]:.3f}" for ct in charge_sums.keys()]))
        
        # Check if charges are significantly unbalanced (optional warning)
        for charge_type, total in charge_sums.items():
            if abs(total) > 2.0:  # Adjust threshold as needed
                print(f"  Warning: Large total charge for {charge_type}: {total:.3f}")


class ChlorophyllB(ChlorophyllA):
    """
    Chlorophyll B - inherits most behavior from ChlorophyllA but with different vacuum energy.
    """
    
    def __init__(self, residue, params_attached=False, cdc_dict_name=None):
        """Initialize ChlorophyllB with proper parameter defaults."""
        # Default to CHL parameters for ChlorophyllB
        if cdc_dict_name is None:
            cdc_dict_name = 'CHL'
        super().__init__(residue, params_attached, cdc_dict_name)
    
    def get_vacuum_energy(self) -> float:
        """Vacuum energy for chlorophyll B."""
        return 15674.0  # cm⁻¹


class Pheophytin(AbstractPigment):
    """
    Pheophytin implementation - similar to chlorophyll but without Mg center.
    """
    
    def get_site_energy_atoms(self) -> List[Tuple[str, float, float]]:
        """Get atoms and charges for CDC site energy calculations."""
        # Implementation similar to ChlorophyllA but without Mg
        if not self.params.get('cdc_loaded'):
            return []
        
        atoms_charges = []
        atom_names = self.params['q_atoms']
        q_00_charges = self.params['q_00']
        q_11_charges = self.params['q_11']
        
        for i, atom_name in enumerate(atom_names):
            if i < len(q_00_charges) and i < len(q_11_charges):
                if self.has_atom(atom_name):
                    atoms_charges.append((
                        atom_name,
                        q_00_charges[i],
                        q_11_charges[i]
                    ))
        
        return atoms_charges
    
    def get_coupling_atoms(self) -> List[Tuple[str, float]]:
        """Get atoms and transition charges for TrEsp coupling calculations."""
        if not self.params.get('tresp_loaded'):
            return []
        
        atoms_charges = []
        atom_names = self.params['tresp_atoms']
        transition_charges = self.params['tresp_pc']
        
        for i, atom_name in enumerate(atom_names):
            if i < len(transition_charges):
                if self.has_atom(atom_name):
                    atoms_charges.append((atom_name, transition_charges[i]))
        
        return atoms_charges
    
    def get_transition_dipole_vector(self, transition: str = 'qy') -> np.ndarray:
        """Same dipole calculation as chlorophyll."""
        return ChlorophyllA.get_transition_dipole_vector(self, transition)
    
    def get_vacuum_energy(self) -> float:
        """Vacuum energy for pheophytin."""
        return 14700.0  # cm⁻¹ (typically slightly lower than chlorophyll)
    
    def get_coord(self) -> np.ndarray:
        """
        Get center position - for pheophytins, use nitrogen center (no Mg).
        """
        try:
            nitrogen_patterns = [
                ['N1A', 'N1B', 'N1C', 'N1D'],
                ['NA', 'NB', 'NC', 'ND']
            ]
            
            for pattern in nitrogen_patterns:
                try:
                    coords = [self.get_atom_coord(n) for n in pattern]
                    return np.mean(coords, axis=0)
                except ValueError:
                    continue
        except ValueError:
            pass
        
        # Fallback to geometric center
        return super().get_coord()
    
    def validate_for_calculations(self) -> Dict[str, bool]:
        """Validate pheophytin - similar to chlorophyll but no Mg requirement."""
        validation = {
            'site_energy_ready': False,
            'coupling_ready': False,
            'dipole_ready': False,
            'structure_complete': False
        }
        
        # Check site energy readiness
        site_atoms = self.get_site_energy_atoms()
        validation['site_energy_ready'] = len(site_atoms) > 0
        
        # Check coupling readiness
        coupling_atoms = self.get_coupling_atoms()
        validation['coupling_ready'] = len(coupling_atoms) > 0
        
        # Check dipole calculation readiness
        qy_dipole = self.get_transition_dipole_vector('qy')
        qx_dipole = self.get_transition_dipole_vector('qx')
        validation['dipole_ready'] = (np.linalg.norm(qy_dipole) > 0 and 
                                    np.linalg.norm(qx_dipole) > 0)
        
        # Check key structural atoms (no Mg for pheophytin)
        key_atoms = ['N1A', 'N1B', 'N1C', 'N1D', 'CHA', 'CHB', 'CHC', 'CHD']
        missing_key_atoms = [atom for atom in key_atoms if not self.has_atom(atom)]
        validation['structure_complete'] = len(missing_key_atoms) == 0
        
        return validation
    
    def _load_calculation_params(self, tresp_dict_name: Optional[str] = None, 
                               cdc_dict_name: Optional[str] = None, 
                               **kwargs) -> None:
        """Load calculation parameters for pheophytin."""
        # Same implementation as ChlorophyllA
        ChlorophyllA._load_calculation_params(self, tresp_dict_name, cdc_dict_name, **kwargs)


# Factory function for creating pigments
def create_pigment(residue: Residue, pigment_type: Optional[str] = None, **kwargs) -> AbstractPigment:
    """
    Factory function to create appropriate pigment objects based on residue type.
    
    Args:
        residue: BioPython residue object
        pigment_type: Override pigment type (optional)
        **kwargs: Parameters for pigment creation (including params_attached)
        
    Returns:
        Appropriate pigment instance
    """
    # Determine pigment type from residue name if not specified
    if pigment_type is None:
        resname = residue.resname.strip().upper()
        
        # Map residue names to pigment classes
        pigment_mapping = {
            'CLA': ChlorophyllA,
            'CHL': ChlorophyllB,
        }
        
        pigment_class = pigment_mapping.get(resname)
        if pigment_class is None:
            raise ValueError(f"Unknown pigment type for residue {resname}. "
                           f"Available types: {list(pigment_mapping.keys())}")
    else:
        # Direct mapping from string to class
        type_mapping = {
            'chlorophyll_a': ChlorophyllA,
            'chlorophyll_b': ChlorophyllB,
        }
        
        pigment_class = type_mapping.get(pigment_type.lower())
        if pigment_class is None:
            raise ValueError(f"Unknown pigment type '{pigment_type}'. "
                           f"Available types: {list(type_mapping.keys())}")
    
    # FIXED: Pass all kwargs including params_attached to the pigment constructor
    return pigment_class(residue, **kwargs)


# Backward compatibility - alias for existing code
Chlorophyll = ChlorophyllA
