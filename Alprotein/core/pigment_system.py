"""
Enhanced pigment system that manages collections of AbstractPigment objects.

This class provides BioPython-style methods for working with pigment collections
while leveraging the abstract pigment interface for calculations.
"""

from typing import Dict, Set, List, Optional, Any, Type, Union, Tuple
import numpy as np
from .protein_structure import ProteinStructure
from .abstract_pigments import AbstractPigment, ChlorophyllA, ChlorophyllB, Pheophytin, create_pigment


class PigmentSystem:
    """
    Enhanced pigment system that manages collections of AbstractPigment objects.
    
    This class provides BioPython-style methods for working with pigment collections
    while leveraging the abstract pigment interface for calculations.
    """
    
    def __init__(self, protein_structure):
        """Initialize pigment system from protein structure.
    
        Args:
            protein_structure: ProteinStructure object with optional enhanced atoms
        """
        self.structure = protein_structure
        self.enhanced_atoms_available = self._check_enhanced_atoms()
        self._pigment_types = {}  # Track pigment types
        self.pigments, self.all_pigment_atoms = self._create_pigments()
    
    def _check_enhanced_atoms(self):
        """Check if structure contains enhanced atoms with parameters."""
        for model in self.structure.pdb:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if hasattr(atom, 'calculation_ready'):
                            return True
        return False

    def _create_pigments(self):
        """Create pigment objects, optimized for enhanced atoms."""
        pigments = {}
        all_pigment_atoms = set()
        for residue in self.structure._identify_pigment_residues():
            pigment_type = self.structure._get_pigment_type(residue.resname)
            if pigment_type:
                # FIXED: Always provide CDC parameters for pigment creation
                # Map residue names to parameter keys in SITE_ENERGY_DATA
                resname = residue.resname.strip().upper()
                cdc_dict_name = 'CLA' if resname in ['CLA'] else 'CHL' if resname in ['CHL'] else None
                
                pigment = create_pigment(
                    residue, 
                    pigment_type=pigment_type.lower().replace('chlorophyll', 'chlorophyll_'),
                    params_attached=self.enhanced_atoms_available,
                    cdc_dict_name=cdc_dict_name  # Provide CDC parameters
                )
                pigment_id = pigment.get_id()
                pigments[pigment_id] = pigment
                self._pigment_types[pigment_id] = type(pigment)
                for atom in pigment.get_atoms():
                    all_pigment_atoms.add(atom)
        return pigments, all_pigment_atoms
    
    def add_pigments_by_residue(self, 
                              resname: str, 
                              pigment_class: Optional[Type[AbstractPigment]] = None,
                              **kwargs) -> int:
        """
        Add pigments by residue name with automatic type detection.
        
        Args:
            resname: Residue name to search for (e.g., 'CLA', 'CHL', 'PHE')
            pigment_class: Optional pigment class override
            **kwargs: Parameters passed to pigment constructor
            
        Returns:
            Number of pigments added
        """
        added_count = 0
        residues = self.structure.get_residues_by_name(resname)
        
        for res in residues:
            try:
                # Use factory function if no specific class provided
                if pigment_class is None:
                    pigment = create_pigment(res, **kwargs)
                else:
                    pigment = pigment_class(res, **kwargs)
                
                if pigment.get_id() not in self.pigments:
                    self.pigments[pigment.get_id()] = pigment
                    self._pigment_types[pigment.get_id()] = type(pigment)
                    
                    # Add atoms to the exclusion set
                    for atom in pigment.get_atoms():
                        self.all_pigment_atoms.add(atom)
                    
                    added_count += 1
                else:
                    print(f"Warning: Pigment {pigment.get_id()} already exists")
                    
            except Exception as e:
                print(f"Warning: Failed to create pigment from residue {res}: {e}")
        
        print(f"Added {added_count} pigments of type {resname}")
        return added_count
    
    def add_pigment(self, pigment: AbstractPigment) -> bool:
        """
        Add a single pigment to the system.
        
        Args:
            pigment: AbstractPigment instance
            
        Returns:
            True if added successfully, False if already exists
        """
        if not isinstance(pigment, AbstractPigment):
            raise TypeError("Pigment must be an instance of AbstractPigment")
        
        pigment_id = pigment.get_id()
        if pigment_id in self.pigments:
            print(f"Warning: Pigment {pigment_id} already exists")
            return False
        
        self.pigments[pigment_id] = pigment
        self._pigment_types[pigment_id] = type(pigment)
        
        # Add atoms to the exclusion set
        for atom in pigment.get_atoms():
            self.all_pigment_atoms.add(atom)
        
        return True
    
    # ==========================================
    # BioPython-style methods
    # ==========================================
    
    def get_pigments(self) -> List[AbstractPigment]:
        """Get all pigments (BioPython-style method)."""
        return list(self.pigments.values())
    
    def get_residues(self) -> List:
        """Get all pigment residues (BioPython-style method)."""
        return [pigment.residue for pigment in self.pigments.values()]
    
    def get_atoms(self) -> List:
        """Get all pigment atoms (BioPython-style method)."""
        atoms = []
        for pigment in self.pigments.values():
            atoms.extend(pigment.get_atoms())
        return atoms
    
    def get_chains(self) -> Set[str]:
        """Get all chain IDs containing pigments."""
        chains = set()
        for pigment in self.pigments.values():
            chain_id = pigment.residue.get_parent().id
            chains.add(chain_id)
        return chains
    
    # ==========================================
    # Pigment-specific methods
    # ==========================================
    
    def get_pigment_names(self) -> List[str]:
        """Get sorted list of pigment names."""
        return sorted(list(self.pigments.keys()))
    
    def get_pigment_by_id(self, pigment_id: str) -> Optional[AbstractPigment]:
        """Get a pigment by its ID."""
        return self.pigments.get(pigment_id)
    
    def get_pigments_by_type(self, pigment_type: Type[AbstractPigment]) -> List[AbstractPigment]:
        """Get all pigments of a specific type."""
        return [pig for pig in self.pigments.values() if isinstance(pig, pigment_type)]
    
    def get_pigments_by_resname(self, resname: str) -> List[AbstractPigment]:
        """Get all pigments with a specific residue name."""
        return [pig for pig in self.pigments.values() if pig.get_resname() == resname]
    
    def get_pigment_positions(self) -> Dict[str, np.ndarray]:
        """
        Get positions of all pigments.
        
        Returns:
            Dictionary mapping pigment names to their positions
        """
        return {name: pig.get_coord() for name, pig in self.pigments.items()}
    
    def get_pigment_distances(self) -> Dict[str, Dict[str, float]]:
        """
        Calculate distances between all pigment pairs.
        
        Returns:
            Nested dictionary with distances between pigments
        """
        positions = self.get_pigment_positions()
        distances = {}
        
        for name1, pos1 in positions.items():
            distances[name1] = {}
            for name2, pos2 in positions.items():
                if name1 != name2:
                    distance = np.linalg.norm(pos1 - pos2)
                    distances[name1][name2] = distance
        
        return distances
    
    # ==========================================
    # Calculation-oriented methods
    # ==========================================
    
    def get_site_energy_data(self) -> Dict[str, List[Tuple[str, float, float]]]:
        """
        Get site energy data for all pigments.
        
        Returns:
            Dictionary mapping pigment names to lists of (atom_name, q_00, q_11) tuples
        """
        site_data = {}
        for name, pigment in self.pigments.items():
            site_data[name] = pigment.get_site_energy_atoms()
        return site_data
    
    def get_coupling_data(self) -> Dict[str, List[Tuple[str, float]]]:
        """
        Get coupling data for all pigments.
        
        Returns:
            Dictionary mapping pigment names to lists of (atom_name, charge) tuples
        """
        coupling_data = {}
        for name, pigment in self.pigments.items():
            coupling_data[name] = pigment.get_coupling_atoms()
        return coupling_data
    
    def get_transition_dipoles(self, transition: str = 'qy') -> Dict[str, np.ndarray]:
        """
        Get transition dipole vectors for all pigments.
        
        Args:
            transition: Type of transition ('qy', 'qx', etc.)
            
        Returns:
            Dictionary mapping pigment names to dipole vectors
        """
        dipoles = {}
        for name, pigment in self.pigments.items():
            dipoles[name] = pigment.get_transition_dipole_vector(transition)
        return dipoles
    
    def get_vacuum_energies(self) -> Dict[str, float]:
        """
        Get vacuum energies for all pigments.
        
        Returns:
            Dictionary mapping pigment names to vacuum energies
        """
        energies = {}
        for name, pigment in self.pigments.items():
            energies[name] = pigment.get_vacuum_energy()
        return energies
    
    def set_all_charge_states(self, state: str = 'ground') -> None:
        """
        Set the charge state for all pigments in the system.
        
        Args:
            state: 'ground', 'excited', or 'transition'
        """
        for pigment in self.pigments.values():
            pigment.set_charge_state(state)
    
    # ==========================================
    # Validation and analysis methods
    # ==========================================
    
    def validate_system(self) -> Dict[str, Any]:
        """
        Validate the entire pigment system for calculations.
        
        Returns:
            Dictionary with comprehensive validation results
        """
        validation = {
            'total_pigments': len(self.pigments),
            'pigment_types': {},
            'calculation_readiness': {
                'site_energy_ready': 0,
                'coupling_ready': 0,
                'dipole_ready': 0,
                'fully_ready': 0
            },
            'pigment_validations': {},
            'background_atoms': len(self.get_protein_background_atoms()),
            'system_valid': True,
            'warnings': []
        }
        
        # Count pigment types
        for pigment_type in self._pigment_types.values():
            type_name = pigment_type.__name__
            validation['pigment_types'][type_name] = validation['pigment_types'].get(type_name, 0) + 1
        
        # Validate each pigment
        for name, pigment in self.pigments.items():
            pig_validation = pigment.validate_for_calculations()
            validation['pigment_validations'][name] = pig_validation
            
            # Count readiness
            if pig_validation.get('site_energy_ready', False):
                validation['calculation_readiness']['site_energy_ready'] += 1
            if pig_validation.get('coupling_ready', False):
                validation['calculation_readiness']['coupling_ready'] += 1
            if pig_validation.get('dipole_ready', False):
                validation['calculation_readiness']['dipole_ready'] += 1
            
            # FIXED: For system validity, only require site_energy_ready for site energy calculations
            # Don't require coupling_ready unless we're doing coupling calculations
            site_energy_ready = pig_validation.get('site_energy_ready', False)
            structure_complete = pig_validation.get('structure_complete', False)
            
            if site_energy_ready and structure_complete:
                validation['calculation_readiness']['fully_ready'] += 1
            else:
                # Only mark system as invalid if site energy requirements are not met
                if not site_energy_ready:
                    validation['system_valid'] = False
                    validation['warnings'].append(f"{name}: site_energy_ready failed")
                if not structure_complete:
                    validation['warnings'].append(f"{name}: structure_complete failed")
                
                # Add coupling warning only if coupling is actually needed
                # For now, we'll just warn but not fail the system validation
                if not pig_validation.get('coupling_ready', False):
                    validation['warnings'].append(f"{name}: coupling_ready failed (warning only)")
        
        return validation
    
    def get_system_summary(self) -> Dict[str, Any]:
        """
        Get a comprehensive summary of the pigment system.
        
        Returns:
            Dictionary with system statistics and properties
        """
        validation = self.validate_system()
        positions = self.get_pigment_positions()
        
        # Calculate system geometry
        if positions:
            coords = np.array(list(positions.values()))
            center_of_mass = np.mean(coords, axis=0)
            
            # Calculate system extent
            distances_from_com = [np.linalg.norm(pos - center_of_mass) for pos in coords]
            max_extent = max(distances_from_com) if distances_from_com else 0.0
            
            # Calculate pairwise distances
            pairwise_distances = []
            coords_list = list(coords)
            for i, pos1 in enumerate(coords_list):
                for pos2 in coords_list[i+1:]:
                    pairwise_distances.append(np.linalg.norm(pos1 - pos2))
            
            geometry_stats = {
                'center_of_mass': center_of_mass.tolist(),
                'max_extent': max_extent,
                'mean_pairwise_distance': np.mean(pairwise_distances) if pairwise_distances else 0.0,
                'min_pairwise_distance': np.min(pairwise_distances) if pairwise_distances else 0.0,
                'max_pairwise_distance': np.max(pairwise_distances) if pairwise_distances else 0.0
            }
        else:
            geometry_stats = {}
        
        summary = {
            'validation': validation,
            'geometry': geometry_stats,
            'energy_range': self._get_energy_range(),
            'dipole_alignment': self._analyze_dipole_alignment(),
        }
        
        return summary
    
    def _get_energy_range(self) -> Dict[str, float]:
        """Get the range of vacuum energies in the system."""
        energies = list(self.get_vacuum_energies().values())
        if energies:
            return {
                'min_energy': min(energies),
                'max_energy': max(energies),
                'energy_spread': max(energies) - min(energies)
            }
        return {}
    
    def _analyze_dipole_alignment(self) -> Dict[str, Any]:
        """Analyze the alignment of transition dipoles."""
        qy_dipoles = self.get_transition_dipoles('qy')
        
        if not qy_dipoles:
            return {}
        
        # Calculate pairwise angles between Qy dipoles
        qy_vectors = [vec for vec in qy_dipoles.values() if np.linalg.norm(vec) > 0]
        angles = []
        
        for i, vec1 in enumerate(qy_vectors):
            for vec2 in qy_vectors[i+1:]:
                cos_angle = np.clip(np.dot(vec1, vec2), -1.0, 1.0)
                angle = np.degrees(np.arccos(cos_angle))
                angles.append(angle)
        
        if angles:
            return {
                'mean_qy_angle': np.mean(angles),
                'std_qy_angle': np.std(angles),
                'min_qy_angle': np.min(angles),
                'max_qy_angle': np.max(angles),
                'num_comparisons': len(angles)
            }
        
        return {}
    
    # ==========================================
    # Background environment methods
    # ==========================================
    
    def get_protein_background_atoms(self) -> List:
        """
        Get all protein atoms that are not part of any pigment and have charges.
        
        Returns:
            List of background atoms with charges
        """
        return self.structure.get_background_atoms(self.all_pigment_atoms)
    
    def analyze_electrostatic_environment(self, radius: float = 15.0) -> Dict[str, Dict[str, Any]]:
        """
        Analyze the electrostatic environment around each pigment.
        
        Args:
            radius: Radius for considering nearby atoms (Å)
            
        Returns:
            Dictionary with environmental analysis for each pigment
        """
        background_atoms = self.get_protein_background_atoms()
        environments = {}
        
        for name, pigment in self.pigments.items():
            pigment_pos = pigment.get_coord()
            nearby_atoms = []
            
            for atom in background_atoms:
                distance = np.linalg.norm(pigment_pos - atom.coord)
                if distance <= radius:
                    charge = getattr(atom, 'charge', 0.0)
                    nearby_atoms.append((distance, charge, atom))
            
            # Calculate electrostatic potential
            potential = 0.0
            total_charge = 0.0
            for dist, charge, atom in nearby_atoms:
                if dist > 0.1:  # Avoid division by zero
                    potential += charge / dist
                    total_charge += charge
            
            environments[name] = {
                'nearby_atoms': len(nearby_atoms),
                'total_nearby_charge': total_charge,
                'electrostatic_potential': potential,
                'mean_distance': np.mean([dist for dist, _, _ in nearby_atoms]) if nearby_atoms else 0.0,
                'charge_distribution': {
                    'positive': sum(1 for _, charge, _ in nearby_atoms if charge > 0.1),
                    'negative': sum(1 for _, charge, _ in nearby_atoms if charge < -0.1),
                    'neutral': sum(1 for _, charge, _ in nearby_atoms if abs(charge) <= 0.1)
                }
            }
        
        return environments
    
    # ==========================================
    # BioPython-style iteration and access
    # ==========================================
    
    def __getitem__(self, key: str) -> AbstractPigment:
        """Get a pigment by name (BioPython-style)."""
        return self.pigments[key]
    
    def __contains__(self, key: str) -> bool:
        """Check if pigment exists (BioPython-style)."""
        return key in self.pigments
    
    def __len__(self) -> int:
        """Number of pigments in the system."""
        return len(self.pigments)
    
    def __iter__(self):
        """Iterate over pigments (BioPython-style)."""
        return iter(self.pigments.values())
    
    def keys(self):
        """Get pigment names (dict-like interface)."""
        return self.pigments.keys()
    
    def values(self):
        """Get pigment objects (dict-like interface)."""
        return self.pigments.values()
    
    def items(self):
        """Get (name, pigment) pairs (dict-like interface)."""
        return self.pigments.items()
    
    # ==========================================
    # Export and utility methods
    # ==========================================
    
    def export_pigment_data(self) -> Dict[str, Any]:
        """
        Export comprehensive pigment data for analysis or storage.
        
        Returns:
            Dictionary with all pigment information
        """
        export_data = {
            'system_info': {
                'total_pigments': len(self.pigments),
                'pigment_types': {name: type(pig).__name__ for name, pig in self.pigments.items()},
                'structure_name': self.structure.name
            },
            'pigments': {},
            'system_summary': self.get_system_summary()
        }
        
        for name, pigment in self.pigments.items():
            export_data['pigments'][name] = {
                'type': type(pigment).__name__,
                'resname': pigment.get_resname(),
                'position': pigment.get_coord().tolist(),
                'vacuum_energy': pigment.get_vacuum_energy(),
                'site_energy_atoms': pigment.get_site_energy_atoms(),
                'coupling_atoms': pigment.get_coupling_atoms(),
                'qy_dipole': pigment.get_transition_dipole_vector('qy').tolist(),
                'qx_dipole': pigment.get_transition_dipole_vector('qx').tolist(),
                'validation': pigment.validate_for_calculations(),
                'charge_summary': pigment.get_charge_summary()
            }
        
        return export_data
    
    def __repr__(self) -> str:
        """Return a human-readable summary of the pigment system."""
        type_counts = {}
        for pigment_type in self._pigment_types.values():
            type_name = pigment_type.__name__
            type_counts[type_name] = type_counts.get(type_name, 0) + 1
        
        type_str = ", ".join([f"{count} {ptype}" for ptype, count in type_counts.items()])
        return f"<PigmentSystem with {len(self.pigments)} pigments: {type_str}>"


# Backward compatibility - keep old class name as alias
# This maintains compatibility with existing imports
Chlorophyll = ChlorophyllA

