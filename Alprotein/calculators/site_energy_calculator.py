"""
Site energy calculations using CDC method.
"""

import logging
import numpy as np
from typing import TYPE_CHECKING, Dict, List, Tuple, Any
from ..core.constants import TRESP_CC

logger = logging.getLogger(__name__)
if not logger.hasHandlers():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")

if TYPE_CHECKING:
    from ..core.abstract_pigments import AbstractPigment
    from ..core.pigment_system import PigmentSystem


class SiteEnergyCalculator:
    """
    Calculate site energies for pigments using the CDC (Charge Density Coupling) method.
    
    This class handles the electrostatic interactions between pigments and protein
    environment to compute site-specific transition energies.
    """
    
    def __init__(self, dielectric_constant=1.0, e0a=0.0, e0b=0.0):
        """Initialize calculator with automatic optimization detection."""
        self.dielectric_constant = dielectric_constant
        self.dielectric = dielectric_constant  # Backward compatibility
        self.e0a = e0a
        self.e0b = e0b
        self.use_direct_access = True  # Will be set automatically

    def calculate_site_energy(self, target_pigment, pigment_system):
        """Calculate site energy with automatic optimization.
        
        This method automatically detects if enhanced atoms are available
        and uses the optimized calculation path when possible.
        """
        # Detect if we can use optimized path
        self.use_direct_access = self._can_use_direct_access(target_pigment)
        
        if self.use_direct_access:
            return self._calculate_site_energy_optimized(target_pigment, pigment_system)
        else:
            return self._calculate_site_energy_legacy(target_pigment, pigment_system)

    def _can_use_direct_access(self, pigment):
        """Check if pigment has enhanced atoms with attached parameters."""
        for atom in pigment.atoms:
            if hasattr(atom, 'calculation_ready') and atom.calculation_ready:
                return True
        return False

    def _calculate_site_energy_optimized(self, target_pigment, pigment_system):
        """Optimized site energy calculation using direct atom access."""
        from ..data.parameters import get_vacuum_energy
        
        # Get the appropriate vacuum energy based on pigment type
        pigment_resname = target_pigment.get_resname().upper()
        if pigment_resname == 'CLA' and self.e0a != 0.0:  # Use E₀a as total vacuum energy for CLA
            vacuum_energy = self.e0a
        elif pigment_resname == 'CHL' and self.e0b != 0.0:  # Use E₀b as total vacuum energy for CHL
            vacuum_energy = self.e0b
        else:
            # Fallback to pigment's default vacuum energy
            vacuum_energy = target_pigment.get_vacuum_energy()
        
        # Calculate protein interaction using direct access
        delta_E_protein = self._pigment_background_interaction_optimized(
            target_pigment, pigment_system.get_protein_background_atoms()
        )
        
        # Calculate pigment-pigment interactions using direct access
        delta_E_pigments = 0.0
        for other_pigment in pigment_system.pigments.values():
            if other_pigment != target_pigment:
                delta_E_pigments += self._pigment_pigment_interaction_optimized(
                    target_pigment, other_pigment
                )
        
        return vacuum_energy + delta_E_protein + delta_E_pigments

    def _pigment_background_interaction_optimized(self, pigment, background_atoms):
        """Optimized pigment-background interaction using correct CDC equation and direct access."""
        # Get site energy atoms with direct access
        site_atoms_data = pigment.get_site_energy_atoms()
        if not site_atoms_data:
            return 0.0
        
        try:
            # Extract data for pigment atoms
            atom_names = [data[0] for data in site_atoms_data]
            q_00_charges = np.array([data[1] for data in site_atoms_data])
            q_11_charges = np.array([data[2] for data in site_atoms_data])
            
            # Calculate charge differences
            delta_q = q_11_charges - q_00_charges
            
            # Get positions
            target_pos = np.array([pigment.get_atom_coord(name) for name in atom_names])
            bg_pos = np.array([atom.coord for atom in background_atoms])
            bg_charges = np.array([getattr(atom, 'charge', 0.0) for atom in background_atoms])
            
            # Create tiled arrays for vectorized CDC calculation
            Q2_pigment_tile = np.tile(delta_q.reshape(-1, 1), (1, len(background_atoms)))
            Q2_background_tile = np.tile(bg_charges.reshape(1, -1), (len(delta_q), 1))
            
            # Calculate distance matrix
            R2_ab_norm = np.linalg.norm(
                target_pos[:, np.newaxis, :] - bg_pos[np.newaxis, :, :], 
                axis=2
            )
            
            # Avoid division by zero
            R2_ab_norm = np.where(R2_ab_norm > 0, R2_ab_norm, np.inf)
            
            # Apply correct CDC equation
            total_interaction = (TRESP_CC / self.dielectric_constant) * np.sum(
                Q2_pigment_tile * Q2_background_tile / R2_ab_norm
            )
            
            return total_interaction
            
        except Exception as e:
            logger.warning(f"Optimized background interaction failed: {e}")
            return 0.0

    def _pigment_pigment_interaction_optimized(self, target_pigment, other_pigment):
        """Optimized pigment-pigment interaction using correct CDC equation and direct atom access."""
        # Get site energy atoms for target pigment
        target_site_data = target_pigment.get_site_energy_atoms()
        other_site_data = other_pigment.get_site_energy_atoms()
        
        if not target_site_data or not other_site_data:
            return 0.0
        
        try:
            # Target pigment data (transition charges)
            target_names = [data[0] for data in target_site_data]
            target_q00 = np.array([data[1] for data in target_site_data])
            target_q11 = np.array([data[2] for data in target_site_data])
            target_delta_q = target_q11 - target_q00
            target_pos = np.array([target_pigment.get_atom_coord(name) for name in target_names])
            
            # Other pigment data (ground state charges)
            other_names = [data[0] for data in other_site_data]
            other_q00 = np.array([data[1] for data in other_site_data])
            other_pos = np.array([other_pigment.get_atom_coord(name) for name in other_names])
            
            # Create tiled arrays for vectorized CDC calculation
            Q2_target_tile = np.tile(target_delta_q.reshape(-1, 1), (1, len(other_q00)))
            Q2_other_tile = np.tile(other_q00.reshape(1, -1), (len(target_delta_q), 1))
            
            # Calculate distance matrix
            R2_ab_norm = np.linalg.norm(
                target_pos[:, np.newaxis, :] - other_pos[np.newaxis, :, :], 
                axis=2
            )
            
            # Avoid division by zero
            R2_ab_norm = np.where(R2_ab_norm > 0, R2_ab_norm, np.inf)
            
            # Apply correct CDC equation
            total_interaction = (TRESP_CC / self.dielectric_constant) * np.sum(
                Q2_target_tile * Q2_other_tile / R2_ab_norm
            )
            
            return total_interaction
            
        except Exception as e:
            logger.warning(f"Optimized pigment-pigment interaction failed: {e}")
            return 0.0

    def _calculate_site_energy_legacy(self, target_pig: 'AbstractPigment', 
                                    pigment_system: 'PigmentSystem') -> float:
        """
        Legacy site energy calculation method with corrected vacuum energy usage."""
        # Get the appropriate vacuum energy based on pigment type
        pigment_resname = target_pig.get_resname().upper()
        if pigment_resname == 'CLA' and self.e0a != 0.0:  # Use E₀a as total vacuum energy for CLA
            energy = self.e0a
        elif pigment_resname == 'CHL' and self.e0b != 0.0:  # Use E₀b as total vacuum energy for CHL
            energy = self.e0b
        else:
            # Fallback to pigment's default vacuum energy
            energy = target_pig.get_vacuum_energy()
        
        # Get protein background
        protein_bkg_atoms = pigment_system.get_protein_background_atoms()
        protein_bkg_charges = [getattr(atom, 'charge', 0.0) for atom in protein_bkg_atoms]
        
        # 1. Pigment-Protein interaction
        energy += self._pigment_background_interaction(
            target_pig, protein_bkg_atoms, protein_bkg_charges
        )
        
        # 2. Pigment-Pigment interactions
        for other_pig in pigment_system.values():
            if target_pig.get_id() == other_pig.get_id():
                continue
            
            energy += self._pigment_pigment_interaction(target_pig, other_pig)
        
        return energy
    
    def _pigment_background_interaction_detailed(self, target_pig: 'AbstractPigment', 
                                               background_atoms: List, 
                                               background_charges: List[float]) -> Tuple[float, Dict[str, float]]:
        """
        Calculate electrostatic interaction with detailed atom-by-atom breakdown using correct CDC equation.
        
        CDC Equation: (TRESP_CC / dielectric_eff) * np.sum(Q2_pigA_tile * Q2_background_tile / R2_ab_norm)
        where TRESP_CC = 1.1615E5 # Units: Angstrom*cm^-1)/e^2
        
        Args:
            target_pig: Target pigment for site energy calculation
            background_atoms: List of background atoms
            background_charges: List of background charges
            
        Returns:
            Tuple of (total_interaction_energy, atom_contributions_dict)
        """
        # Get site energy atoms from the abstract interface
        site_atoms_data = target_pig.get_site_energy_atoms()
        if not site_atoms_data:
            return 0.0, {}
        
        try:
            # Extract atom names, positions, and charge differences
            atom_names = [data[0] for data in site_atoms_data]
            q_00_charges = np.array([data[1] for data in site_atoms_data])
            q_11_charges = np.array([data[2] for data in site_atoms_data])
            
            # Calculate charge differences (delta_q = q11 - q00)
            delta_q = q_11_charges - q_00_charges
            
            # Get pigment atom positions
            target_pos = np.array([target_pig.get_atom_coord(name) for name in atom_names])
            
            # Calculate interaction energy with detailed breakdown
            atom_contributions = {}
            total_interaction = 0.0
            
            for j, (bg_atom, q_j) in enumerate(zip(background_atoms, background_charges)):
                bg_coord = bg_atom.coord
                atom_name = bg_atom.name
                
                # Calculate distances from this background atom to all pigment atoms
                distances = np.linalg.norm(target_pos - bg_coord, axis=1)
                
                # Avoid division by zero
                valid_distances = distances > 0
                
                if np.any(valid_distances):
                    # Apply CDC equation for this background atom
                    # Sum over all pigment atoms interacting with this background atom
                    atom_interaction = np.sum(
                        (delta_q[valid_distances] * q_j) / distances[valid_distances]
                    )
                    
                    # Convert to cm⁻¹ and store
                    atom_contribution = (TRESP_CC / self.dielectric) * atom_interaction
                    
                    if abs(atom_contribution) > 1e-6:  # Only store significant contributions
                        atom_contributions[atom_name] = atom_contribution
                        total_interaction += atom_interaction
            
            # Convert total to cm⁻¹
            total_energy = (TRESP_CC / self.dielectric) * total_interaction
            return total_energy, atom_contributions
            
        except (KeyError, ValueError) as e:
            logger.warning(
                "Detailed CDC calculation failed for %s: %s",
                target_pig.get_id(),
                e,
            )
            return 0.0, {}
    
    def _pigment_background_interaction(self, target_pig: 'AbstractPigment', 
                                      background_atoms: List, 
                                      background_charges: List[float]) -> float:
        """
        Calculate electrostatic interaction between a pigment and background charges using correct CDC equation.
        
        CDC Equation: (TRESP_CC / dielectric_eff) * np.sum(Q2_pigA_tile * Q2_background_tile / R2_ab_norm)
        where TRESP_CC = 1.1615E5 # Units: Angstrom*cm^-1)/e^2
        
        Args:
            target_pig: Target pigment for site energy calculation
            background_atoms: List of background atoms
            background_charges: List of background charges
            
        Returns:
            Interaction energy in cm⁻¹
        """
        # Get site energy atoms from the abstract interface
        site_atoms_data = target_pig.get_site_energy_atoms()
        if not site_atoms_data:
            return 0.0
        
        try:
            # Extract atom names, positions, and charge differences
            atom_names = [data[0] for data in site_atoms_data]
            q_00_charges = np.array([data[1] for data in site_atoms_data])
            q_11_charges = np.array([data[2] for data in site_atoms_data])
            
            # Calculate charge differences (delta_q = q11 - q00)
            delta_q = q_11_charges - q_00_charges
            
            # Get pigment atom positions
            target_pos = np.array([target_pig.get_atom_coord(name) for name in atom_names])
            
            # Get background positions and charges
            background_pos = np.array([atom.coord for atom in background_atoms])
            background_charges_array = np.array(background_charges)
            
            # Create tiled arrays for vectorized calculation
            # Q2_pigA_tile: (n_pig_atoms, n_bg_atoms) - each row has the same delta_q value
            Q2_pigA_tile = np.tile(delta_q.reshape(-1, 1), (1, len(background_atoms)))
            
            # Q2_background_tile: (n_pig_atoms, n_bg_atoms) - each column has the same bg charge
            Q2_background_tile = np.tile(background_charges_array.reshape(1, -1), (len(delta_q), 1))
            
            # Calculate distance matrix R2_ab_norm
            # Shape: (n_pig_atoms, n_bg_atoms)
            R2_ab_norm = np.linalg.norm(
                target_pos[:, np.newaxis, :] - background_pos[np.newaxis, :, :], 
                axis=2
            )
            
            # Avoid division by zero
            R2_ab_norm = np.where(R2_ab_norm > 0, R2_ab_norm, np.inf)
            
            # Apply correct CDC equation
            interaction_energy = (TRESP_CC / self.dielectric) * np.sum(
                Q2_pigA_tile * Q2_background_tile / R2_ab_norm
            )
            
            return interaction_energy
            
        except (KeyError, ValueError) as e:
            logger.warning(
                "CDC calculation failed for %s: %s",
                target_pig.get_id(),
                e,
            )
            return 0.0
    
    def _pigment_pigment_interaction(self, target_pig: 'AbstractPigment',
                                   other_pig: 'AbstractPigment') -> float:
        """
        Calculate pigment-pigment electrostatic interaction.
        
        Args:
            target_pig: Target pigment
            other_pig: Other pigment providing the electrostatic field
            
        Returns:
            Interaction energy in cm⁻¹
        """
        # Get other pigment's ground state charges
        other_site_data = other_pig.get_site_energy_atoms()
        if not other_site_data:
            return 0.0
        
        try:
            # Create pseudo-atoms list from other pigment's charged atoms
            other_atoms = []
            other_charges = []
            
            for atom_name, q_00, q_11 in other_site_data:
                try:
                    atom = other_pig.residue[atom_name]
                    other_atoms.append(atom)
                    other_charges.append(q_00)  # Use ground state charges
                except KeyError:
                    continue
            
            # Calculate interaction using the background interaction method
            return self._pigment_background_interaction(target_pig, other_atoms, other_charges)
            
        except (KeyError, ValueError) as e:
            logger.warning(
                "Pigment-pigment interaction failed for %s-%s: %s",
                target_pig.get_id(),
                other_pig.get_id(),
                e,
            )
            return 0.0
    

    def calculate_all_site_energies(self, pigment_system: 'PigmentSystem',
                                  vacuum_energies: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate site energies for all pigments in the system.
        
        Args:
            pigment_system: Complete pigment system
            vacuum_energies: Dictionary mapping residue names to vacuum energies (DEPRECATED - now uses E₀a/E₀b)
            
        Returns:
            Dictionary mapping pigment names to site energies
        """
        site_energies = {}
        
        logger.info("--- Calculating Site Energies (CDC) ---")
        logger.info("  Using E₀a = %.1f cm⁻¹ for CLA pigments", self.e0a)
        logger.info("  Using E₀b = %.1f cm⁻¹ for CHL pigments", self.e0b)
        
        for name, pigment in pigment_system.items():
            # Calculate site energy using the corrected method
            site_energy = self.calculate_site_energy(pigment, pigment_system)
            site_energies[name] = site_energy
            
            # Get the vacuum energy that was actually used for logging
            pigment_resname = pigment.get_resname().upper()
            if pigment_resname == 'CLA' and self.e0a != 0.0:
                vacuum_energy = self.e0a
                energy_source = "E₀a"
            elif pigment_resname == 'CHL' and self.e0b != 0.0:
                vacuum_energy = self.e0b
                energy_source = "E₀b"
            else:
                vacuum_energy = pigment.get_vacuum_energy()
                energy_source = "default"
            
            # Calculate the shift for logging (site energy - vacuum energy)
            energy_shift = site_energy - vacuum_energy
            
            logger.info(
                "  %s (%s): %.1f cm⁻¹ (vacuum: %.1f [%s], shift: %+.1f)",
                name,
                pigment_resname,
                site_energy,
                vacuum_energy,
                energy_source,
                energy_shift
            )
        
        return site_energies
    
    def calculate_detailed_site_energy_contributions(self, pigment_system: 'PigmentSystem') -> Tuple[Dict[str, float], Dict[str, Dict[str, float]]]:
        """
        Calculate detailed site energy contributions for each pigment, breaking down
        contributions from individual residues and other pigments.
        
        This method provides the core CDC (Chromophore-Dependent Contribution) analysis
        by calculating how each component of the protein environment contributes to
        the site energy shift of each pigment.
        
        Args:
            pigment_system: Complete pigment system
            
        Returns:
            Tuple of:
            - Dictionary mapping pigment names to total site energy shifts
            - Dictionary mapping pigment names to individual contribution breakdowns
        """
        dict_total_site_energy = {}
        dict_total_contribution_site_energy = {}
        
        logger.info("--- Calculating Detailed CDC Contributions ---")
        
        for target_pigment in pigment_system:
            pigment_name = target_pigment.get_id()
            logger.info("Processing %s...", pigment_name)
            
            dict_background_contribution = {}
            total_shift = 0.0
            
            # Get vacuum energy
            vacuum_energy = target_pigment.get_vacuum_energy()
            dict_background_contribution['vacuum'] = vacuum_energy
            
            # 1. Calculate contributions from other pigments
            for other_pigment in pigment_system:
                if target_pigment.get_id() == other_pigment.get_id():
                    continue
                
                contribution = self._pigment_pigment_interaction(target_pigment, other_pigment)
                if abs(contribution) > 1e-6:  # Only store non-zero contributions
                    dict_background_contribution[other_pigment.get_id()] = contribution
                    total_shift += contribution
            
            # 2. Calculate contributions from protein residues
            protein_atoms = pigment_system.get_protein_background_atoms()
            residue_groups = self._group_protein_atoms_by_residue(protein_atoms)
            
            for residue_id, atoms in residue_groups.items():
                charges = [getattr(atom, 'charge', 0.0) for atom in atoms]
                contribution = self._pigment_background_interaction(target_pigment, atoms, charges)
                
                if abs(contribution) > 1e-6:  # Only store non-zero contributions
                    dict_background_contribution[residue_id] = contribution
                    total_shift += contribution
            
            # Store results
            dict_total_site_energy[pigment_name] = total_shift
            dict_total_contribution_site_energy[pigment_name] = dict_background_contribution
            
            logger.info("  Total shift: %+.1f cm⁻¹", total_shift)
            logger.info(
                "  Components: %d significant contributors",
                len(
                    [
                        k
                        for k, v in dict_background_contribution.items()
                        if k != "vacuum" and abs(v) > 1.0
                    ]
                ),
            )
        
        return dict_total_site_energy, dict_total_contribution_site_energy
    
    def calculate_detailed_site_energy_contributions_with_atoms(self, pigment_system: 'PigmentSystem') -> Tuple[Dict[str, float], Dict[str, Dict[str, Any]]]:
        """
        Calculate detailed site energy contributions with atom-level breakdown for each residue.
        
        This method provides comprehensive CDC analysis including which specific atoms
        within each residue contribute most to the site energy shifts.
        
        Args:
            pigment_system: Complete pigment system
            
        Returns:
            Tuple of:
            - Dictionary mapping pigment names to total site energy shifts
            - Dictionary mapping pigment names to detailed contribution breakdowns with atom info
        """
        dict_total_site_energy = {}
        dict_detailed_contributions = {}
        
        logger.info(
            "--- Calculating Detailed CDC Contributions with Atom Breakdown ---"
        )
        
        for target_pigment in pigment_system:
            pigment_name = target_pigment.get_id()
            logger.info("Processing %s...", pigment_name)
            
            dict_background_contribution = {}
            total_shift = 0.0
            
            # Get vacuum energy
            vacuum_energy = target_pigment.get_vacuum_energy()
            dict_background_contribution['vacuum'] = {
                'total_contribution': vacuum_energy,
                'atom_breakdown': {'vacuum_energy': vacuum_energy}
            }
            
            # 1. Calculate contributions from other pigments
            for other_pigment in pigment_system:
                if target_pigment.get_id() == other_pigment.get_id():
                    continue
                
                contribution = self._pigment_pigment_interaction(target_pigment, other_pigment)
                if abs(contribution) > 1e-6:
                    dict_background_contribution[other_pigment.get_id()] = {
                        'total_contribution': contribution,
                        'atom_breakdown': {'pigment_interaction': contribution}
                    }
                    total_shift += contribution
            
            # 2. Calculate contributions from protein residues with atom breakdown
            protein_atoms = pigment_system.get_protein_background_atoms()
            residue_groups = self._group_protein_atoms_by_residue(protein_atoms)
            
            for residue_id, atoms in residue_groups.items():
                charges = [getattr(atom, 'charge', 0.0) for atom in atoms]
                
                # Get detailed atom-by-atom breakdown
                total_contribution, atom_contributions = self._pigment_background_interaction_detailed(
                    target_pigment, atoms, charges
                )
                
                if abs(total_contribution) > 1e-6:
                    # Find dominant atoms (contributing more than 10% of total or > 5 cm⁻¹)
                    significant_atoms = {}
                    dominant_threshold = max(abs(total_contribution) * 0.1, 5.0)
                    
                    for atom_name, atom_contrib in atom_contributions.items():
                        if abs(atom_contrib) > dominant_threshold:
                            significant_atoms[atom_name] = atom_contrib
                    
                    # Sort atoms by contribution magnitude
                    sorted_atoms = dict(sorted(significant_atoms.items(), 
                                             key=lambda x: abs(x[1]), reverse=True))
                    
                    dict_background_contribution[residue_id] = {
                        'total_contribution': total_contribution,
                        'atom_breakdown': atom_contributions,
                        'dominant_atoms': sorted_atoms,
                        'num_atoms': len(atoms),
                        'num_significant_atoms': len(significant_atoms)
                    }
                    total_shift += total_contribution
            
            # Store results
            dict_total_site_energy[pigment_name] = total_shift
            dict_detailed_contributions[pigment_name] = dict_background_contribution
            
            logger.info("  Total shift: %+.1f cm⁻¹", total_shift)
            logger.info(
                "  Residues: %d",
                len(
                    [
                        k
                        for k, v in dict_background_contribution.items()
                        if k != "vacuum"
                        and "total_contribution" in v
                        and abs(v["total_contribution"]) > 1.0
                    ]
                ),
            )
        
        return dict_total_site_energy, dict_detailed_contributions
    
    def _group_protein_atoms_by_residue(self, atoms: List) -> Dict[str, List]:
        """
        Group protein atoms by their residue identifier.
        
        Args:
            atoms: List of BioPython atoms
            
        Returns:
            Dictionary mapping residue IDs to lists of atoms
        """
        residue_groups = {}
        
        for atom in atoms:
            try:
                residue = atom.get_parent()
                chain = residue.get_parent()
                residue_id = f"{chain.id}_{residue.resname}_{residue.id[1]}"
                
                if residue_id not in residue_groups:
                    residue_groups[residue_id] = []
                residue_groups[residue_id].append(atom)
                
            except AttributeError:
                # Handle atoms without proper parent structure
                continue
        
        return residue_groups
    
    def analyze_energy_contributions(self, target_pig: 'AbstractPigment',
                                   pigment_system: 'PigmentSystem',
                                   vacuum_energy: float) -> Dict[str, float]:
        """
        Legacy method: Analyze individual contributions to site energy for a single pigment.
        
        This method is kept for backward compatibility. For comprehensive analysis,
        use calculate_detailed_site_energy_contributions() instead.
        
        Args:
            target_pig: Target pigment
            pigment_system: Complete pigment system
            vacuum_energy: Vacuum transition energy
            
        Returns:
            Dictionary with energy contribution breakdown
        """
        # Use the new detailed method and extract data for the target pigment
        _, all_contributions = self.calculate_detailed_site_energy_contributions(pigment_system)
        
        target_contributions = all_contributions.get(target_pig.get_id(), {})
        
        # Restructure to match legacy format
        contributions = {
            'vacuum': vacuum_energy,
            'protein_background': 0.0,
            'pigment_interactions': {},
            'total': 0.0
        }
        
        # Separate protein background and pigment interactions
        for contributor, value in target_contributions.items():
            if contributor == 'vacuum':
                continue
            elif any(pig.get_id() in contributor for pig in pigment_system):
                contributions['pigment_interactions'][contributor] = value
            else:
                contributions['protein_background'] += value
        
        # Calculate total
        contributions['total'] = (
            contributions['vacuum'] + 
            contributions['protein_background'] + 
            sum(contributions['pigment_interactions'].values())
        )
        
        return contributions
    
    def validate_site_energy_readiness(self, pigment_system: 'PigmentSystem') -> Dict[str, Dict[str, Any]]:
        """
        Validate that all pigments are ready for site energy calculations.
        
        Args:
            pigment_system: Pigment system to validate
            
        Returns:
            Dictionary with validation results for each pigment
        """
        validation_results = {}
        
        for name, pigment in pigment_system.items():
            result = {
                'ready': False,
                'site_atoms_available': 0,
                'vacuum_energy': 0.0,
                'issues': []
            }
            
            # Check site energy atoms
            site_atoms = pigment.get_site_energy_atoms()
            result['site_atoms_available'] = len(site_atoms)
            
            if len(site_atoms) == 0:
                result['issues'].append("No site energy atoms available")
            
            # Check vacuum energy
            vacuum_energy = pigment.get_vacuum_energy()
            result['vacuum_energy'] = vacuum_energy
            
            if vacuum_energy == 0.0:
                result['issues'].append("No vacuum energy specified")
            
            # Check overall validation from pigment
            pig_validation = pigment.validate_for_calculations()
            if not pig_validation.get('site_energy_ready', False):
                result['issues'].append("Pigment validation failed")
            
            result['ready'] = len(result['issues']) == 0
            validation_results[name] = result
        
        return validation_results
    
    def get_calculation_summary(self, pigment_system: 'PigmentSystem') -> Dict[str, Any]:
        """
        Get a summary of the site energy calculation setup.
        
        Args:
            pigment_system: Pigment system
            
        Returns:
            Summary dictionary
        """
        validation = self.validate_site_energy_readiness(pigment_system)
        background_atoms = pigment_system.get_protein_background_atoms()
        
        ready_count = sum(1 for v in validation.values() if v['ready'])
        total_site_atoms = sum(v['site_atoms_available'] for v in validation.values())
        
        summary = {
            'calculator_setup': {
                'dielectric': self.dielectric,
                'method': 'CDC (Charge Density Coupling)'
            },
            'system_readiness': {
                'total_pigments': len(pigment_system),
                'ready_pigments': ready_count,
                'total_site_atoms': total_site_atoms,
                'background_atoms': len(background_atoms)
            },
            'vacuum_energies': pigment_system.get_vacuum_energies(),
            'validation_details': validation
        }
        
        return summary
    
    def set_parameters(self, dielectric: float, e0a: float, e0b: float) -> None:
        """
        Update dielectric constant and E0a, E0b for calculations.
        
        Args:
            dielectric: New dielectric constant
            e0a: New E0a value
            e0b: New E0b value
        """
        self.dielectric = dielectric
        self.dielectric_constant = dielectric  # Keep both in sync
        self.e0a = e0a
        self.e0b = e0b
    
    def perform_cdc_analysis(self, pigment_system: 'PigmentSystem',
                           pdb_path: str,
                           output_dir: str,
                           cutoff: float = 20.0,
                           font_size: int = 8,
                           show_values: bool = True) -> Dict[str, Dict[str, float]]:
        """
        Perform comprehensive CDC (Chromophore-Dependent Contribution) analysis.
        
        This method provides a convenient interface to the CDC analysis tools
        directly from the site energy calculator.
        
        Args:
            pigment_system: The pigment system to analyze
            pdb_path: Path to the PDB file for spatial visualization
            output_dir: Directory to save results
            cutoff: Minimum contribution threshold (cm⁻¹)
            font_size: Font size for plots
            show_values: Whether to show values on spatial plots
            
        Returns:
            Dictionary of all calculated contributions
        """
        from ..utils.cdc_analysis import analyze_site_energy_contributions
        
        return analyze_site_energy_contributions(
            pigment_system=pigment_system,
            site_energy_calculator=self,
            pdb_path=pdb_path,
            output_dir=output_dir,
            cutoff=cutoff,
            font_size=font_size,
            show_values=show_values
        )
    
    def __repr__(self) -> str:
        """Return a brief representation of the site energy calculator."""
        return f"SiteEnergyCalculator(dielectric={self.dielectric}, e0a={self.e0a}, e0b={self.e0b})"
