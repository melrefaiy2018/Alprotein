"""
Coupling calculations between pigments.
"""

import logging
import numpy as np
from typing import TYPE_CHECKING, Dict, Tuple
from ..core.constants import TRESP_CC, DIPOLE_CC
from ..utils.calculations import scale_charges

logger = logging.getLogger(__name__)
if not logger.hasHandlers():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")

if TYPE_CHECKING:
    from ..core.abstract_pigments import AbstractPigment


class CouplingCalculator:
    """
    Calculate electronic coupling between pigments using different methods.
    
    This class provides methods for TrEsp and dipole-dipole coupling calculations
    using the abstract pigment interface.
    """
    
    def __init__(self, dielectric: float = 1.0, f_val: float = 0.72):
        """
        Initialize coupling calculator.
        
        Args:
            dielectric: Dielectric constant for coupling calculations
            f_val: Oscillator strength factor for TrEsp calculations
        """
        self.dielectric = dielectric
        self.f_val = f_val
    
    def calculate_tresp_coupling(self, pigA: 'AbstractPigment', pigB: 'AbstractPigment') -> float:
        """
        Calculate TrEsp coupling between two pigments.
        
        Args:
            pigA: First pigment
            pigB: Second pigment
            
        Returns:
            Coupling strength in cm⁻¹
        """
        # Get coupling atoms from the abstract interface
        couplingA_data = pigA.get_coupling_atoms()
        couplingB_data = pigB.get_coupling_atoms()
        
        if not couplingA_data or not couplingB_data:
            return 0.0
        
        try:
            # Extract atom names and charges
            atomsA, chargesA = zip(*couplingA_data)
            atomsB, chargesB = zip(*couplingB_data)
            
            # Get atomic positions
            posA = np.array([pigA.get_atom_coord(name) for name in atomsA])
            posB = np.array([pigB.get_atom_coord(name) for name in atomsB])
            
            chargesA = np.array(chargesA)
            chargesB = np.array(chargesB)
            
            # Get vacuum dipole magnitudes for scaling
            # Try to get from params, fallback to reasonable defaults
            magA = getattr(pigA, 'params', {}).get('vacuum_mag', 4.5 * 0.20819)  # Default for CLA
            magB = getattr(pigB, 'params', {}).get('vacuum_mag', 4.5 * 0.20819)
            
            # Calculate scaling factors
            scaleA = scale_charges(chargesA.tolist(), posA, magA)
            scaleB = scale_charges(chargesB.tolist(), posB, magB)
            
            # Scale partial charges
            scaled_qA = chargesA / scaleA
            scaled_qB = chargesB / scaleB
            
            # Calculate coupling
            V_ab = 0.0
            for i, q_i in enumerate(scaled_qA):
                for j, q_j in enumerate(scaled_qB):
                    dist = np.linalg.norm(posA[i] - posB[j])
                    if dist > 1e-10:  # Avoid very small distances
                        V_ab += (q_i * q_j) / dist
            
            coupling = TRESP_CC * (self.f_val / self.dielectric) * V_ab
            
            # Check for NaN or infinite values
            if not np.isfinite(coupling):
                return 0.0
            
            return coupling
            
        except (KeyError, ValueError) as e:
            logger.warning(
                "TrEsp coupling calculation failed for %s-%s: %s",
                pigA.get_id(),
                pigB.get_id(),
                e,
            )
            return 0.0
    
    def calculate_dipole_coupling(self, pigA: 'AbstractPigment', pigB: 'AbstractPigment') -> float:
        """
        Calculate dipole-dipole coupling between two pigments.
        
        Args:
            pigA: First pigment
            pigB: Second pigment
            
        Returns:
            Coupling strength in cm⁻¹
        """
        try:
            # Get transition dipole moments using the abstract interface
            dipA_dir = pigA.get_transition_dipole_vector('qy')
            dipB_dir = pigB.get_transition_dipole_vector('qy')
            
            # Check if dipole directions are valid
            if np.linalg.norm(dipA_dir) == 0 or np.linalg.norm(dipB_dir) == 0:
                return 0.0
            
            # Get dipole magnitudes - try from params, use defaults based on pigment type
            magA = getattr(pigA, 'params', {}).get('dipole_mag')
            magB = getattr(pigB, 'params', {}).get('dipole_mag')
            
            # Fallback to reasonable defaults based on pigment type
            if magA is None:
                if 'CHL' in pigA.get_resname():
                    magA = 3.6  # Chlorophyll B
                else:
                    magA = 4.5  # Chlorophyll A, Pheophytin
            
            if magB is None:
                if 'CHL' in pigB.get_resname():
                    magB = 3.6  # Chlorophyll B
                else:
                    magB = 4.5  # Chlorophyll A, Pheophytin
            
            # Create full dipole vectors
            dipA = dipA_dir * magA
            dipB = dipB_dir * magB
            
            # Get center-to-center vector
            R_vec = pigA.get_coord() - pigB.get_coord()
            R_norm = np.linalg.norm(R_vec)
            
            if R_norm < 1e-10:  # Avoid very small distances
                return 0.0
            
            # Calculate dipole-dipole interaction
            term1 = np.dot(dipA, dipB) / R_norm**3
            term2 = 3 * (np.dot(dipA, R_vec) * np.dot(dipB, R_vec)) / R_norm**5
            J = term1 - term2
            
            coupling = (DIPOLE_CC / self.dielectric) * J
            
            # Check for NaN or infinite values
            if not np.isfinite(coupling):
                return 0.0
            
            return coupling
            
        except (ValueError, AttributeError) as e:
            logger.warning(
                "Dipole coupling calculation failed for %s-%s: %s",
                pigA.get_id(),
                pigB.get_id(),
                e,
            )
            return 0.0
    
    def calculate_coupling(self, pigA: 'AbstractPigment', pigB: 'AbstractPigment', 
                          method: str = 'tresp') -> float:
        """
        Calculate coupling using specified method.
        
        Args:
            pigA: First pigment
            pigB: Second pigment
            method: Coupling method ('tresp' or 'dipole')
            
        Returns:
            Coupling strength in cm⁻¹
        """
        if method.lower() == 'tresp':
            return self.calculate_tresp_coupling(pigA, pigB)
        elif method.lower() == 'dipole':
            return self.calculate_dipole_coupling(pigA, pigB)
        else:
            raise ValueError(f"Unknown coupling method: {method}")
    
    def calculate_all_couplings(self, pigment_system, method: str = 'tresp') -> Dict[Tuple[str, str], float]:
        """
        Calculate all pairwise couplings in a pigment system.
        
        Args:
            pigment_system: PigmentSystem instance
            method: Coupling method ('tresp' or 'dipole')
            
        Returns:
            Dictionary with coupling results
        """
        pigment_names = pigment_system.get_pigment_names()
        couplings = {}
        
        for i, nameA in enumerate(pigment_names):
            for j, nameB in enumerate(pigment_names):
                if i >= j:
                    continue
                
                pigA = pigment_system[nameA]
                pigB = pigment_system[nameB]
                
                coupling = self.calculate_coupling(pigA, pigB, method)
                
                # Store both directions for convenience
                couplings[(nameA, nameB)] = coupling
                couplings[(nameB, nameA)] = coupling
        
        return couplings
    
    def validate_coupling_readiness(self, pigment_system, method: str = 'tresp') -> Dict[str, Dict[str, any]]:
        """
        Validate that pigments are ready for coupling calculations.
        
        Args:
            pigment_system: PigmentSystem instance
            method: Coupling method to validate for
            
        Returns:
            Dictionary with validation results for each pigment
        """
        validation_results = {}
        
        for name, pigment in pigment_system.items():
            result = {
                'ready': False,
                'coupling_atoms': 0,
                'dipole_ready': False,
                'issues': []
            }
            
            if method.lower() == 'tresp':
                # Check TrEsp readiness
                coupling_atoms = pigment.get_coupling_atoms()
                result['coupling_atoms'] = len(coupling_atoms)
                
                if len(coupling_atoms) == 0:
                    result['issues'].append("No coupling atoms available for TrEsp")
                
                result['ready'] = len(coupling_atoms) > 0
                
            elif method.lower() == 'dipole':
                # Check dipole readiness
                qy_dipole = pigment.get_transition_dipole_vector('qy')
                dipole_magnitude = np.linalg.norm(qy_dipole)
                result['dipole_ready'] = dipole_magnitude > 0
                
                if dipole_magnitude == 0:
                    result['issues'].append("No valid dipole vector available")
                
                result['ready'] = dipole_magnitude > 0
            
            validation_results[name] = result
        
        return validation_results
    
    def get_calculation_summary(self, pigment_system, method: str = 'tresp') -> Dict[str, any]:
        """
        Get a summary of the coupling calculation setup.
        
        Args:
            pigment_system: PigmentSystem instance
            method: Coupling method
            
        Returns:
            Summary dictionary
        """
        validation = self.validate_coupling_readiness(pigment_system, method)
        
        ready_count = sum(1 for v in validation.values() if v['ready'])
        total_pairs = len(pigment_system) * (len(pigment_system) - 1) // 2
        
        summary = {
            'calculator_setup': {
                'method': method.upper(),
                'dielectric': self.dielectric,
                'f_val': self.f_val if method.lower() == 'tresp' else 'N/A'
            },
            'system_readiness': {
                'total_pigments': len(pigment_system),
                'ready_pigments': ready_count,
                'total_coupling_pairs': total_pairs,
                'ready_pairs': ready_count * (ready_count - 1) // 2
            },
            'validation_details': validation
        }
        
        if method.lower() == 'tresp':
            total_coupling_atoms = sum(v['coupling_atoms'] for v in validation.values())
            summary['system_readiness']['total_coupling_atoms'] = total_coupling_atoms
        
        return summary
    
    def set_parameters(self, dielectric: float = None, f_val: float = None) -> None:
        """
        Update calculation parameters.
        
        Args:
            dielectric: New dielectric constant
            f_val: New oscillator strength factor
        """
        if dielectric is not None:
            self.dielectric = dielectric
        if f_val is not None:
            self.f_val = f_val
    
    def __repr__(self) -> str:
        """Return a short representation of the calculator parameters."""
        return f"CouplingCalculator(dielectric={self.dielectric}, f_val={self.f_val})"
