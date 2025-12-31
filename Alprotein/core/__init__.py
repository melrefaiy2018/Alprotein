"""
Core module initialization.
"""

from .constants import *
from .protein_structure import ProteinStructure
from .abstract_pigments import (
    AbstractPigment, ChlorophyllA, ChlorophyllB, Pheophytin, 
    create_pigment, Chlorophyll  # Backward compatibility alias
)
from .pigment_system import PigmentSystem

__all__ = [
    # Constants
    'TRESP_CC', 'DIPOLE_CC', 'DEBYE_TO_E_A',
    'DEFAULT_DIELECTRIC_CDC', 'DEFAULT_DIELECTRIC_TRESP', 'DEFAULT_F_VAL', 'DEFAULT_E_VAC',
    
    # Core classes
    'ProteinStructure', 'PigmentSystem', 
    
    # Abstract pigment interface
    'AbstractPigment',
    
    # Concrete pigment implementations
    'ChlorophyllA', 'ChlorophyllB', 'Pheophytin',
    
    # Factory function
    'create_pigment',
    
    # Backward compatibility
    'Chlorophyll'  # Alias for ChlorophyllA
]
