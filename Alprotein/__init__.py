"""
Alprotein: A fast and modern package for calculating optical spectra of proteins from PDB files.

This package provides tools for:
- Loading and processing protein structures from PDB files
- Setting up pigment systems with various chlorophyll types
- Calculating site energies using charge density coupling (CDC)
- Computing electronic couplings using transition charge (TrEsp) or dipole methods
- Constructing and diagonalizing Hamiltonian matrices
- Analyzing electronic spectra of pigment-protein complexes

Main Classes:
    ProteinStructure: Handle protein structure and charge information
    PigmentSystem: Manage collection of pigments in the protein
    HamiltonianCalculator: Main calculator for spectra calculations
    
Example:
    >>> from Alprotein import ProteinStructure, PigmentSystem, HamiltonianCalculator
    >>> from Alprotein.core.abstract_pigments import ChlorophyllA
    >>> 
    >>> # Load protein structure
    >>> protein = ProteinStructure.from_file('structure.pdb', name='MyComplex')
    >>> 
    >>> # Set up pigment system
    >>> pigment_system = PigmentSystem(protein)
    >>> pigment_system.add_pigments_by_residue(
    ...     resname='CLA', 
    ...     pigment_class=ChlorophyllA,
    ...     tresp_dict_name='CLA_IPPC',
    ...     cdc_dict_name='CLA'
    ... )
    >>> 
    >>> # Calculate Hamiltonian
    >>> calculator = HamiltonianCalculator(pigment_system)
    >>> hamiltonian = calculator.construct_hamiltonian()
    >>> eigenvalues, eigenvectors = calculator.diagonalize_hamiltonian(hamiltonian)
"""

# Core functionality
from .core.protein_structure import ProteinStructure
from .core.pigment_system import PigmentSystem
from .core.abstract_pigments import (
    AbstractPigment, 
    ChlorophyllA, 
    ChlorophyllB, 
    Pheophytin,
    create_pigment,
    Chlorophyll  # Backward compatibility alias
)
from .core.constants import (
    TRESP_CC, 
    DIPOLE_CC, 
    DEBYE_TO_E_A,
    DEFAULT_DIELECTRIC_CDC, 
    DEFAULT_DIELECTRIC_TRESP, 
    DEFAULT_F_VAL, 
    DEFAULT_E_VAC
)

# Calculators
from .calculators.hamiltonian_calculator import HamiltonianCalculator
from .calculators.site_energy_calculator import SiteEnergyCalculator
from .calculators.coupling_calculator import CouplingCalculator

# Version
__version__ = "0.2.0"

# Author information
__author__ = "Mohamed Elrefaiy, Bailey Raber, Doran Raccah"
__email__ = "melrefaiy@example.edu, braber@example.edu, draccah@example.edu"

# Package metadata
__description__ = "A fast and modern package for calculating optical spectra of proteins from PDB files"
__url__ = "https://github.com/melrefaiy2018/Alprotein-Alpha"

# GUI functionality (optional import)
try:
    from .gui import launch_gui
    _GUI_AVAILABLE = True
except ImportError:
    _GUI_AVAILABLE = False
    launch_gui = None

# Define what gets imported with "from Alprotein import *"
__all__ = [
    # Core classes
    'ProteinStructure',
    'PigmentSystem', 
    'HamiltonianCalculator',
    
    # Pigment classes
    'AbstractPigment',
    'ChlorophyllA',
    'ChlorophyllB', 
    'Pheophytin',
    'Chlorophyll',  # Alias
    'create_pigment',
    
    # Calculators
    'SiteEnergyCalculator',
    'CouplingCalculator',
    
    # Constants
    'TRESP_CC',
    'DIPOLE_CC', 
    'DEBYE_TO_E_A',
    'DEFAULT_DIELECTRIC_CDC',
    'DEFAULT_DIELECTRIC_TRESP',
    'DEFAULT_F_VAL',
    'DEFAULT_E_VAC',
    
    # Metadata
    '__version__',
    '__author__',
    '__email__',
    '__description__',
    '__url__'
]

# Add GUI to __all__ if available
if _GUI_AVAILABLE:
    __all__.append('launch_gui')

# Package-level functions for convenience
def quick_calculation(pdb_file: str, pigment_configs: list = None, **kwargs):
    """
    Quick calculation function for simple use cases.
    
    Args:
        pdb_file: Path to PDB file
        pigment_configs: List of pigment configurations
        **kwargs: Additional parameters for calculation
        
    Returns:
        Dictionary with calculation results
    """
    if pigment_configs is None:
        pigment_configs = [
            {
                'resname': 'CLA',
                'pigment_class': ChlorophyllA,
                'tresp_dict_name': 'CLA_IPPC',
                'cdc_dict_name': 'CLA'
            }
        ]
    
    # Load structure
    protein = ProteinStructure.from_file(pdb_file, name='QuickCalc')
    
    # Set up pigment system
    pigment_system = PigmentSystem(protein)
    for config in pigment_configs:
        pigment_system.add_pigments_by_residue(**config)
    
    # Calculate
    calculator = HamiltonianCalculator(pigment_system)
    hamiltonian = calculator.construct_hamiltonian(kwargs)
    eigenvalues, eigenvectors = calculator.diagonalize_hamiltonian(hamiltonian)
    
    return {
        'hamiltonian': hamiltonian,
        'eigenvalues': eigenvalues,
        'eigenvectors': eigenvectors,
        'pigment_names': calculator.pigment_names
    }


def list_available_pigments():
    """Return list of available pigment types."""
    return ['ChlorophyllA', 'ChlorophyllB', 'Pheophytin']


def get_default_parameters():
    """Return default calculation parameters."""
    return {
        'dielectric_cdc': DEFAULT_DIELECTRIC_CDC,
        'dielectric_tresp': DEFAULT_DIELECTRIC_TRESP,
        'f_val': DEFAULT_F_VAL,
        'e_vac': DEFAULT_E_VAC.copy(),
        'coupling_calc': 'tresp'
    }


# Print welcome message when imported (optional, can be disabled)
def _print_welcome():
    """Print welcome message with basic info."""
    print(f"Alprotein v{__version__}")
    print("Package for calculating optical spectra of pigment-protein complexes")
    if _GUI_AVAILABLE:
        print("GUI available - use launch_gui() to start the graphical interface")
    print("For help, try: help(Alprotein) or visit the documentation")


def gui():
    """Convenience function to launch the GUI."""
    if _GUI_AVAILABLE:
        return launch_gui()
    else:
        print("GUI not available. Please install PyQt5 and other GUI dependencies:")
        print("pip install PyQt5 matplotlib")
        return None


# Add convenience function to __all__
if _GUI_AVAILABLE:
    __all__.append('gui')


# Uncomment the line below to show welcome message on import
# _print_welcome()
