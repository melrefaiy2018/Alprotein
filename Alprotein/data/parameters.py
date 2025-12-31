"""
Parameter data storage and access functions.

This module contains the parameter dictionaries for coupling calculations
(TrEsp) and site energy calculations (CDC), along with functions to access them.
"""

from ..core.constants import DEBYE_TO_E_A

# Coupling parameters for TrEsp calculations
COUPLING_DATA = {
    'CLA_IPPC': {
        'tresp_atoms': ['CAA', 'CAB', 'CAC', 'CAD', 'N1A', 'CBA', 'CBB', 'CBC', 'CBD', 'N1B',
                        'OBD', 'CGA', 'CGD', 'N1D', 'CED', 'CHA', 'CHB', 'CHC', 'CHD',
                        'CMA', 'CMB', 'CMC', 'CMD', 'C1A', 'C1B', 'C1C', 'C1D', 'O1A',
                        'O1D', 'C2A', 'C2B', 'C2C', 'C2D', 'O2A', 'O2D', 'C3A', 'C3B',
                        'C3C', 'C3D', 'C4A', 'C4B', 'C4C', 'C4D', 'N1C', 'MG'],
        'tresp_pc': [-0.0010500000000000002, 0.010681, 0.00802, -0.01925, 0.031683,
                     0.000735, 0.034947000000000006, 0.001012, -0.011238, -0.062297,
                     -0.020044, -0.008039999999999999, 0.005979, 0.108292, -0.005256,
                     0.106779, -0.048695999999999996, -0.098725, 0.072726, 0.005556,
                     0.016963000000000002, -0.005161, -0.025155999999999998, -0.13082,
                     0.081122, 0.083646, -0.110812, -0.00133, -0.0053879999999999996,
                     0.010048, 0.004777, -0.0074199999999999995, -0.011981, 0.005667,
                     0.001811, 0.0023889999999999996, -0.009204, -0.0011259999999999998,
                     0.008799, 0.07798300000000001, 0.106271, -0.044008000000000005,
                     -0.125044, -0.012166, -0.021674],
        'vacuum_mag': 4.5 * DEBYE_TO_E_A,
        'dipole_mag': 4.5,
    },
    'CHL_IPPC': {
        'tresp_atoms': ['MG', 'CHA', 'CHB', 'CHC', 'CHD', 'N1A', 'C1A', 'C2A', 'C3A',
                        'C4A', 'CMA', 'CAA', 'CBA', 'CGA', 'O1A', 'O2A', 'N1B', 'C1B',
                        'C2B', 'C3B', 'C4B', 'CMB', 'CAB', 'CBB', 'N1C', 'C1C', 'C2C',
                        'C3C', 'C4C', 'CMC', 'OMC', 'CAC', 'CBC', 'N1D', 'C1D', 'C2D', 'C3D',
                        'C4D', 'CMD', 'CAD', 'OBD', 'CBD', 'CGD', 'O1D', 'O2D', 'CED'],
        'tresp_pc': [-0.015512, 0.075628, -0.057844, -0.107433, 0.086081, 0.011574, -0.09967, 0.007454, 0.007749,
                     0.077389, 0.003351, -0.00272, -0.005754, 0.004859, 0.00065, -0.003916, -0.07399, 0.085402,
                     -0.00214, -0.005618, 0.108571, 0.014838, 0.001978, 0.024357, -0.02165, 0.079538, 0.005163, -0.005803,
                     -0.031957, -0.013413, -0.001389, 0.006975, 0.002317, 0.090881, -0.08817, -0.013166, 0.015953, -0.098277,
                     -0.016435, -0.028234, -0.014273, -3.5e-05, 0.007743, -0.007166, 0.000641, -0.004526],
        'vacuum_mag': 3.6 * DEBYE_TO_E_A,
        'dipole_mag': 3.6,
    },
}

# Site energy parameters for CDC calculations
SITE_ENERGY_DATA = {
    'CLA': {
        'atom': ['MG', 'CHA', 'CHB', 'HHB', 'CHC', 'HHC', 'CHD', 'HHD', 'N1A',
                 'C1A', 'C2A', 'H2A', 'C3A', 'H3A', 'C4A', 'CMA', 'HMA1', 'HMA2',
                 'HMA3', 'CAA', 'HAA1', 'HAA2', 'CBA', 'HBA2', 'HBA1', 'CGA', 'O1A',
                 'O2A', 'N1B', 'C1B', 'C2B', 'C3B', 'C4B', 'CMB', 'HMB1', 'HMB2',
                 'HMB3', 'CAB', 'HBB', 'CBB', 'HBB1', 'HBB3', 'N1C', 'C1C', 'C2C',
                 'C3C', 'C4C', 'CMC', 'HMC1', 'HMC2', 'HMC3', 'CAC', 'HAC2', 'HAC1',
                 'CBC', 'HBC1', 'HBC2', 'HBC3', 'N1D', 'C1D', 'C2D', 'C3D', 'C4D', 'CMD',
                 'HMD1', 'HMD2', 'HMD3', 'CAD', 'OBD', 'CBD', 'HBD', 'CGD', 'O1D', 'O2D',
                 'CED', 'HED1', 'HED2', 'HED3'],
        'q_00': [1.035, 0.416, -0.557, 0.155, -0.037, 0.091, -0.075, 0.148, -0.539,
                 -0.098, -0.266, 0.121, 0.396, 0.008, 0.419, -0.62, 0.153, 0.153, 0.153,
                 0.238, -0.018, -0.018, -0.377, 0.105, 0.105, 0.732, -0.439, -0.239,
                 -0.483, 0.282, 0.056, 0.056, 0.058, -0.366, 0.116, 0.116, 0.116,
                 -0.191, 0.164, -0.36, 0.177, 0.177, -0.33, -0.122, 0.321, -0.232,
                 0.031, -0.6, 0.164, 0.164, 0.164, 0.155, 0.001, 0.001, -0.196, 0.052,
                 0.052, 0.052, -0.41, -0.011, 0.248, -0.248, -0.069, -0.488, 0.148,
                 0.148, 0.148, 0.75, -0.555, -0.869, 0.273, 0.967, -0.612, -0.338,
                 -0.251, 0.143, 0.143, 0.143],
        'q_11': [1.034, 0.442, -0.502, 0.145, -0.093, 0.095, -0.289, 0.165, -0.532,
                 -0.108, -0.264, 0.123, 0.397, 0.008, 0.399, -0.618, 0.153, 0.153,
                 0.153, 0.235, -0.017, -0.017, -0.381, 0.107, 0.107, 0.736, -0.44,
                 -0.24, -0.46, 0.205, 0.044, 0.012, 0.117, -0.353, 0.113, 0.113, 0.113,
                 -0.184, 0.164, -0.374, 0.179, 0.179, -0.416, -0.079, 0.342, -0.253,
                 0.234, -0.602, 0.166, 0.166, 0.166, 0.156, 0.0, 0.0, -0.196, 0.051,
                 0.051, 0.051, -0.451, 0.162, 0.19, -0.247, -0.046, -0.481, 0.147,
                 0.147, 0.147, 0.749, -0.555, -0.897, 0.279, 0.975, -0.614, -0.338,
                 -0.252, 0.143, 0.143, 0.143]
    },
    'CHL': {
        'atom': ['MG', 'CHA', 'CHB', 'HHB', 'CHC', 'HHC', 'CHD', 'HHD', 'N1A',
                 'C1A', 'C2A', 'H2A', 'C3A', 'H3A', 'C4A', 'CMA', 'HMA1', 'HMA2',
                 'HMA3', 'CAA', 'HAA2', 'HAA1', 'CBA', 'HBA2', 'HBA1', 'CGA', 'O1A',
                 'O2A', 'N1B', 'C1B', 'C2B', 'C3B', 'C4B', 'CMB', 'HMB1', 'HMB2',
                 'HMB3', 'CAB', 'HBB', 'CBB', 'HBB1', 'HBB3', 'N1C', 'C1C', 'C2C',
                 'C3C', 'C4C', 'CMC', 'OMC', 'HMC1', 'CAC', 'HAC2', 'HAC1', 'CBC',
                 'HBC1', 'HBC2', 'HBC3', 'N1D', 'C1D', 'C2D', 'C3D', 'C4D', 'CMD',
                 'HMD1', 'HMD2', 'HMD3', 'CAD', 'OBD', 'CBD', 'HBD', 'CGD', 'O1D',
                 'O2D', 'CED', 'HED1', 'HED2', 'HED3'],
        'q_00': [1.19, -0.051, -0.729, 0.179, -0.291, 0.223, -0.212, 0.175, -0.694, 0.32,
                 -0.593, 0.273, 0.422, -0.032, 0.573, -0.527, 0.134, 0.134, 0.134, 0.011,
                 0.075, 0.075, -0.32, 0.11, 0.11, 0.765, -0.416, -0.286, -0.594, 0.374, 0.157, 0.036,
                 0.15, -0.548, 0.159, 0.159, 0.159, -0.218, 0.153, -0.311, 0.173, 0.173,
                 -0.547, 0.241, -0.06, -0.299, 0.28, 0.452, -0.558, 0.044, 0.172, 0.01,
                 0.01, -0.308, 0.085, 0.085, 0.085, -0.613, 0.107, 0.19, -0.278, 0.263,
                 -0.477, 0.151, 0.151, 0.151, 0.763, -0.558, -0.633, 0.181, 0.895,
                 -0.604, -0.341, -0.18, 0.122, 0.122, 0.122],
        'q_11': [1.194, 0.004, -0.657, 0.168, -0.305, 0.222, -0.438, 0.192, -0.68, 0.31,
                 -0.603, 0.277, 0.426, -0.033, 0.539, -0.526, 0.134, 0.134, 0.134, 0.01,
                 0.076, 0.076, -0.316, 0.109, 0.109, 0.764, -0.415, -0.286, -0.569,
                 0.296, 0.147, -0.006, 0.199, -0.534, 0.156, 0.156, 0.156, -0.211, 0.153,
                 -0.322, 0.175, 0.175, -0.616, 0.247, -0.039, -0.304, 0.472, 0.453,
                 -0.554, 0.043, 0.172, 0.009, 0.009, -0.304, 0.084, 0.084, 0.084,
                 -0.658, 0.288, 0.112, -0.279, 0.265, -0.467, 0.149, 0.149, 0.149,
                 0.774, -0.56, -0.696, 0.194, 0.905, -0.604, -0.337, -0.183, 0.123,
                 0.123, 0.123],
    },
}


def get_coupling_parameters(param_name: str) -> dict:
    """
    Get coupling parameters for TrEsp calculations.
    
    Args:
        param_name: Name of the parameter set (e.g., 'CLA_IPPC', 'CHL_IPPC')
        
    Returns:
        Dictionary containing coupling parameters
        
    Raises:
        KeyError: If parameter set is not found
    """
    if param_name not in COUPLING_DATA:
        available = list(COUPLING_DATA.keys())
        raise KeyError(f"Coupling parameters '{param_name}' not found. Available: {available}")
    
    return COUPLING_DATA[param_name].copy()


def get_site_energy_parameters(param_name: str) -> dict:
    """
    Get site energy parameters for CDC calculations.
    
    Args:
        param_name: Name of the parameter set (e.g., 'CLA', 'CHL')
        
    Returns:
        Dictionary containing site energy parameters
        
    Raises:
        KeyError: If parameter set is not found
    """
    if param_name not in SITE_ENERGY_DATA:
        available = list(SITE_ENERGY_DATA.keys())
        raise KeyError(f"Site energy parameters '{param_name}' not found. Available: {available}")
    
    return SITE_ENERGY_DATA[param_name].copy()


def list_available_parameters():
    """
    List all available parameter sets.
    
    Returns:
        Dictionary with available coupling and site energy parameters
    """
    return {
        'coupling_parameters': list(COUPLING_DATA.keys()),
        'site_energy_parameters': list(SITE_ENERGY_DATA.keys())
    }


def add_coupling_parameters(name: str, parameters: dict) -> None:
    """
    Add new coupling parameters (for extensibility).
    
    Args:
        name: Name for the new parameter set
        parameters: Dictionary containing the parameters
    """
    required_keys = ['tresp_atoms', 'tresp_pc', 'vacuum_mag', 'dipole_mag']
    for key in required_keys:
        if key not in parameters:
            raise ValueError(f"Missing required parameter: {key}")
    
    COUPLING_DATA[name] = parameters.copy()


def add_site_energy_parameters(name: str, parameters: dict) -> None:
    """
    Add new site energy parameters (for extensibility).
    
    Args:
        name: Name for the new parameter set
        parameters: Dictionary containing the parameters
    """
    required_keys = ['atom', 'q_00', 'q_11']
    for key in required_keys:
        if key not in parameters:
            raise ValueError(f"Missing required parameter: {key}")
    
    SITE_ENERGY_DATA[name] = parameters.copy()

def get_vacuum_energy(pigment_type):
    """Get vacuum energy for a pigment type.
    
    Args:
        pigment_type (str): Type of pigment (e.g., 'ChlorophyllA')
        
    Returns:
        float: Vacuum energy in cm^-1
    """
    vacuum_energies = {
        'ChlorophyllA': 15100.0,  # Example value - replace with actual
        'ChlorophyllB': 15200.0,  # Example value - replace with actual
        # Add more pigment types as needed
    }
    return vacuum_energies.get(pigment_type, 15100.0)  # Default fallback

def create_parameter_config(pigment_residues):
    """Create default parameter configuration from residue list.
    
    Args:
        pigment_residues (list): List of residue names found in structure
        
    Returns:
        dict: Configuration mapping pigment types to parameter names
    """
    # FIXED: Use the correct parameter keys that exist in SITE_ENERGY_DATA
    default_mapping = {
        'ChlorophyllA': 'CLA',  # This maps to the actual key in SITE_ENERGY_DATA
        'ChlorophyllB': 'CHL',  # This maps to the actual key in SITE_ENERGY_DATA
        # Add more defaults as needed
    }
    
    config = {}
    for resname in pigment_residues:
        pigment_type = _map_resname_to_type(resname)
        if pigment_type in default_mapping:
            config[pigment_type] = default_mapping[pigment_type]
    
    return config

def _map_resname_to_type(resname):
    """Map PDB residue name to pigment type."""
    mapping = {
        'CLA': 'ChlorophyllA',
        'CHL': 'ChlorophyllA',
        'BCL': 'ChlorophyllB', 
        'CHD': 'ChlorophyllB',
        # Add more mappings
    }
    return mapping.get(resname, 'ChlorophyllA')
