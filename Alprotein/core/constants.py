"""
Physical constants and conversion factors used in the calculations.
"""

import numpy as np
from numpy import pi

# Physical constants for energy calculations in cm^-1
TRESP_CC = 1.1615E5  # Units: (Angstrom*cm^-1)/e^2
DIPOLE_CC = 1 / (4 * pi * 1.5812E-5)  # Units: (cm^-1 * Angstrom^3)/D^2
DEBYE_TO_E_A = 0.208194  # Converts Debye to electron-Angstroms

# Default calculation parameters
DEFAULT_DIELECTRIC_CDC = 2.0
DEFAULT_DIELECTRIC_TRESP = 1.0
DEFAULT_F_VAL = 0.72

# Default vacuum energies (cm^-1)
DEFAULT_E_VAC = {
    'CLA': 14900,  # Chlorophyll a
    'CHL': 15674,  # Chlorophyll b  
}
