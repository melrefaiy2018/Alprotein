"""
Utility functions for calculations and data processing.
"""

import numpy as np
from typing import List, Tuple


def scale_charges(list_partial_charges: List[float], 
                 list_atom_positions: List[np.ndarray], 
                 vacuum_mag: float) -> float:
    """
    Calculate the scaling factor for transition charges.
    
    Args:
        list_partial_charges: List of partial charges
        list_atom_positions: List of atomic positions
        vacuum_mag: Vacuum magnitude for normalization
        
    Returns:
        Scaling factor
    """
    if vacuum_mag == 0:
        return 1.0  # Avoid division by zero
    
    dipole_mag = np.linalg.norm(
        sum([q * pos for q, pos in zip(list_partial_charges, list_atom_positions)])
    )
    
    # Avoid division by zero if calculated dipole is zero
    if dipole_mag == 0:
        return 1.0  # Return unity scaling to avoid division by zero
    
    return dipole_mag / vacuum_mag


def validate_atom_lists(atoms1: List[str], atoms2: List[str], 
                       charges1: List[float], charges2: List[float]) -> bool:
    """
    Validate that atom and charge lists have consistent lengths.
    
    Args:
        atoms1: First atom list
        atoms2: Second atom list
        charges1: First charge list
        charges2: Second charge list
        
    Returns:
        True if all lists are consistent
        
    Raises:
        ValueError: If lists have inconsistent lengths
    """
    if len(atoms1) != len(charges1):
        raise ValueError(f"Atom list length ({len(atoms1)}) != charge list length ({len(charges1)})")
    
    if len(atoms2) != len(charges2):
        raise ValueError(f"Atom list length ({len(atoms2)}) != charge list length ({len(charges2)})")
    
    return True


def calculate_distance_matrix(positions1: np.ndarray, positions2: np.ndarray) -> np.ndarray:
    """
    Calculate distance matrix between two sets of positions.
    
    Args:
        positions1: First set of positions (N x 3)
        positions2: Second set of positions (M x 3)
        
    Returns:
        Distance matrix (N x M)
    """
    # Ensure inputs are numpy arrays
    pos1 = np.asarray(positions1)
    pos2 = np.asarray(positions2)
    
    # Calculate distance matrix using broadcasting
    diff = pos1[:, np.newaxis, :] - pos2[np.newaxis, :, :]
    distances = np.linalg.norm(diff, axis=2)
    
    return distances


def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """
    Safely divide two numbers, returning default if denominator is zero.
    
    Args:
        numerator: Numerator
        denominator: Denominator
        default: Default value if denominator is zero
        
    Returns:
        Result of division or default value
    """
    return numerator / denominator if abs(denominator) > 1e-10 else default


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    """
    Normalize a vector, returning zero vector if norm is zero.
    
    Args:
        vector: Input vector
        
    Returns:
        Normalized vector
    """
    norm = np.linalg.norm(vector)
    return vector / norm if norm > 1e-10 else np.zeros_like(vector)


def rotation_matrix_from_vectors(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    """
    Create rotation matrix to rotate vec1 to vec2.
    
    Args:
        vec1: Source vector
        vec2: Target vector
        
    Returns:
        3x3 rotation matrix
    """
    # Normalize vectors
    a = normalize_vector(vec1)
    b = normalize_vector(vec2)
    
    # Check if vectors are parallel or antiparallel
    dot_product = np.dot(a, b)
    if abs(dot_product) > 0.9999:
        if dot_product > 0:
            return np.eye(3)  # Same direction
        else:
            # Opposite direction - rotate 180 degrees around any perpendicular axis
            perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
            perp = normalize_vector(np.cross(a, perp))
            return 2 * np.outer(perp, perp) - np.eye(3)
    
    # General case using Rodrigues' rotation formula
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = dot_product
    
    vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])
    
    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s * s))
    return R


def compute_center_of_mass(positions: np.ndarray, masses: np.ndarray = None) -> np.ndarray:
    """
    Compute center of mass for a set of positions.
    
    Args:
        positions: Array of positions (N x 3)
        masses: Array of masses (N,). If None, assumes equal masses.
        
    Returns:
        Center of mass position (3,)
    """
    positions = np.asarray(positions)
    
    if masses is None:
        return np.mean(positions, axis=0)
    
    masses = np.asarray(masses)
    if len(masses) != len(positions):
        raise ValueError("Number of masses must match number of positions")
    
    total_mass = np.sum(masses)
    if total_mass == 0:
        return np.mean(positions, axis=0)
    
    return np.sum(positions * masses[:, np.newaxis], axis=0) / total_mass


def format_energy_matrix(matrix: np.ndarray, precision: int = 2) -> str:
    """
    Format an energy matrix for display.
    
    Args:
        matrix: Energy matrix
        precision: Number of decimal places
        
    Returns:
        Formatted string representation
    """
    return np.array2string(matrix, precision=precision, suppress_small=True)
