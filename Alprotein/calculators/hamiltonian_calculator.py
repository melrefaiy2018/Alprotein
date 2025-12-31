"""
Main Hamiltonian calculator that orchestrates site energy and coupling calculations.
"""

import logging
import time
import pandas as pd
import numpy as np
from typing import TYPE_CHECKING, Dict, Any, Optional, Tuple, List
from .site_energy_calculator import SiteEnergyCalculator
from .coupling_calculator import CouplingCalculator
from ..core.constants import (
    DEFAULT_F_VAL,
    DEFAULT_E_VAC,
)

logger = logging.getLogger(__name__)
if not logger.hasHandlers():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")

if TYPE_CHECKING:
    from ..core.pigment_system import PigmentSystem


class HamiltonianCalculator:
    """
    Main calculator for constructing Hamiltonian matrices.
    
    This class orchestrates the site energy and coupling calculations to build
    the complete electronic Hamiltonian for a pigment system.
    """
    
    def __init__(self, pigment_system: 'PigmentSystem', dielectric_cdc: int, dielectric_tresp:int, f_val:int, E_0a:int, E_0b:int):
        """
        Initialize Hamiltonian calculator.

        Args:
            pigment_system: PigmentSystem containing all pigments
            dielectric_cdc: Dielectric constant for CDC calculations
            dielectric_tresp: Dielectric constant for TRESP calculations
            f_val: F value for coupling calculations
            E_0a: Vacuum energy for CLA pigment type
            E_0b: Vacuum energy for CHL pigment type
        """
        self.system = pigment_system
        self.pigment_names = pigment_system.get_pigment_names()
        self.site_energy_calc = None
        self.coupling_calc = None
        self.hamiltonian = None

        # Define default dielectric constants for CDC and TRESP
        self.dielectric_cdc = dielectric_cdc
        self.dielectric_tresp = dielectric_tresp
        self.f_val = f_val 
        self.e0_a = E_0a
        self.e0_b = E_0b

    def construct_hamiltonian(self, params: Optional[Dict[str, Any]] = None) -> pd.DataFrame:
        """
        Calculate and assemble the full Hamiltonian matrix.
        
        Args:
            params: Dictionary with calculation parameters
            
        Returns:
            Hamiltonian matrix as pandas DataFrame
        """
        # Set default parameters if not provided
        if params is None:
            params = {}
        
        # Extract parameters with defaults
        dielectric_cdc = self.dielectric_cdc
        dielectric_tresp = self.dielectric_tresp
        f_val = self.f_val
        coupling_calc = params.get('coupling_calc', 'tresp')
        
        # Initialize calculators
        self.site_energy_calc = SiteEnergyCalculator(dielectric_constant=dielectric_cdc, e0a=self.e0_a, e0b=self.e0_b)
        self.coupling_calc = CouplingCalculator(dielectric=dielectric_tresp, f_val=f_val)
        
        logger.info("--- Starting Hamiltonian Construction ---")
        start_time = time.time()
        
        # Initialize Hamiltonian matrix
        hamiltonian = pd.DataFrame(
            index=self.pigment_names, 
            columns=self.pigment_names, 
            dtype=float
        )
        hamiltonian.fillna(0.0, inplace=True)
        
        # Calculate site energies (diagonal elements)
        logger.info("Calculating site energies...")
        # Create vacuum energies dictionary for fallback
        vacuum_energies = {'CLA': self.e0_a, 'CHL': self.e0_b}
        site_energies = self.site_energy_calc.calculate_all_site_energies(
            self.system, vacuum_energies
        )
        
        # Fill diagonal elements
        for name, energy in site_energies.items():
            hamiltonian.loc[name, name] = energy
        
        # Calculate couplings (off-diagonal elements)
        logger.info("Calculating couplings...")
        coupling_count = 0
        total_pairs = len(self.pigment_names) * (len(self.pigment_names) - 1) // 2
        
        for i, pigA_name in enumerate(self.pigment_names):
            for j, pigB_name in enumerate(self.pigment_names):
                if i >= j:
                    continue

                try:
                    pigA = self.system[pigA_name]
                    pigB = self.system[pigB_name]

                    coupling = self.coupling_calc.calculate_coupling(
                        pigA, pigB, method=coupling_calc
                    )

                    # Fill both upper and lower triangular elements
                    hamiltonian.loc[pigA_name, pigB_name] = coupling
                    hamiltonian.loc[pigB_name, pigA_name] = coupling

                    coupling_count += 1
                    if coupling_count % 10 == 0 or coupling_count == total_pairs:
                        logger.info(
                            "  Calculated %d/%d couplings",
                            coupling_count,
                            total_pairs,
                        )
                except Exception as e:
                    error_msg = (
                        f"Failed to calculate coupling between {pigA_name} and {pigB_name} "
                        f"(indices {i}, {j}):\n"
                        f"Error: {str(e)}\n"
                        f"Pigment A: {pigA_name}\n"
                        f"Pigment B: {pigB_name}"
                    )
                    logger.error(error_msg)
                    raise RuntimeError(error_msg) from e
        
        end_time = time.time()
        logger.info(
            "--- Hamiltonian Construction Complete (%.2fs) ---",
            end_time - start_time,
        )
        self.hamiltonian = hamiltonian
        return hamiltonian

    def build_domains(
        self,
        cutoff: float,
        hamiltonian: Optional[Any] = None,
        use_absolute: bool = True,
    ) -> Dict[int, List[int]]:
        """
        Build domains by clustering pigments with couplings above a cutoff.

        Args:
            cutoff: Coupling cutoff threshold (cm-1)
            hamiltonian: Hamiltonian matrix (DataFrame or ndarray). If None,
                use the last constructed Hamiltonian.
            use_absolute: If True, compare |coupling| to cutoff.

        Returns:
            dict: {domain_id: [pigment_indices]}
        """
        if hamiltonian is None:
            if self.hamiltonian is None:
                raise ValueError("Hamiltonian not provided and no cached Hamiltonian is available.")
            hamiltonian = self.hamiltonian

        if hasattr(hamiltonian, "values"):
            matrix = hamiltonian.values
        else:
            matrix = np.asarray(hamiltonian)

        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError("Hamiltonian must be a square matrix.")

        n = matrix.shape[0]
        graph = [set() for _ in range(n)]

        for i in range(n):
            for j in range(i + 1, n):
                value = matrix[i, j]
                if np.isnan(value):
                    continue
                compare_value = abs(value) if use_absolute else value
                if compare_value >= cutoff:
                    graph[i].add(j)
                    graph[j].add(i)

        visited = [False] * n
        domains = []
        for node in range(n):
            if visited[node]:
                continue
            stack = [node]
            domain = []
            while stack:
                current = stack.pop()
                if visited[current]:
                    continue
                visited[current] = True
                domain.append(current)
                stack.extend(graph[current])
            domains.append(sorted(domain))

        return {domain_id: domain for domain_id, domain in enumerate(domains)}
    

    def __repr__(self) -> str:
        """Return a concise representation of the calculator."""
        return f"HamiltonianCalculator(system with {len(self.pigment_names)} pigments)"
