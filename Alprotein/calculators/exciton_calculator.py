"""
Exciton distribution calculator with ensemble averaging.

Provides the ExcitonCalculator class for computing exciton state probability
distributions across pigments using disorder ensemble averaging.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any


class ExcitonCalculator:
    """
    Calculator for exciton state distributions with ensemble averaging.

    Calculates the probability distribution of exciton states across pigments
    using disorder ensemble averaging, similar to SpectraCalculator. This helps
    understand exciton localization and energy distribution across the system.

    The calculator delegates to proven existing implementations for core
    calculations while providing a clean interface consistent with other
    AlProtein calculators.
    """

    def __init__(self, disorder_sigma: float, n_ensemble: int, temperature: float = 300.0):
        """
        Initialize exciton distribution calculator.

        Args:
            disorder_sigma: Standard deviation of Gaussian disorder (cm⁻¹).
                          This should match the sigma from SpectraCalculator
                          (already in σ form, not FWHM).
            n_ensemble: Number of disorder realizations for ensemble averaging.
            temperature: Temperature in Kelvin (for future use, currently not
                        used in absorption calculations).
        """
        self.disorder_sigma = disorder_sigma
        self.n_ensemble = n_ensemble
        self.temperature = temperature

        # Storage for calculation results
        self.distributions = None
        self.site_labels = None

    def calculate_distributions(
        self,
        hamiltonian: pd.DataFrame,
        domains: Dict[int, List[int]],
        pigment_system: Any
    ) -> Tuple[Dict[str, Tuple[List[float], List[float]]], List[str]]:
        """
        Calculate exciton distributions for all pigments.

        This method applies ensemble averaging with Gaussian disorder to compute
        the probability distribution of exciton states for each pigment site.

        Args:
            hamiltonian: pd.DataFrame with Hamiltonian matrix. Rows and columns
                        are indexed by pigment names. Values are in cm⁻¹.
                        (from HamiltonianCalculator.construct_hamiltonian())
            domains: Dict[int, List[int]] mapping domain IDs to lists of pigment
                    indices. (from HamiltonianCalculator.build_domains())
            pigment_system: PigmentSystem instance containing pigment data.

        Returns:
            Tuple containing:
                - distributions: Dict mapping pigment names to tuples of
                                (energies, probabilities). Energies are in cm⁻¹.
                - site_labels: List of pigment names in order from Hamiltonian.

        Example:
            >>> exciton_calc = ExcitonCalculator(disorder_sigma=55.0, n_ensemble=300)
            >>> distributions, labels = exciton_calc.calculate_distributions(
            ...     hamiltonian=hamiltonian_df,
            ...     domains=domains_dict,
            ...     pigment_system=pigment_system
            ... )
            >>> # Access distribution for a specific pigment
            >>> energies, probs = distributions['A_CLA_501']
        """
        # Convert domains from index-based to name-based format
        pigment_names = hamiltonian.columns.tolist()
        list_pigment_domains = self._convert_domains_to_names(domains, pigment_names)

        # Import and delegate to existing proven implementation
        from ..utils.calculate_exciton_distribution import calculate_exciton_distribution

        distributions, site_labels = calculate_exciton_distribution(
            H_data=hamiltonian,
            list_pigment_domains=list_pigment_domains,
            N_ens=self.n_ensemble,
            sigma_e=self.disorder_sigma,
            protein_atomic=pigment_system,  # PigmentSystem has .dict_pigments attribute
            calc_type='abs'  # Absorption mode only
        )

        # Store results for later use in plotting/export
        self.distributions = distributions
        self.site_labels = site_labels

        return distributions, site_labels

    def _convert_domains_to_names(
        self,
        domains: Dict[int, List[int]],
        pigment_names: List[str]
    ) -> List[List[str]]:
        """
        Convert domain format from indices to names.

        HamiltonianCalculator.build_domains() returns Dict[int, List[int]]
        where values are pigment indices. The underlying calculation function
        expects List[List[str]] where values are pigment names. This method
        performs the conversion.

        Args:
            domains: Dict mapping domain IDs to lists of pigment indices.
                    E.g., {0: [0, 1, 2], 1: [3, 4, 5]}
            pigment_names: Ordered list of pigment names from Hamiltonian columns.

        Returns:
            List of domains where each domain is a list of pigment names.
            E.g., [['A_CLA_501', 'A_CLA_502'], ['A_CLA_503']]
        """
        list_pigment_domains = []

        # Sort domain IDs for deterministic ordering
        for domain_id in sorted(domains.keys()):
            indices = domains[domain_id]

            # Skip empty domains
            if not indices:
                continue

            # Convert indices to names
            domain_names = [pigment_names[idx] for idx in indices]
            list_pigment_domains.append(domain_names)

        return list_pigment_domains

    def plot_distributions(
        self,
        output_path: str,
        show_labels: bool = True,
        x_min: int = 600,
        x_max: int = 720,
        x_interval: int = 20,
        legend: bool = False,
        **kwargs
    ) -> None:
        """
        Plot exciton distributions for all pigments.

        Creates a plot showing the kernel density estimate (KDE) of exciton
        probability distributions for each pigment. Each pigment is shown in
        a different color with optional peak annotations.

        Args:
            output_path: Directory path where plot will be saved.
                        Plot filename will be 'exciton_distribution.png'.
            show_labels: If True, annotate the highest peak for each pigment
                        with its residue number.
            x_min: Minimum wavelength for x-axis (nm). Default 600.
            x_max: Maximum wavelength for x-axis (nm). Default 720.
            x_interval: Interval for x-axis ticks (nm). Default 20.
            legend: If True, show legend with all pigment labels.
            **kwargs: Additional plotting parameters:
                     - fontsize: Font size for annotations (default 16)
                     - fontweight: Font weight ('bold', 'normal', etc.)
                     - font_style: Font family ('serif', 'sans-serif', etc.)

        Raises:
            ValueError: If calculate_distributions() has not been called yet.

        Example:
            >>> exciton_calc.plot_distributions(
            ...     output_path='/path/to/results/',
            ...     show_labels=True,
            ...     fontsize=16,
            ...     fontweight='bold'
            ... )
        """
        if self.distributions is None:
            raise ValueError(
                "No distributions calculated. "
                "Run calculate_distributions() first."
            )

        from ..utils.calculate_exciton_distribution import plot_exciton_distribution

        plot_exciton_distribution(
            exciton_distribution=self.distributions,
            list_site_label=self.site_labels,
            path_saving=output_path,
            show_labels=show_labels,
            x_min=x_min,
            x_max=x_max,
            x_interval=x_interval,
            legend=legend,
            **kwargs
        )

    def plot_combined_with_absorption(
        self,
        wavelengths_abs: np.ndarray,
        absorption: np.ndarray,
        exp_absorption: Optional[np.ndarray] = None,
        output_path: str = '.',
        temp: int = 300,
        show_labels: bool = True,
        x_min: Optional[int] = None,
        x_max: Optional[int] = None,
        **kwargs
    ) -> None:
        """
        Plot combined absorption spectrum and exciton distributions.

        Creates a two-panel plot with absorption spectrum on top and exciton
        distributions on bottom, useful for correlating spectral features with
        exciton localization patterns.

        Args:
            wavelengths_abs: Wavelength axis from SpectraCalculator (nm).
            absorption: Absorption spectrum from SpectraCalculator (normalized).
            exp_absorption: Optional experimental absorption data. Should be a
                           2-column array with wavelengths in column 0 and
                           absorption in column 1.
            output_path: Directory path for saving plot.
                        Filename will be 'exciton_distribution_with_absorption.png'.
            temp: Temperature in Kelvin (used for plot labeling).
            show_labels: If True, annotate exciton distribution peaks.
            x_min: Minimum wavelength (nm). Default 600.
            x_max: Maximum wavelength (nm). Default 720.
            **kwargs: Additional plotting parameters (fontsize, fontweight, font_style).

        Raises:
            ValueError: If calculate_distributions() has not been called yet.

        Example:
            >>> exciton_calc.plot_combined_with_absorption(
            ...     wavelengths_abs=wavelengths,
            ...     absorption=abs_spectrum,
            ...     exp_absorption=exp_data,
            ...     output_path='/path/to/results/',
            ...     temp=300
            ... )
        """
        if self.distributions is None:
            raise ValueError(
                "No distributions calculated. "
                "Run calculate_distributions() first."
            )

        from ..utils.calculate_exciton_distribution import plot_combined_absorption_and_exciton

        plot_combined_absorption_and_exciton(
            sliced_lambda_axis_a=wavelengths_abs,
            absorption_data=absorption,
            exp_absorption=exp_absorption,
            exciton_distribution=self.distributions,
            list_site_label=self.site_labels,
            temp=temp,
            path_saving=output_path,
            show_labels=show_labels,
            x_min=x_min or 600,
            x_max=x_max or 720,
            **kwargs
        )

    def export_distributions(self, output_path: str) -> None:
        """
        Export exciton distributions to CSV file.

        Creates a CSV file with columns: pigment, energy_cm-1, wavelength_nm, probability.
        Each row represents one (energy, probability) pair for a pigment.

        Args:
            output_path: Full path to output CSV file (including filename).

        Raises:
            ValueError: If calculate_distributions() has not been called yet.

        Example:
            >>> exciton_calc.export_distributions('/path/to/exciton_distributions.csv')
        """
        if self.distributions is None:
            raise ValueError(
                "No distributions calculated. "
                "Run calculate_distributions() first."
            )

        # Build rows for CSV export
        rows = []
        for site in self.site_labels:
            energies, probabilities = self.distributions[site]

            for energy, prob in zip(energies, probabilities):
                # Convert energy to wavelength, handling potential division by zero
                wavelength = 1e7 / energy if energy > 0 else np.nan

                rows.append({
                    'pigment': site,
                    'energy_cm-1': energy,
                    'wavelength_nm': wavelength,
                    'probability': prob
                })

        # Create DataFrame and save
        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)

    def __repr__(self) -> str:
        """Return a concise representation of the calculator."""
        return (
            f"ExcitonCalculator("
            f"disorder_sigma={self.disorder_sigma:.1f} cm⁻¹, "
            f"n_ensemble={self.n_ensemble}, "
            f"temperature={self.temperature:.1f} K)"
        )
