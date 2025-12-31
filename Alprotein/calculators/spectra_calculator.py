"""
Absorption and fluorescence spectra using Renger lineshape theory.

Provides the SpectraCalculator class for computing spectra with a
double-overdamped Brownian oscillator spectral density, vibronic 0-0 and
0-1 transitions, and ensemble averaging under Gaussian disorder.

Units: energies in cm^-1, time in fs, wavelength in nm, dipoles in Debye.
"""

import numpy as np
from typing import Dict, Tuple, Optional
from scipy.constants import hbar, c, k


class SpectraCalculator:
    """
    Advanced absorption spectra calculator using Renger lineshape theory

    References:
    - Renger et al., J. Phys. Chem. B 1997, 101, 7232-7242
    - Renger & Marcus, J. Chem. Phys. 2002, 116, 9997-10019
    """

    def __init__(self,
                 temperature,
                 disorder_sigma,
                 n_ensemble,
                 dt,
                 t_max):
        """
        Initialize the spectra calculator and precompute core functions.

        Args:
            temperature: Temperature in Kelvin.
            disorder_sigma: FWHM of Gaussian disorder in cm^-1.
            n_ensemble: Number of ensemble realizations.
            dt: Time axis step in fs.
            t_max: Maximum time in fs.
        """
        self.temperature = temperature
        self.disorder_sigma = disorder_sigma / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to σ
        self.n_ensemble = n_ensemble
        self.dt = dt
        self.t_max = t_max

        # Physical constants (CGS units for consistency)
        self.hbar_cgs = 1.05457e-27  # erg⋅s
        self.c_cgs = 2.998e10  # cm/s
        self.kB = 0.695  # cm⁻¹/K (Boltzmann constant)

        # Vibrational parameters from experimental data
        # Frequencies in cm⁻¹, Huang-Rhys factors dimensionless
        self.vib_frequencies = np.array([100, 175, 250, 300, 375, 500, 600, 725, 800, 875])
        self.huang_rhys_factors = np.array([0.2, 0.1, 0.06, 0.04, 0.06, 0.04, 0.015, 0.04, 0.02, 0.02])

        # Spectral density parameters (Renger double-overdamped Brownian oscillator)
        self.S0 = 0.5        # Reorganization energy prefactor
        self.s1 = 0.8        # First oscillator weight
        self.s2 = 0.5        # Second oscillator weight
        self.w1 = 0.069 * 8.0656  # First oscillator frequency (meV → cm⁻¹)
        self.w2 = 0.24 * 8.0656   # Second oscillator frequency (meV → cm⁻¹)

        # Setup calculation axes
        self.setup_axes()

        # Precalculate spectral density and line broadening function
        self.J_w = self.calculate_spectral_density()
        self.reorganization_energy = np.trapz(self.w_axis * self.J_w, x=self.w_axis)
        self.g_t = self.calculate_g_function()

    def setup_axes(self):
        """
        Initialize time and frequency axes used by the calculator.

        Sets:
            self.t_axis: Time axis in fs.
            self.w_axis: Frequency axis in cm^-1.
            self.dw: Frequency step in cm^-1.
            self.w_max: Maximum frequency in cm^-1.
        """
        # Time axis for Fourier transforms
        self.t_axis = np.arange(0, self.t_max, self.dt)  # fs

        # Frequency axis for spectral density
        self.dw = 0.1  # cm⁻¹
        self.w_max = 2000  # cm⁻¹
        self.w_axis = np.arange(self.dw, self.w_max, self.dw)

    def calculate_spectral_density(self):
        """
        Calculate spectral density J(w) using Renger's functional form.

        Uses the exponential decay model from Renger et al.:
        J(w) = (S0/(s1+s2)) * w^3/10080 * [s1/w1^4 exp(-sqrt(w/w1)) +
        s2/w2^4 exp(-sqrt(w/w2))]

        Returns:
            Spectral density array in cm^-1.
        """
        w = self.w_axis

        # Prefactor
        prefactor = self.S0 / (self.s1 + self.s2) * np.power(w, 3) / 10080

        # First oscillator contribution with exponential decay
        sd_1 = self.s1 / self.w1**4 * np.exp(-np.sqrt(w / self.w1))

        # Second oscillator contribution
        sd_2 = self.s2 / self.w2**4 * np.exp(-np.sqrt(w / self.w2))

        # Combined spectral density
        J_w = prefactor * (sd_1 + sd_2)

        return J_w

    def bose_einstein(self, w):
        """
        Compute the Bose-Einstein distribution n(w, T).

        Args:
            w: Frequency array in cm^-1.

        Returns:
            Bose-Einstein occupation number array.
        """
        if self.temperature == 0:
            return np.zeros_like(w)

        # Avoid division by zero
        x = w / (self.kB * self.temperature)
        x = np.clip(x, 1e-10, 100)  # Prevent overflow

        n = 1.0 / (np.exp(x) - 1.0)
        return n

    def calculate_g_function(self):
        """
        Calculate the phonon-exciton coupling function g(t).

        Implements Renger's formulation:
        g(t) = integral d(w) J(w) [(1+n(w)) exp(-i w t) + n(w) exp(i w t)]

        Returns:
            Complex line broadening function array.
        """
        t = self.t_axis
        w = self.w_axis
        J = self.J_w

        # Bose-Einstein distribution
        n_w = self.bose_einstein(w)

        # Calculate g(t) using Renger's formulation
        # Note: hbar is in cgs units (erg·s), w is in cm⁻¹
        # Need to convert: ω(cm⁻¹) → ω(rad/s) = ω * 2πc
        # Then ω*t where t is in fs → needs factor of 1e-15

        g_t = np.array([
            np.trapz(
                (1 + n_w) * J * np.exp(-1j * w * 2 * np.pi * self.c_cgs * ti * 1e-15)
                + n_w * J * np.exp(1j * w * 2 * np.pi * self.c_cgs * ti * 1e-15),
                x=w
            ) for ti in t
        ])

        return g_t

    def calculate_domain_exciton_properties(self, hamiltonian, domains, pigment_system,
                                           list_site_label=None, dict_dipole_by_site=None):
        """
        Calculate exciton properties for each domain.

        Args:
            hamiltonian: Full Hamiltonian matrix in cm^-1.
            domains: Mapping {domain_id: [pigment_indices]}.
            pigment_system: PigmentSystem instance with pigment data.
            list_site_label: Optional list of site labels.
            dict_dipole_by_site: Optional mapping of dipole magnitudes by site.

        Returns:
            Tuple containing:
                - dict_domain_coefficients: {domain_id: eigenvectors}
                - dict_domain_energies: {domain_id: eigenvalues}
                - dict_domain_dipoles: {domain_id: transition_dipole_vectors}
                - dict_domain_lifetimes: {domain_id: lifetime_broadenings}
                - dict_exciton_gamma: {domain_id: gamma_factors}
        """
        dict_coefficients = {}
        dict_energies = {}
        dict_dipoles = {}
        dict_lifetimes = {}
        dict_exciton_gamma = {}

        # Get site labels if not provided
        if list_site_label is None:
            list_site_label = list(pigment_system.pigments.keys())

        # Get dipole magnitudes if not provided
        if dict_dipole_by_site is None:
            dict_dipole_by_site = {}
            for name, pigment in pigment_system.pigments.items():
                resname = pigment.get_resname()
                # Default dipole moments for CLA and CHL
                if 'CHL' in resname:
                    dict_dipole_by_site[name] = 4.61  # Debye
                else:  # CLA
                    dict_dipole_by_site[name] = 5.47  # Debye

        # Build domain index mapping
        dict_domain_index = {}
        dict_domain_labels = {}
        for domain_id, pig_indices in domains.items():
            dict_domain_index[domain_id] = pig_indices
            dict_domain_labels[domain_id] = [list_site_label[idx] for idx in pig_indices]

        for domain_id, pig_indices in domains.items():
            if len(pig_indices) == 0:
                continue

            # Extract domain submatrix
            H_domain = hamiltonian[np.ix_(pig_indices, pig_indices)]

            # Diagonalize domain Hamiltonian
            eigenvalues, eigenvectors = np.linalg.eigh(H_domain)

            dict_coefficients[domain_id] = eigenvectors
            dict_energies[domain_id] = eigenvalues

            # Calculate gamma factors (Σ|c_n|^4 for each exciton)
            gamma_factors = np.sum(np.abs(eigenvectors)**4, axis=0)
            dict_exciton_gamma[domain_id] = gamma_factors

            # Calculate lifetime broadening using Renger's approach
            # This includes energy transfer between exciton states
            list_tau = []
            for exciton in range(len(eigenvalues)):
                tau_exciton = 0
                for exciton_N in range(len(eigenvalues)):
                    if exciton_N != exciton:
                        omega_MN = eigenvalues[exciton] - eigenvalues[exciton_N]
                        gamma_MN = np.dot(
                            np.abs(eigenvectors[:, exciton])**2,
                            np.abs(eigenvectors[:, exciton_N])**2
                        )

                        # Calculate J(ω) at this frequency
                        if omega_MN > 0:
                            # Find closest frequency in w_axis
                            idx = np.argmin(np.abs(self.w_axis - omega_MN))
                            J_omega = self.J_w[idx]
                            n_omega = self.bose_einstein(np.array([omega_MN]))[0]

                            tau_exciton += np.pi * gamma_MN * omega_MN**2 * (1 + n_omega) * J_omega
                        else:
                            omega_abs = np.abs(omega_MN)
                            idx = np.argmin(np.abs(self.w_axis - omega_abs))
                            J_omega = self.J_w[idx]
                            n_omega = self.bose_einstein(np.array([omega_abs]))[0]

                            tau_exciton += np.pi * gamma_MN * omega_MN**2 * n_omega * J_omega

                list_tau.append(tau_exciton)

            dict_lifetimes[domain_id] = list_tau

            # Calculate exciton dipole moments using actual dipole vectors
            list_dipole_by_exciton = []
            for exciton in range(len(eigenvalues)):
                mu_exc = np.zeros(3)
                for (index_site, label_site) in enumerate(dict_domain_labels[domain_id]):
                    try:
                        dipole_dir = pigment_system.pigments[label_site].get_transition_dipole_moment()
                        if np.linalg.norm(dipole_dir) > 0:
                            dipole_dir = dipole_dir / np.linalg.norm(dipole_dir)
                    except:
                        # Fallback to z-direction
                        dipole_dir = np.array([0., 0., 1.])

                    mu_exc += (eigenvectors[index_site, exciton] * dipole_dir
                              * dict_dipole_by_site.get(label_site, 5.0))

                list_dipole_by_exciton.append(mu_exc)

            dict_dipoles[domain_id] = list_dipole_by_exciton

        return dict_coefficients, dict_energies, dict_dipoles, dict_lifetimes, dict_exciton_gamma

    def calculate_spectra_renger(self, hamiltonian, domains, pigment_system,
                                list_site_label=None, dict_dipole_by_site=None):
        """
        Calculate absorption and fluorescence spectra using Renger's formulation.

        This follows the approach used in the CP29_Doran reference implementation.

        Args:
            hamiltonian: Hamiltonian matrix in cm^-1.
            domains: Mapping {domain_id: [pigment_indices]}.
            pigment_system: PigmentSystem instance.
            list_site_label: Optional list of site labels.
            dict_dipole_by_site: Optional mapping of dipole magnitudes by site.

        Returns:
            Tuple:
                F_w: Fluorescence spectrum array.
                A_w: Absorption spectrum array.
                Fw_site: Site-level fluorescence spectrum.
                dict_exciton_fl_contribution: Per-exciton fluorescence contributions.
                Aw_site: Site-level absorption spectrum.
                dict_exciton_a_contribution: Per-exciton absorption contributions.
        """
        # Get site labels if not provided
        if list_site_label is None:
            list_site_label = list(pigment_system.pigments.keys())

        # Get dipole magnitudes if not provided
        if dict_dipole_by_site is None:
            dict_dipole_by_site = {}
            for name in list_site_label:
                pigment = pigment_system.pigments[name]
                resname = pigment.get_resname()
                if 'CHL' in resname:
                    dict_dipole_by_site[name] = 4.61  # Debye
                else:  # CLA
                    dict_dipole_by_site[name] = 5.47  # Debye

        # Build domain labels
        dict_domain_labels = {}
        for domain_id, pig_indices in domains.items():
            dict_domain_labels[domain_id] = [list_site_label[idx] for idx in pig_indices]

        # Calculate central frequency
        site_energies = np.diag(hamiltonian)
        central_freq = np.mean(site_energies)

        # Get domain exciton properties
        dict_domain_coefficients, dict_domain_energies, dict_domain_dipoles, \
        dict_domain_lifetimes, dict_exciton_gamma = \
            self.calculate_domain_exciton_properties(hamiltonian, domains, pigment_system,
                                                     list_site_label, dict_dipole_by_site)

        # Calculate exciton transition energies with reorganization energy
        dict_exciton_etrans = {}
        for domain_id in domains.keys():
            if domain_id not in dict_domain_energies:
                continue
            dict_exciton_etrans[domain_id] = [
                dict_domain_energies[domain_id][exciton] -
                (dict_exciton_gamma[domain_id][exciton] * self.reorganization_energy)
                for exciton in range(len(dict_exciton_gamma[domain_id]))
            ]

        # Calculate Boltzmann thermalization
        kB_T = self.kB * self.temperature
        bolt_denom = 0.
        for domain_id in domains.keys():
            if domain_id not in dict_exciton_etrans:
                continue
            for exciton in range(len(dict_exciton_gamma[domain_id])):
                bolt_denom += np.exp((-dict_exciton_etrans[domain_id][exciton]) / kB_T)

        dict_exciton_therm = {}
        for domain_id in domains.keys():
            if domain_id not in dict_exciton_etrans:
                continue
            dict_exciton_therm[domain_id] = [
                np.exp((-dict_exciton_etrans[domain_id][exciton]) / kB_T) / bolt_denom
                for exciton in range(len(dict_exciton_gamma[domain_id]))
            ]

        # === Construct vibronic time-domain functions ===

        # For fluorescence (d_ti_pre)
        d_ti_pre = np.zeros(len(self.t_axis), dtype=np.complex128)
        for (w_vib, s_vib) in zip(self.vib_frequencies, self.huang_rhys_factors):
            fc_01_sq = np.exp(-s_vib) * s_vib
            fc_00_sq = np.exp(-s_vib)
            w_mi = -w_vib - self.reorganization_energy - central_freq

            # Convert to angular frequency
            omega_mi = w_mi * 2 * np.pi * self.c_cgs * 1e-15  # rad/fs

            g_mi = self.g_t - self.g_t[0]
            d_ti_pre += np.exp(1j * omega_mi * self.t_axis) * np.exp(g_mi) * fc_01_sq / fc_00_sq

        # For absorption (Dt_pre)
        Dt_pre = np.zeros(len(self.t_axis), dtype=np.complex128)
        for (w_vib, s_vib) in zip(self.vib_frequencies, self.huang_rhys_factors):
            fc_01_sq = np.exp(-s_vib) * s_vib
            fc_00_sq = np.exp(-s_vib)
            w_mi = w_vib - self.reorganization_energy - central_freq

            omega_mi = w_mi * 2 * np.pi * self.c_cgs * 1e-15  # rad/fs

            Dt_pre += fc_01_sq / fc_00_sq * np.exp(-1j * omega_mi * self.t_axis) * np.exp(self.g_t - self.g_t[0])

        # === Calculate exciton contributions ===
        dict_exciton_fl_contribution = {}
        dict_exciton_a_contribution = {}
        Fw_site = np.zeros(len(self.t_axis), dtype=np.float64)

        for domain_id in domains.keys():
            if domain_id not in dict_domain_energies:
                continue

            list_fl_contribution = []
            list_a_contribution = []

            for exciton in range(len(dict_domain_energies[domain_id])):
                w_Md = dict_exciton_etrans[domain_id][exciton] - central_freq
                omega_Md = w_Md * 2 * np.pi * self.c_cgs * 1e-15  # rad/fs

                G_Md = dict_exciton_gamma[domain_id][exciton] * self.g_t

                # Absorption 0-0
                Dt_a = (np.exp(-1j * omega_Md * self.t_axis) *
                       np.exp(G_Md - G_Md[0]) *
                       np.exp(-self.t_axis * dict_domain_lifetimes[domain_id][exciton] *
                             2 * np.pi * self.c_cgs * 1e-15))
                Dw_a = 1 / (2 * np.pi) * np.fft.fft(Dt_a)

                mu_ex = dict_domain_dipoles[domain_id][exciton]
                mu_ex_sq = np.linalg.norm(mu_ex)**2

                A_cont = mu_ex_sq * np.real(Dw_a)

                # Fluorescence 0-0
                Dt_f = (np.exp(1j * omega_Md * self.t_axis) *
                       np.exp(G_Md - G_Md[0]) *
                       np.exp(-self.t_axis * dict_domain_lifetimes[domain_id][exciton] *
                             2 * np.pi * self.c_cgs * 1e-15))
                Dw_f = 1 / (2 * np.pi) * np.fft.fft(Dt_f)

                F_cont = dict_exciton_therm[domain_id][exciton] * mu_ex_sq * np.real(Dw_f)

                # Fluorescence 0-1 site contributions
                Ft_cont = np.zeros(len(self.t_axis), dtype=np.complex128)
                dict_domain_index = {d: list(inds) for d, inds in domains.items()}

                for site_index in dict_domain_index[domain_id]:
                    site_label = list_site_label[site_index]
                    coef_index = np.where(np.array(dict_domain_index[domain_id]) == site_index)[0][0]
                    coef_sq = dict_domain_coefficients[domain_id][:, exciton][coef_index]**2

                    omega_site = site_energies[site_index] * 2 * np.pi * self.c_cgs * 1e-15

                    d_ti_site = (dict_exciton_therm[domain_id][exciton] *
                                dict_dipole_by_site[site_label]**2 *
                                coef_sq *
                                np.exp(1j * omega_site * self.t_axis) *
                                d_ti_pre)
                    Ft_cont += d_ti_site

                Fw_site += np.real(np.fft.fft(Ft_cont) / (2 * np.pi))

                list_fl_contribution.append(F_cont)
                list_a_contribution.append(A_cont)

            dict_exciton_fl_contribution[domain_id] = list_fl_contribution
            dict_exciton_a_contribution[domain_id] = list_a_contribution

        # === Calculate site-level absorption 0-1 ===
        Dt_site = np.zeros(len(self.t_axis), dtype=np.complex128)
        for site_index in range(len(list_site_label)):
            site_label = list_site_label[site_index]
            omega_site = site_energies[site_index] * 2 * np.pi * self.c_cgs * 1e-15

            Dt_site += (dict_dipole_by_site[site_label]**2 *
                       np.exp(-1j * omega_site * self.t_axis) *
                       Dt_pre)

        Aw_site = np.real(np.fft.fft(Dt_site) / (2 * np.pi))

        # === Sum up total fluorescence ===
        F_w = Fw_site.copy()
        for domain_id in domains.keys():
            if domain_id not in dict_exciton_fl_contribution:
                continue
            for exciton in range(len(dict_exciton_fl_contribution[domain_id])):
                F_w += dict_exciton_fl_contribution[domain_id][exciton]

        # === Sum up total absorption ===
        A_w = Aw_site.copy()
        for domain_id in domains.keys():
            if domain_id not in dict_exciton_a_contribution:
                continue
            for exciton in range(len(dict_exciton_a_contribution[domain_id])):
                A_w += dict_exciton_a_contribution[domain_id][exciton]

        return F_w, A_w, Fw_site, dict_exciton_fl_contribution, Aw_site, dict_exciton_a_contribution

    def calculate_spectrum(self, hamiltonian, domains, pigment_system,
                          include_vibronic=True, use_ensemble=True,
                          list_site_label=None, dict_dipole_by_site=None,
                          progress_callback=None):
        """
        Calculate absorption and fluorescence spectra with optional disorder.

        Args:
            hamiltonian: NxN Hamiltonian matrix in cm^-1.
            domains: Domain clustering mapping.
            pigment_system: PigmentSystem instance.
            include_vibronic: Include 0-1 vibronic transitions (Renger method).
            use_ensemble: Enable ensemble averaging with disorder.
            list_site_label: Optional list of site labels.
            dict_dipole_by_site: Optional mapping of dipole magnitudes.
            progress_callback: Optional callback for progress updates.

        Returns:
            Tuple containing:
                wavelengths_abs: Wavelength axis for absorption in nm.
                absorption: Normalized absorption spectrum.
                wavelengths_fl: Wavelength axis for fluorescence in nm.
                fluorescence: Normalized fluorescence spectrum.
        """
        # Convert to numpy array if it's a DataFrame
        if hasattr(hamiltonian, 'values'):
            hamiltonian = hamiltonian.values
        else:
            hamiltonian = np.asarray(hamiltonian)

        # Get central frequency for wavelength conversion
        central_freq = np.mean(np.diag(hamiltonian))

        if use_ensemble:
            # Ensemble averaging with Gaussian disorder
            F_w_total = np.zeros(len(self.t_axis))
            A_w_total = np.zeros(len(self.t_axis))

            for n in range(self.n_ensemble):
                if progress_callback:
                    progress = 10 + int((n / max(self.n_ensemble, 1)) * 80)
                    progress_callback(f"Spectra ensemble {n + 1}/{self.n_ensemble}", progress)
                # Add diagonal disorder to site energies
                rng = np.random.RandomState(seed=n)
                site_energies = np.diag(hamiltonian)
                disorder = rng.normal(0, self.disorder_sigma, size=len(site_energies))

                H_disorder = hamiltonian.copy()
                np.fill_diagonal(H_disorder, site_energies + disorder)

                # Calculate spectra using Renger method
                F_w, A_w, _, _, _, _ = self.calculate_spectra_renger(
                    H_disorder, domains, pigment_system,
                    list_site_label, dict_dipole_by_site
                )

                F_w_total += F_w
                A_w_total += A_w

            # Average over ensemble
            F_w_total /= self.n_ensemble
            A_w_total /= self.n_ensemble
            if progress_callback:
                progress_callback("Finalizing spectra average...", 90)

        else:
            # Single realization without disorder
            F_w_total, A_w_total, _, _, _, _ = self.calculate_spectra_renger(
                hamiltonian, domains, pigment_system,
                list_site_label, dict_dipole_by_site
            )
            if progress_callback:
                progress_callback("Finalizing spectra...", 90)

        # Convert frequency domain to wavelength axis
        # Following the standard conversion: λ(nm) = 10^7 / ν(cm^-1)

        # Get FFT frequency axis
        freq_fft = np.fft.fftfreq(len(self.t_axis), self.dt)  # Frequency in 1/fs

        # Convert to wavenumber (cm^-1)
        # freq_fft is in 1/fs, need to convert to cm^-1
        # ν(cm^-1) = freq(Hz) / c(cm/s) = freq(1/s) / c(cm/s)
        # Since freq is in 1/fs, multiply by 1e15 to get 1/s
        w_fft = freq_fft * 1e15 / self.c_cgs  # Now in cm^-1

        # Add central frequency to get absolute wavenumbers
        # For absorption (negative frequency convention)
        w_axis_a = central_freq - w_fft

        # For fluorescence (positive frequency convention)
        w_axis_f = central_freq + w_fft

        # Convert wavenumber to wavelength: λ(nm) = 10^7 / ν(cm^-1)
        with np.errstate(divide='ignore', invalid='ignore'):
            lambda_axis_a = np.where(w_axis_a > 0, 1e7 / w_axis_a, np.inf)
            lambda_axis_f = np.where(w_axis_f > 0, 1e7 / w_axis_f, np.inf)

        # Normalize spectra
        if np.max(F_w_total) > 0:
            F_w_norm = F_w_total / np.max(F_w_total)
        else:
            F_w_norm = F_w_total

        if np.max(A_w_total) > 0:
            A_w_norm = A_w_total / np.max(A_w_total)
        else:
            A_w_norm = A_w_total

        # Shift for proper ordering
        wavelengths_fl = np.fft.fftshift(lambda_axis_f)
        fluorescence = np.fft.fftshift(F_w_norm)

        wavelengths_abs = np.fft.fftshift(lambda_axis_a)
        absorption = np.fft.fftshift(A_w_norm)

        # Filter to reasonable wavelength range (550-750 nm for photosynthesis)
        mask_fl = (wavelengths_fl >= 550) & (wavelengths_fl <= 800)
        mask_abs = (wavelengths_abs >= 550) & (wavelengths_abs <= 800)

        wavelengths_fl = wavelengths_fl[mask_fl]
        fluorescence = fluorescence[mask_fl]

        wavelengths_abs = wavelengths_abs[mask_abs]
        absorption = absorption[mask_abs]

        # Ensure monotonic wavelength axes for clean plotting
        if wavelengths_fl.size > 0:
            fl_sort = np.argsort(wavelengths_fl)
            wavelengths_fl = wavelengths_fl[fl_sort]
            fluorescence = fluorescence[fl_sort]

        if wavelengths_abs.size > 0:
            abs_sort = np.argsort(wavelengths_abs)
            wavelengths_abs = wavelengths_abs[abs_sort]
            absorption = absorption[abs_sort]

        return wavelengths_abs, absorption, wavelengths_fl, fluorescence
