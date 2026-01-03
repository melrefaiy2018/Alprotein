#!/usr/bin/env python3
"""
Test Core Workflow: PDB → Site Energies → Hamiltonian → Spectrum

This script tests the core scientific calculations without GUI overhead.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add repo root to path (ensures local package is used)
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Alprotein.core.protein_structure import ProteinStructure
from Alprotein.core.pigment_system import PigmentSystem
from Alprotein.calculators.site_energy_calculator import SiteEnergyCalculator
from Alprotein.calculators.hamiltonian_calculator import HamiltonianCalculator
from Alprotein.calculators.spectra_calculator import SpectraCalculator


# Configuration
# PDB_FILE_PATH = "/Users/mohamed/Documents/Research/Projects/IsiA/manuscript/IsiA_monomer/mcce/mcce_Jun27/most_occ/most_occ_pH8.pdb"

# PDB_FILE_PATH = "/Users/mohamed/Documents/Research/Projects/GitHub_repo/AlProtein/examples/data/most_occ_pH7_IsiA.pdb"
# PDB_FILE_PATH = "/Users/mohamed/Documents/Research/Projects/GitHub_repo/AlProtein/examples/data/extended_most_occ_pH7_IsiA.pdb"
PDB_FILE_PATH ="/Users/mohamed/Documents/Research/Projects/SuperComplex_Kajwal/correct_calculation/CP24/CP24_LYS_gold/fit_HIS_183/CDC/2025/raw/extended_most_occ_pH8.pdb"

EXP_DATA_BASE_PATH = '/Users/mohamed/Documents/Research/Projects/IsiA/IsiA_monomer/MCCE_run/Amber_topology/IsiA_BCR_SQD/CDC/Exp_gabriela/'
PATH_EXP_ABS = f'{EXP_DATA_BASE_PATH}Abs_IsiA_monomer_300K_Gabriela.npy'
PATH_EXP_FLUO = f'{EXP_DATA_BASE_PATH}Fl_IsiA_monomer_300k_Gabriela.npy'

print("="*80)
print("ALPROTEIN CORE WORKFLOW TEST (DEBUG SCRIPT)")
print("="*80)

# Create results directory
RESULTS_DIR = Path(__file__).parent / "results_isia"
RESULTS_DIR.mkdir(exist_ok=True)
print(f"\n[0] Results directory: {RESULTS_DIR}")

# Check file exists
pdb_path = Path(PDB_FILE_PATH)
if not pdb_path.exists():
    print(f"Error: File not found: {pdb_path}")
    sys.exit(1)

# ============================================================================
# STEP 1: Load PDB File
# ============================================================================
print("\n[1] Loading PDB file...")
print(f"    File: {pdb_path}")

# Use from_file_extended to ensure charges are loaded
protein_structure = ProteinStructure.from_file_extended(str(pdb_path), "IsiA")
pigment_system = PigmentSystem(protein_structure)

print(f"    ✓ Loaded successfully")
print(f"    ✓ Found {len(pigment_system.pigments)} pigments")

# Show pigment types
pigment_types = {}
for pig in pigment_system.pigments.values():
    res = pig.get_resname()
    pigment_types[res] = pigment_types.get(res, 0) + 1

print(f"    ✓ Pigment types:")
for res, count in pigment_types.items():
    print(f"      - {res}: {count}")

bg_atoms = pigment_system.get_protein_background_atoms()
print(f"    ✓ Background atoms: {len(bg_atoms)}")

# ============================================================================
# STEP 1b: Load TrEsp Parameters for Coupling Calculations
# ============================================================================
print("\n[1b] Loading TrEsp parameters for coupling calculations...")

from Alprotein.data.parameters import get_coupling_parameters

# Map residue names to TrEsp parameter names
tresp_mapping = {
    'CLA': 'CLA_IPPC',
    'CHL': 'CHL_IPPC',
    # Add more mappings as needed
}

n_loaded = 0
for name, pigment in pigment_system.pigments.items():
    resname = pigment.get_resname()
    if resname in tresp_mapping:
        tresp_param_name = tresp_mapping[resname]
        # Keep the CDC params (already loaded) and add TrEsp params
        cdc_param_name = resname  # CLA or CHL
        try:
            # Load BOTH CDC and TrEsp parameters for this pigment
            pigment._load_calculation_params(tresp_dict_name=tresp_param_name, cdc_dict_name=cdc_param_name)
            n_loaded += 1
        except Exception as e:
            print(f"    ! Warning: Could not load TrEsp params for {name}: {e}")

print(f"    ✓ Loaded TrEsp parameters for {n_loaded}/{len(pigment_system.pigments)} pigments")

# ============================================================================
# STEP 2: Calculate Site Energies
# ============================================================================
print("\n[2] Calculating site energies...")
e0a = 14900.0  # Vacuum energy for CLA
e0b = 15350.0  # Vacuum energy for CHL
dielectric_constant = 2.0  # Protein dielectric constant
site_energy_calculator = SiteEnergyCalculator(
    dielectric_constant=dielectric_constant,
    e0a=e0a,
    # e0b=15674.0
    e0b =e0b
)

site_energies = {}
for name, pigment in pigment_system.pigments.items():
    site_energies[name] = site_energy_calculator.calculate_site_energy(pigment, pigment_system)

# vacuum_energies = pigment_system.get_vacuum_energies()

print(f"    ✓ Calculated {len(site_energies)} site energies")

# Statistics
energies = np.array(list(site_energies.values()))
print(f"    ✓ Statistics:")
print(f"      - Mean: {np.mean(energies):.1f} cm⁻¹")
print(f"      - Std:  {np.std(energies):.1f} cm⁻¹")
print(f"      - Min:  {np.min(energies):.1f} cm⁻¹")
print(f"      - Max:  {np.max(energies):.1f} cm⁻¹")

# ============================================================================
# STEP 3: Construct Hamiltonian
# ============================================================================
print("\n[3] Constructing Hamiltonian...")

hamiltonian_calculator = HamiltonianCalculator(
    pigment_system,
    dielectric_cdc=dielectric_constant,
    dielectric_tresp=1.0,
    f_val=0.72,
    E_0a=e0a,
    E_0b=e0b
)

# Use TrEsp coupling with loaded transition charges
# hamiltonian_calculator.use_dipole_coupling = True

hamiltonian_df = hamiltonian_calculator.construct_hamiltonian(params={'coupling_calc': 'tresp'})
hamiltonian = hamiltonian_df.values

# Verify it's a numpy array
print(f"    ✓ Hamiltonian constructed")
print(f"    ✓ Type: {type(hamiltonian)}")
print(f"    ✓ Shape: {hamiltonian.shape}")
print(f"    ✓ Is numpy array: {isinstance(hamiltonian, np.ndarray)}")

# Debug: recompute site energies using the calculator inside HamiltonianCalculator
ham_site_energies = hamiltonian_calculator.site_energy_calc.calculate_all_site_energies(
    pigment_system, {'CLA': e0a, 'CHL': e0b}
)
diff_by_name = (pd.Series(ham_site_energies) - pd.Series(site_energies)).sort_values(
    key=lambda s: s.abs(), ascending=False
)
print("    ✓ Hamiltonian calculator site energies (top diffs vs Step 2):")
print(diff_by_name.head(5).to_string())

# Statistics on couplings (off-diagonal elements)
n = hamiltonian.shape[0]
off_diag = []
for i in range(n):
    for j in range(i+1, n):
        off_diag.append(abs(hamiltonian[i, j]))

if off_diag:
    print(f"    ✓ Coupling statistics:")
    print(f"      - Mean: {np.mean(off_diag):.1f} cm⁻¹")
    print(f"      - Max:  {np.max(off_diag):.1f} cm⁻¹")
    print(f"      - Min:  {np.min(off_diag):.1f} cm⁻¹")

# Diagonalize
eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
print(f"    ✓ Diagonalized successfully")
print(f"    ✓ Eigenvalue range: {eigenvalues[0]:.1f} to {eigenvalues[-1]:.1f} cm⁻¹")

# Save Hamiltonian as CSV
hamiltonian_csv_path = RESULTS_DIR / "hamiltonian.csv"
hamiltonian_df.to_csv(hamiltonian_csv_path)
pigment_names = list(hamiltonian_df.index)
print(f"    ✓ Saved Hamiltonian to: {hamiltonian_csv_path}")

# Sanity check: diagonal should match site energies by pigment name
diag_by_name = pd.Series(np.diag(hamiltonian), index=hamiltonian_df.index)
site_by_name = pd.Series(site_energies)
site_by_name = site_by_name.reindex(hamiltonian_df.index)
if not np.allclose(diag_by_name.values, site_by_name.values, rtol=1e-7, atol=1e-6):
    diff = (diag_by_name - site_by_name).abs().sort_values(ascending=False)
    top = diff.head(5).to_string()
    raise AssertionError(
        "Hamiltonian diagonal does not match site energies by pigment name.\n"
        f"Top mismatches (abs):\n{top}"
    )
print("    ✓ Hamiltonian diagonal matches site energies by pigment name")

# Also save eigenvalues
eigenvalues_csv_path = RESULTS_DIR / "eigenvalues.csv"
eig_df = pd.DataFrame({'eigenvalue_cm-1': eigenvalues, 'wavelength_nm': 1e7/eigenvalues})
eig_df.to_csv(eigenvalues_csv_path, index_label='state')
print(f"    ✓ Saved eigenvalues to: {eigenvalues_csv_path}")

# ============================================================================
# STEP 4: Domain Definition (Automatic)
# ============================================================================
print("\n[4] Defining domains automatically...")
# manual_domains_config = {
#     0: ["A_CLA_501"],
#     1: ["A_CLA_502"],
#     2: ["A_CLA_503"],
#     3: ["A_CLA_504"],
#     4: ["A_CLA_508", "A_CLA_507", "A_CLA_509", "A_CLA_512", 
#         "A_CLA_511", "A_CLA_516", "A_CLA_510", "A_CLA_518", 
#         "A_CLA_506", "A_CLA_505"],
#     5: ["A_CLA_513"],
#     6: ["A_CLA_517"],
#     7: ["A_CLA_519"]
# }
domain_cutoff = 20.0  # cm-1
domains_dict = hamiltonian_calculator.build_domains(
    cutoff=domain_cutoff,
    hamiltonian=hamiltonian_df,
)

print(f"    Using coupling cutoff: {domain_cutoff} cm-1")
print(f"    ✓ Found {len(domains_dict)} domains")
for domain_id, indices in domains_dict.items():
    domain_names = [pigment_names[idx] for idx in indices]
    print(f"      - Domain {domain_id}: {len(indices)} pigments")
    print(f"        {', '.join(domain_names)}")

# ============================================================================
# STEP 5: Calculate Spectra using Renger's Approach
# ============================================================================
print("\n[5] Calculating absorption and fluorescence spectra...")

# === Use Proper SpectraCalculator with Renger Lineshape Theory ===
use_proper_spectra_calculator = True  # Use the correct Renger implementation

if use_proper_spectra_calculator:
    # === Use Proper SpectraCalculator with Renger Lineshape Theory ===
    print("    Using SpectraCalculator with Renger's phonon-exciton coupling...")

    # Parameters matching CP29_Doran reference
    temperature = 300.0  # K (can be adjusted to match experimental conditions)
    disorder_fwhm = 130.0  # cm⁻¹ - FWHM for inhomogeneous broadening
    n_ensemble = 300  # Number of disorder realizations
    dt = 0.1  # Time step in fs
    t_max = 2000.0  # Maximum time in fs

    print(f"    Parameters:")
    print(f"      - Temperature: {temperature} K")
    print(f"      - Disorder FWHM (inhomogeneous): {disorder_fwhm} cm⁻¹")
    print(f"      - Disorder sigma: {disorder_fwhm/(2*np.sqrt(2*np.log(2))):.1f} cm⁻¹")
    print(f"      - Ensemble size: {n_ensemble}")
    print(f"      - Time axis: 0 to {t_max} fs, step {dt} fs")

    # Initialize SpectraCalculator
    spectra_calc = SpectraCalculator(
        temperature=temperature,
        disorder_sigma=disorder_fwhm,  # Will be converted to sigma internally
        n_ensemble=n_ensemble,
        dt=dt,
        t_max=t_max
    )

    print(f"    ✓ SpectraCalculator initialized")
    print(f"      - Reorganization energy: {spectra_calc.reorganization_energy:.1f} cm⁻¹")

    # Calculate spectra (both absorption and fluorescence)
    wavelengths_abs, absorption, wavelengths_fl, fluorescence = spectra_calc.calculate_spectrum(
        hamiltonian=hamiltonian,
        domains=domains_dict,
        pigment_system=pigment_system,
        include_vibronic=True,  # Include 0-1 vibronic transitions
        use_ensemble=True  # Use ensemble averaging
    )

    # Find absorption peak
    mask_abs = (wavelengths_abs >= 600) & (wavelengths_abs <= 720)
    if np.any(mask_abs):
        peak_idx_abs = np.argmax(absorption[mask_abs])
        peak_wl = wavelengths_abs[mask_abs][peak_idx_abs]
    else:
        peak_idx_abs = np.argmax(absorption)
        peak_wl = wavelengths_abs[peak_idx_abs]

    # Find fluorescence peak
    mask_fl = (wavelengths_fl >= 600) & (wavelengths_fl <= 720)
    if np.any(mask_fl):
        peak_idx_fl = np.argmax(fluorescence[mask_fl])
        fl_peak_wl = wavelengths_fl[mask_fl][peak_idx_fl]
    else:
        peak_idx_fl = np.argmax(fluorescence)
        fl_peak_wl = wavelengths_fl[peak_idx_fl]

    # Use absorption spectrum for backward compatibility
    wavelengths = wavelengths_abs
    spectrum = absorption

else:
    # === Simplified Gaussian Lineshape Approach ===
    print("    Using simplified Gaussian broadening approach...")

    # Parameters
    temperature = 300.0  # K
    disorder_fwhm = 250.0  # cm⁻¹ - FWHM for inhomogeneous broadening
    homogeneous_fwhm = 500.0  # cm⁻¹ - FWHM for homogeneous broadening (increased to merge peaks)
    n_ensemble = 300  # Number of disorder realizations

    # Convert FWHM to sigma
    disorder_sigma = disorder_fwhm / (2 * np.sqrt(2 * np.log(2)))
    homogeneous_sigma = homogeneous_fwhm / (2 * np.sqrt(2 * np.log(2)))

    print(f"    Parameters:")
    print(f"      - Temperature: {temperature} K")
    print(f"      - Disorder FWHM (inhomogeneous): {disorder_fwhm} cm⁻¹ (σ = {disorder_sigma:.1f} cm⁻¹)")
    print(f"      - Homogeneous FWHM: {homogeneous_fwhm} cm⁻¹ (σ = {homogeneous_sigma:.1f} cm⁻¹)")
    print(f"      - Ensemble size: {n_ensemble}")

    # Create wavelength axis (nm)
    wavelengths = np.linspace(550, 750, 2000)

    # Initialize spectrum
    spectrum = np.zeros_like(wavelengths)

    # Get pigment dipole moments (if available)
    pigment_names = list(pigment_system.pigments.keys())
    dipole_vectors = []
    for name in pigment_names:
        pig = pigment_system.pigments[name]
        try:
            # Try to get transition dipole moment vector
            dipole = pig.get_transition_dipole_moment()
            dipole_vectors.append(dipole)
        except:
            # Fallback: assume unit dipole along z-axis
            dipole_vectors.append(np.array([0.0, 0.0, 1.0]))

    dipole_vectors = np.array(dipole_vectors)

    # Ensemble averaging with inhomogeneous broadening (diagonal disorder)
    for n in range(n_ensemble):
        # Add diagonal disorder to site energies (INHOMOGENEOUS BROADENING)
        rng = np.random.RandomState(seed=n)
        site_energies_diag = np.diag(hamiltonian)
        disorder = rng.normal(0, disorder_sigma, size=len(site_energies_diag))

        H_disorder = hamiltonian.copy()
        np.fill_diagonal(H_disorder, site_energies_diag + disorder)

        # Diagonalize to get exciton states
        eigenvalues, eigenvectors = np.linalg.eigh(H_disorder)

        # Calculate oscillator strengths using proper transition dipole moments
        # μ_ex = Σ_n c_n μ_n (coherent sum of site dipoles weighted by exciton coefficients)
        for state_idx in range(len(eigenvalues)):
            E_ex = eigenvalues[state_idx]  # Exciton energy in cm⁻¹
            coeffs = eigenvectors[:, state_idx]  # Exciton coefficients

            # Calculate exciton transition dipole (vectorial sum)
            mu_ex_vector = np.sum(coeffs[:, np.newaxis] * dipole_vectors, axis=0)
            mu_ex_squared = np.dot(mu_ex_vector, mu_ex_vector)  # |μ_ex|²

            # Skip dark states (negligible oscillator strength)
            if mu_ex_squared < 1e-6:
                continue

            # Add HOMOGENEOUS BROADENING: Lorentzian or Gaussian lineshape
            # Using Gaussian for consistency
            # Convert exciton energy to wavelength
            lambda_ex = 1e7 / E_ex  # nm

            # Gaussian lineshape in wavelength space
            # Note: FWHM in cm⁻¹ needs conversion to nm for wavelength-space Gaussian
            # Δλ ≈ λ²/(1e7) * Δν (where Δν is in cm⁻¹)
            sigma_lambda = (lambda_ex**2 / 1e7) * homogeneous_sigma

            lineshape = mu_ex_squared * np.exp(-0.5 * ((wavelengths - lambda_ex) / sigma_lambda)**2)
            spectrum += lineshape

    # Normalize by ensemble size
    spectrum = spectrum / n_ensemble

    # Apply overall normalization
    if np.max(spectrum) > 0:
        spectrum = spectrum / np.max(spectrum)

    # Find peak
    peak_idx = np.argmax(spectrum)
    peak_wl = wavelengths[peak_idx]

print(f"    ✓ Absorption spectrum calculated")
print(f"    ✓ Absorption wavelength range: {wavelengths_abs.min():.1f} - {wavelengths_abs.max():.1f} nm")
print(f"    ✓ Absorption peak wavelength: {peak_wl:.1f} nm")
print(f"    ✓ Absorption peak wavenumber: {1e7/peak_wl:.1f} cm⁻¹")

print(f"    ✓ Fluorescence spectrum calculated")
print(f"    ✓ Fluorescence wavelength range: {wavelengths_fl.min():.1f} - {wavelengths_fl.max():.1f} nm")
print(f"    ✓ Fluorescence peak wavelength: {fl_peak_wl:.1f} nm")
print(f"    ✓ Fluorescence peak wavenumber: {1e7/fl_peak_wl:.1f} cm⁻¹")
print(f"    ✓ Stokes shift: {1e7/peak_wl - 1e7/fl_peak_wl:.1f} cm⁻¹ ({fl_peak_wl - peak_wl:.1f} nm)")

# For comparison, also print the mean site energy
mean_site_energy = np.mean(np.diag(hamiltonian))
mean_wavelength = 1e7 / mean_site_energy
print(f"    ✓ Mean site energy: {mean_site_energy:.1f} cm⁻¹ ({mean_wavelength:.1f} nm)")

# ============================================================================
# STEP 6: Plot Results
# ============================================================================
print("\n[6] Plotting results...")

# === FIGURE 1: Absorption Spectrum ===
fig1, ax1 = plt.subplots(figsize=(10, 6))

# Try to load and plot experimental data
try:
    if Path(PATH_EXP_ABS).exists():
        exp_data = np.load(PATH_EXP_ABS)
        # Check shape - assuming [wavelength, intensity] columns
        if exp_data.ndim == 2 and exp_data.shape[1] >= 2:
            exp_wl = exp_data[:, 0]
            exp_abs = exp_data[:, 1]
            # Normalize
            if np.max(exp_abs) > 0:
                exp_abs = exp_abs / np.max(exp_abs)
            ax1.plot(exp_wl, exp_abs, 'k-', linewidth=2, alpha=0.7, label='Experimental')
            print(f"    ✓ Loaded experimental data from {PATH_EXP_ABS}")
        else:
            print(f"    ! Experimental data shape not recognized: {exp_data.shape}")
    else:
        print(f"    ! Experimental data file not found: {PATH_EXP_ABS}")
except Exception as e:
    print(f"    ! Failed to load experimental data: {e}")

ax1.plot(wavelengths, spectrum, 'b-', linewidth=2.5, label='Simulated (Renger)')
ax1.set_xlabel('Wavelength (nm)', fontsize=12)
ax1.set_ylabel('Absorption (normalized)', fontsize=12)
ax1.set_title(f'IsiA Monomer Absorption Spectrum (Peak: {peak_wl:.1f} nm)', fontsize=14)
ax1.legend(fontsize=11)
ax1.set_xlim(550, 750)
ax1.grid(alpha=0.3)
fig1.tight_layout()

# Save absorption spectrum figure
output_file1 = RESULTS_DIR / "absorption_spectrum.png"
fig1.savefig(output_file1, dpi=150, bbox_inches='tight')
print(f"    ✓ Saved absorption plot to: {output_file1}")

# Also save spectrum data as CSV
spectrum_csv_path = RESULTS_DIR / "absorption_spectrum.csv"
spectrum_df = pd.DataFrame({'wavelength_nm': wavelengths, 'absorption': spectrum})
spectrum_df.to_csv(spectrum_csv_path, index=False)
print(f"    ✓ Saved spectrum data to: {spectrum_csv_path}")

# === FIGURE 2: Fluorescence Spectrum ===
fig2, ax2 = plt.subplots(figsize=(10, 6))

# Try to load and plot experimental fluorescence data
try:
    if Path(PATH_EXP_FLUO).exists():
        exp_fl_data = np.load(PATH_EXP_FLUO)
        if exp_fl_data.ndim == 2 and exp_fl_data.shape[1] >= 2:
            exp_fl_wl = exp_fl_data[:, 0]
            exp_fl_int = exp_fl_data[:, 1]
            # Normalize
            if np.max(exp_fl_int) > 0:
                exp_fl_int = exp_fl_int / np.max(exp_fl_int)
            ax2.plot(exp_fl_wl, exp_fl_int, 'k-', linewidth=2, alpha=0.7, label='Experimental')
            print(f"    ✓ Loaded experimental fluorescence data from {PATH_EXP_FLUO}")
        else:
            print(f"    ! Experimental fluorescence data shape not recognized: {exp_fl_data.shape}")
    else:
        print(f"    ! Experimental fluorescence file not found: {PATH_EXP_FLUO}")
except Exception as e:
    print(f"    ! Failed to load experimental fluorescence data: {e}")

ax2.plot(wavelengths_fl, fluorescence, 'r-', linewidth=2.5, label='Simulated (Renger)')
ax2.set_xlabel('Wavelength (nm)', fontsize=12)
ax2.set_ylabel('Fluorescence (normalized)', fontsize=12)
ax2.set_title(f'IsiA Monomer Fluorescence Spectrum (Peak: {fl_peak_wl:.1f} nm)', fontsize=14)
ax2.legend(fontsize=11)
ax2.set_xlim(550, 750)
ax2.grid(alpha=0.3)
fig2.tight_layout()

# Save fluorescence spectrum figure
output_file2 = RESULTS_DIR / "fluorescence_spectrum.png"
fig2.savefig(output_file2, dpi=150, bbox_inches='tight')
print(f"    ✓ Saved fluorescence plot to: {output_file2}")

# Save fluorescence data as CSV
fl_spectrum_csv_path = RESULTS_DIR / "fluorescence_spectrum.csv"
fl_spectrum_df = pd.DataFrame({'wavelength_nm': wavelengths_fl, 'fluorescence': fluorescence})
fl_spectrum_df.to_csv(fl_spectrum_csv_path, index=False)
print(f"    ✓ Saved fluorescence data to: {fl_spectrum_csv_path}")

# === FIGURE 3: Site Energy Shifts ===
fig3, ax3 = plt.subplots(figsize=(12, 6))

# Calculate site energy shifts (relative to E0a)
pigment_names = list(site_energies.keys())
energies = np.array(list(site_energies.values()))
# if the pigment is CHL, adjust the reference energy
E0a = e0a
E0b = e0b

shifts = np.array([
    energy - (E0b if "CHL" in name else E0a)
    for name, energy in zip(pigment_names, energies)
])

# Create bar colors based on shift direction
colors = ['#d62728' if s > 0 else '#1f77b4' for s in shifts]

# Create x positions and labels
x_pos = np.arange(len(pigment_names))
# Simplify labels (extract just the residue number)
labels = [name.split('_')[-1] for name in pigment_names]

bars = ax3.bar(x_pos, shifts, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)

# Add horizontal line at zero
ax3.axhline(y=0, color='black', linestyle='-', linewidth=1)

# Add value labels on bars
for i, (bar, shift) in enumerate(zip(bars, shifts)):
    height = bar.get_height()
    va = 'bottom' if height >= 0 else 'top'
    offset = 3 if height >= 0 else -3
    ax3.annotate(f'{shift:.0f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, offset),
                    textcoords="offset points",
                    ha='center', va=va, fontsize=9)

ax3.set_xlabel('Pigment (CLA residue number)', fontsize=12)
ax3.set_ylabel(f'Site Energy Shift (cm⁻¹) relative to E₀ = {E0a} cm⁻¹', fontsize=12)
ax3.set_title('Site Energy Shifts from CDC Calculation', fontsize=14)
ax3.set_xticks(x_pos)
ax3.set_xticklabels(labels, rotation=45, ha='right')
ax3.grid(axis='y', alpha=0.3)

# Add statistics annotation
stats_text = f'Mean shift: {np.mean(shifts):.1f} cm⁻¹\nStd: {np.std(shifts):.1f} cm⁻¹\nRange: {np.min(shifts):.0f} to {np.max(shifts):.0f} cm⁻¹'
ax3.text(0.98, 0.98, stats_text, transform=ax3.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

fig3.tight_layout()

# Save site energy shift figure
output_file3 = RESULTS_DIR / "site_energy_shifts.png"
fig3.savefig(output_file3, dpi=150, bbox_inches='tight')
print(f"    ✓ Saved site energy shift plot to: {output_file3}")

# Also save site energies as CSV
site_energies_csv_path = RESULTS_DIR / "site_energies.csv"
site_energies_df = pd.DataFrame({
    'pigment': pigment_names,
    'site_energy_cm-1': energies,
    'shift_cm-1': shifts,
    'wavelength_nm': 1e7 / energies
})
site_energies_df.to_csv(site_energies_csv_path, index=False)
print(f"    ✓ Saved site energies to: {site_energies_csv_path}")

# Show plots
plt.show()

# ============================================================================
# Success!
# ============================================================================
print("\n" + "="*80)
print("✓ ALL TESTS PASSED!")
print("="*80)
print("\nCore implementation is working correctly.")
