<div align="center">
  <img src="docs/assets/logo.png" alt="Alprotein Logo" width="600"/>
</div>

<div align="center">

[![Build Status](https://travis-ci.org/melrefaiy2018/Alprotein-Alpha.svg?branch=main)](https://travis-ci.org/melrefaiy2018/Alprotein-Alpha)
[![codecov](https://codecov.io/gh/melrefaiy2018/Alprotein-Alpha/branch/main/graph/badge.svg)](https://codecov.io/gh/melrefaiy2018/Alprotein-Alpha)
[![PyPI version](https://badge.fury.io/py/Alprotein.svg)](https://badge.fury.io/py/Alprotein)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/melrefaiy2018/Alprotein-Alpha.svg?style=social&label=Star)](https://github.com/melrefaiy2018/Alprotein-Alpha)

</div>

# Alprotein: A modular toolkit for excitonic Hamiltonians and optical spectra of molecular aggregates

Alprotein is a modular Python toolkit for calculating excitonic Hamiltonians, site energies, and optical spectra of pigment protein complexes from structural data.

It is designed for computational biophysics workflows where you start from a PDB structure, identify pigments, compute pigment site energies and couplings, and then derive excitonic properties and spectra.

## Key capabilities

- Load protein structures from PDB files
- Identify pigments and build pigment systems
- Compute pigment site energies using CDC
- Compute pigment to pigment couplings using TrEsp or dipole approximations
- Construct and diagonalize excitonic Hamiltonians
- Export results for analysis and downstream modeling

## Installation

### From source

```bash
git clone https://github.com/melrefaiy2018/Alprotein.git
conda create -n alprotein python=3.11 
cd Alprotein
pip install -e .
```

### Requirements

- Python 3.8 or newer
- NumPy, SciPy, pandas
- BioPython
- matplotlib

## Usage Examples

Here is a complete workflow example to get you started. This demonstrates how to calculate excitonic Hamiltonians and optical spectra from a PDB structure:

```python
import numpy as np
import matplotlib.pyplot as plt
from Alprotein import (
    ProteinStructure,
    PigmentSystem,
    ChlorophyllA,
    HamiltonianCalculator,
    SpectraCalculator
)

# ============================================================================
# STEP 1: Load Protein Structure
# ============================================================================
print("[1] Loading protein structure...")
protein = ProteinStructure.from_file("structure.pdb", name="MyComplex")
print(f"    ✓ Loaded structure with {len(protein.atoms)} atoms")

# ============================================================================
# STEP 2: Build Pigment System
# ============================================================================
print("\n[2] Building pigment system...")
pigment_system = PigmentSystem(protein)

# Add chlorophyll A pigments from the structure
pigment_system.add_pigments_by_residue(
    resname="CLA",              # Residue name in PDB file
    pigment_class=ChlorophyllA, # Pigment type class
    tresp_dict_name="CLA_IPPC", # Transition charge parameter set
    cdc_dict_name="CLA",        # CDC parameter set
)

print(f"    ✓ Identified {len(pigment_system.pigments)} pigments")

# Compute site energies using CDC method
site_energies = pigment_system.compute_site_energies()
print(f"    ✓ Computed site energies (mean: {np.mean(site_energies):.1f} cm⁻¹)")

# Compute pigment-pigment couplings using TrEsp method
couplings = pigment_system.compute_couplings(method="TrEsp")
print(f"    ✓ Computed couplings (max: {np.max(np.abs(couplings)):.1f} cm⁻¹)")

# ============================================================================
# STEP 3: Construct Hamiltonian
# ============================================================================
print("\n[3] Constructing Hamiltonian...")

hamiltonian_calculator = HamiltonianCalculator(
    pigment_system,
    dielectric_cdc=2.0,      # Dielectric constant for CDC calculation
    dielectric_tresp=1.0,    # Dielectric constant for TrEsp coupling
    f_val=0.72,              # Oscillator strength
    E_0a=14900.0,            # Reference energy for Qy transition (cm⁻¹)
    E_0b=15674.0             # Reference energy for Qx transition (cm⁻¹)
)

# Build the Hamiltonian matrix
hamiltonian_df = hamiltonian_calculator.construct_hamiltonian(
    params={'coupling_calc': 'tresp'}  # Use TrEsp for couplings
)
hamiltonian = hamiltonian_df.values
print(f"    ✓ Hamiltonian constructed ({hamiltonian.shape[0]}×{hamiltonian.shape[1]})")

# ============================================================================
# STEP 4: Define Excitonic Domains
# ============================================================================
print("\n[4] Defining excitonic domains...")

domain_cutoff = 20.0  # Coupling threshold in cm⁻¹
domains_dict = hamiltonian_calculator.build_domains(
    cutoff=domain_cutoff,
    hamiltonian=hamiltonian_df,
)
print(f"    ✓ Found {len(domains_dict)} domains")

# ============================================================================
# STEP 5: Calculate Optical Spectra
# ============================================================================
print("\n[5] Calculating optical spectra...")

# Set calculation parameters
temperature = 300.0       # Temperature in Kelvin
disorder_fwhm = 130.0     # Inhomogeneous broadening FWHM (cm⁻¹)
n_ensemble = 300          # Number of disorder realizations
dt = 0.1                  # Time step in femtoseconds
t_max = 2000.0            # Maximum time in femtoseconds

print(f"    Parameters:")
print(f"      - Temperature: {temperature} K")
print(f"      - Disorder FWHM: {disorder_fwhm} cm⁻¹")
print(f"      - Ensemble size: {n_ensemble}")
print(f"      - Time range: 0-{t_max} fs (step: {dt} fs)")

# Initialize spectra calculator
spectra_calc = SpectraCalculator(
    temperature=temperature,
    disorder_sigma=disorder_fwhm,  # Converted to sigma internally
    n_ensemble=n_ensemble,
    dt=dt,
    t_max=t_max
)
print(f"    ✓ SpectraCalculator initialized")
print(f"      - Reorganization energy: {spectra_calc.reorganization_energy:.1f} cm⁻¹")

# Calculate absorption and fluorescence spectra
wavelengths_abs, absorption, wavelengths_fl, fluorescence = spectra_calc.calculate_spectrum(
    hamiltonian=hamiltonian,
    domains=domains_dict,
    pigment_system=pigment_system,
    include_vibronic=True,  # Include 0-1 vibronic transitions
    use_ensemble=True       # Use ensemble averaging
)
print(f"    ✓ Spectra calculated")

# ============================================================================
# STEP 6: Analyze and Visualize Results
# ============================================================================
print("\n[6] Analyzing results...")

# Find absorption peak
mask_abs = (wavelengths_abs >= 600) & (wavelengths_abs <= 720)
if np.any(mask_abs):
    peak_idx_abs = np.argmax(absorption[mask_abs])
    peak_wl_abs = wavelengths_abs[mask_abs][peak_idx_abs]
else:
    peak_idx_abs = np.argmax(absorption)
    peak_wl_abs = wavelengths_abs[peak_idx_abs]

# Find fluorescence peak
mask_fl = (wavelengths_fl >= 600) & (wavelengths_fl <= 720)
if np.any(mask_fl):
    peak_idx_fl = np.argmax(fluorescence[mask_fl])
    peak_wl_fl = wavelengths_fl[mask_fl][peak_idx_fl]
else:
    peak_idx_fl = np.argmax(fluorescence)
    peak_wl_fl = wavelengths_fl[peak_idx_fl]

print(f"    ✓ Absorption peak: {peak_wl_abs:.1f} nm")
print(f"    ✓ Fluorescence peak: {peak_wl_fl:.1f} nm")
print(f"    ✓ Stokes shift: {peak_wl_fl - peak_wl_abs:.1f} nm")

# Plot spectra
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(wavelengths_abs, absorption, label='Absorption', linewidth=2)
ax.plot(wavelengths_fl, fluorescence, label='Fluorescence', linewidth=2)
ax.axvline(peak_wl_abs, color='C0', linestyle='--', alpha=0.5)
ax.axvline(peak_wl_fl, color='C1', linestyle='--', alpha=0.5)
ax.set_xlabel('Wavelength (nm)', fontsize=12)
ax.set_ylabel('Intensity (arb. units)', fontsize=12)
ax.set_title('Optical Spectra', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('spectra.png', dpi=300)
print(f"\n    ✓ Spectra saved to 'spectra.png'")

print("\n✅ Calculation complete!")
```

### Quick Start Guide

For a minimal working example:

```python
import numpy as np
from Alprotein import ProteinStructure, PigmentSystem, ChlorophyllA

# 1. Load structure and build pigment system
protein = ProteinStructure.from_file("structure.pdb", name="MyComplex")
pigment_system = PigmentSystem(protein)
pigment_system.add_pigments_by_residue(
    resname="CLA",
    pigment_class=ChlorophyllA,
    tresp_dict_name="CLA_IPPC",
    cdc_dict_name="CLA"
)

# 2. Compute energies and couplings
site_energies = pigment_system.compute_site_energies()
couplings = pigment_system.compute_couplings(method="TrEsp")

print(f"Found {len(pigment_system.pigments)} pigments")
print(f"Mean site energy: {np.mean(site_energies):.1f} cm⁻¹")
```

## Project structure

- `Alprotein/core/` core system objects such as structure and pigment handling
- `Alprotein/calculators/` site energy, coupling, and Hamiltonian calculators
- `Alprotein/data/` parameter sets for supported pigment types
- `examples/` runnable examples and workflows


## License

MIT License.

## Citation

If you use Alprotein in published work, please cite:

```bibtex
@software{Alprotein,
  title={Alprotein: A modular toolkit for excitonic Hamiltonians and optical spectra of molecular aggregates},  author={Mohamed Elrefaiy and Bailey Raber and Doran Raccah},
  author={Mohamed Elrefaiy and Bailey Raber and Doran Raccah},
  year={2025},
  version={0.2.0},
  url={https://github.com/melrefaiy2018/Alprotein}
}
```

## Contributing

We welcome contributions! Please read our [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests.

## Contact

For any questions or support, please reach out to [melrefaiy@utexas.edu].
