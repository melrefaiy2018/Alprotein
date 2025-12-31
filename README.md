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

Here are some basic usage examples to get you started:

```python
from Alprotein import ProteinStructure, PigmentSystem, ChlorophyllA, HamiltonianCalculator

# Load a protein structure
protein = ProteinStructure.from_file("structure.pdb", name="MyComplex")

# Identify pigments and build pigment systems
pigment_system = PigmentSystem(protein)
pigment_system.add_pigments_by_residue(
    resname="CLA",
    pigment_class=ChlorophyllA,
    tresp_dict_name="CLA_IPPC",
    cdc_dict_name="CLA",
)

# Compute pigment site energies using CDC
site_energies = pigment_system.compute_site_energies()

# Compute pigment to pigment couplings using TrEsp or dipole approximations
couplings = pigment_system.compute_couplings(method="TrEsp")

# Construct and diagonalize excitonic Hamiltonians
calculator = HamiltonianCalculator(pigment_system)
H = calculator.construct_hamiltonian()
eigs, vecs = calculator.diagonalize_hamiltonian(H)

print("Hamiltonian shape:", H.shape)
print("Exciton energies (cm^-1):", eigs)
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
  year={2025},
  version={0.2.0},
  url={https://github.com/melrefaiy2018/Alprotein}
}
```

## Contributing

We welcome contributions! Please read our [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests.

## Contact

For any questions or support, please reach out to [your_email@example.com].
