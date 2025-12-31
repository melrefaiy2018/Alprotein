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

# Alprotein

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
git clone https://github.com/melrefaiy2018/Alprotein-Alpha.git
cd Alprotein-Alpha
pip install -e .
```

### Requirements

- Python 3.8 or newer
- NumPy, SciPy, pandas
- BioPython
- matplotlib

## Quick start

```python
from Alprotein import ProteinStructure, PigmentSystem, ChlorophyllA, HamiltonianCalculator

protein = ProteinStructure.from_file("structure.pdb", name="MyComplex")

pigment_system = PigmentSystem(protein)
pigment_system.add_pigments_by_residue(
    resname="CLA",
    pigment_class=ChlorophyllA,
    tresp_dict_name="CLA_IPPC",
    cdc_dict_name="CLA",
)

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

## Documentation

For a more detailed technical explanation of the computation pipeline, see `agent.md`.

## License

MIT License.

## Citation

If you use Alprotein in published work, please cite:

```bibtex
@software{Alprotein,
  title={Alprotein: A modular toolkit for pigment-protein complex excitonics and spectra},
  author={Mohamed Elrefaiy and Bailey Raber and Doran Raccah},
  year={2025},
  version={0.2.0},
  url={https://github.com/melrefaiy2018/Alprotein-Alpha}
}
```
