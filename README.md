# deltahf

[![Python versions](https://img.shields.io/badge/python-≥3.10-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

## Introduction

**deltahf** estimates gas-phase standard heats of formation (ΔHf°) from atom equivalent energies and semi-empirical quantum chemistry. This approach was described by Cawkwell et al.<sup>1</sup> using DFTB; this implementation replaces DFTB with [GFN2-xTB](https://github.com/grimme-lab/xtb), requiring the atom equivalent energies to be re-parameterized.

The core equation is:

```
ΔHf° = u_xtb − Σ(nl × εl)
```

where `u_xtb` is the xTB total energy (in kcal/mol), `nl` is the count of atom type `l`, and `εl` is the fitted atom equivalent energy. Two models are available: a **4-parameter** model using elemental stoichiometry (C, H, N, O) and a **7-parameter** model that additionally distinguishes multiply-bonded atoms (C', N', O') using a bond order threshold of 1.25.

## Quick Start

```bash
# Install
conda env create -f environment.yml
conda activate deltahf
pip install -e ".[dev]"

# Sanity check with a single molecule
python -c "
from deltahf.pipeline import process_molecule
r = process_molecule('C', n_conformers=1, name='methane')
print(f'xTB energy: {r.xtb_energy:.6f} Eh ({r.xtb_energy_kcal:.2f} kcal/mol)')
print(f'Atom counts: {r.atom_counts_4param}')
"

# Fit atom equivalents on 103 training molecules
python -m deltahf fit -i deltahf/data/training_data.csv --model both --kfold 10 --n-conformers 1 -o params.json
```

## Installation

**Via conda (recommended):**

The environment includes RDKit and xTB, which are most easily obtained through conda-forge.

```bash
conda env create -f environment.yml
conda activate deltahf
pip install -e ".[dev]"
```

**From source (if RDKit and xTB are already available):**

```bash
pip install -e .
```

## Usage

deltahf provides two CLI subcommands: `fit` (parameterize atom equivalents) and `predict` (apply fitted parameters to new molecules).

```bash
python -m deltahf [fit|predict] [options]
```

### Subcommand: `fit`

Fit atom equivalent energies to a training set of molecules with known experimental ΔHf° values.

```bash
python -m deltahf fit -i training.csv --model both --kfold 10 --n-conformers 5 -o params.json
```

| Option | Description | Default |
|--------|-------------|---------|
| `--input, -i` | CSV file with `smiles` and `exp_dhf_kcal_mol` columns (required). | — |
| `--model` | Which model to fit: `4param`, `7param`, or `both`. | `both` |
| `--kfold` | Number of cross-validation folds. | `10` |
| `--n-conformers` | Number of lowest-energy RDKit conformers to optimize with xTB. | `5` |
| `--output, -o` | Output JSON file for fitted epsilon values. | — |

### Subcommand: `predict`

Predict ΔHf° for new molecules using previously fitted atom equivalent energies.

```bash
python -m deltahf predict -i molecules.csv --epsilon params.json --model 4param --n-conformers 5 -o results.csv
```

| Option | Description | Default |
|--------|-------------|---------|
| `--input, -i` | CSV file with a `smiles` column (required). | — |
| `--epsilon` | JSON file with fitted atom equivalent energies (required). | — |
| `--model` | Which model to use: `4param` or `7param`. | `4param` |
| `--n-conformers` | Number of conformers to optimize with xTB. | `5` |
| `--output, -o` | Output CSV with predicted ΔHf° values. | — |

## Examples

### Example 1: Process a Single Molecule

```python
from deltahf.pipeline import process_molecule

result = process_molecule("C[N+](=O)[O-]", n_conformers=1, name="nitromethane")
print(f"xTB energy: {result.xtb_energy:.6f} Eh")
print(f"4-param counts: {result.atom_counts_4param}")
print(f"7-param counts: {result.atom_counts_7param}")
```

### Example 2: Fit Atom Equivalents on the Training Set

The package includes a curated dataset of 103 molecules from Table 1 of the reference paper with experimental ΔHf° values.

```bash
python -m deltahf fit \
    -i deltahf/data/training_data.csv \
    --model both \
    --kfold 10 \
    --n-conformers 1 \
    -o params.json
```

This processes all 103 molecules through the pipeline (SMILES → RDKit conformers → xTB optimization), fits atom equivalents by least squares, and reports RMSD, maximum deviation, and 10-fold cross-validation error for both models:

```
103/103 molecules completed successfully

4-parameter atom equivalents (kcal/mol):
  C: -1350.394263
  H: -314.769803
  N: -1847.189318
  O: -2512.025581
  RMSD: 17.52 kcal/mol
  Max deviation: 77.45 kcal/mol
  CV error (10-fold): 325.99 (kcal/mol)^2

7-parameter atom equivalents (kcal/mol):
  C: -1364.507773
  H: -308.435763
  N: -1856.918408
  O: -2518.577026
  C_prime: -1355.652429
  N_prime: -1843.470587
  O_prime: -2510.955339
  RMSD: 14.71 kcal/mol
  Max deviation: 70.02 kcal/mol
  CV error (10-fold): 262.58 (kcal/mol)^2
```

### Example 3: Predict ΔHf° for New Molecules

Given a CSV file `new_molecules.csv`:

```csv
smiles,name
CC(=O)O,acetic acid
c1ccccc1,benzene
```

```bash
python -m deltahf predict \
    -i new_molecules.csv \
    --epsilon params.json \
    --model 4param \
    --n-conformers 3 \
    -o predictions.csv
```

## Pipeline

For each molecule, deltahf performs the following steps:

1. **SMILES parsing** — Convert SMILES to an RDKit mol object with explicit hydrogens
2. **Conformer generation** — Generate 3D conformers via ETKDG and rank by MMFF94 energy
3. **xTB optimization** — Optimize the lowest *n* conformers with GFN2-xTB
4. **ΔHf° prediction** — Apply the atom equivalent formula using the lowest xTB energy

## Dependencies

- [Python](https://www.python.org/) >= 3.10
- [RDKit](https://www.rdkit.org/) — molecular representation, 3D embedding, ETKDG conformer generation
- [xTB](https://github.com/grimme-lab/xtb) — GFN2 semi-empirical optimization (CLI)
- [NumPy](https://numpy.org/) / [SciPy](https://scipy.org/) — numerical fitting
- [pandas](https://pandas.pydata.org/) — CSV I/O

## Testing

```bash
# Unit tests (no xTB required)
pytest -m "not integration"

# Integration tests (requires xTB on PATH)
pytest -m integration

# All tests
pytest
```

## References

1. Cawkwell, M. J.; Manner, V. W.; Kress, J. D. *J. Chem. Inf. Model.* **2021**, *61*, 3337–3347 [**DOI:** 10.1021/acs.jcim.1c00312](https://doi.org/10.1021/acs.jcim.1c00312)

---
License: [MIT](https://opensource.org/licenses/MIT)
