# deltahf

[![Python versions](https://img.shields.io/badge/python-≥3.10-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

## Introduction

**deltahf** estimates gas-phase standard heats of formation (ΔHf°) from atom equivalent energies and semi-empirical quantum chemistry. This approach was described by Cawkwell et al.<sup>1</sup> using DFTB; this implementation replaces DFTB with [GFN2-xTB](https://github.com/grimme-lab/xtb), requiring the atom equivalent energies to be re-parameterized.

The core equation is:

```
ΔHf° = u_xtb − Σ(nl × εl)
```

where `u_xtb` is the xTB total energy (in kcal/mol), `nl` is the count of atom type `l`, and `εl` is the fitted atom equivalent energy.

### Atom Equivalent Models

Four atom classification schemes are available:

| Model | Parameters | Classification |
|-------|-----------|----------------|
| `4param` | 4 | Elemental stoichiometry: C, H, N, O |
| `7param` | 7 | Adds multiply-bonded variants: C', N', O' (bond order > 1.25) |
| `hybrid` | up to 10 | RDKit hybridization: C/N/O split by sp3, sp2, sp |
| `extended` | up to 15 | Hybridization + H-count for carbon (e.g. C_sp3_3H, C_sp2_1H) |

Parameters with zero training examples are automatically excluded from fitting.

## Quick Start

```bash
# Install
conda env create -f environment.yml
conda activate deltahf
pip install -e ".[dev]"

# Fit atom equivalents on the training set
python -m deltahf fit -i deltahf/data/training_data.csv --model all --kfold 10 --n-conformers 1 -o params.json
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
python -m deltahf fit -i training.csv --model all --kfold 10 --n-conformers 5 -o params.json
```

| Option | Description | Default |
|--------|-------------|---------|
| `--input, -i` | CSV file with `smiles` and `exp_dhf_kcal_mol` columns (required). | — |
| `--model` | Which model(s) to fit: `4param`, `7param`, `hybrid`, `extended`, `both`, or `all`. | `both` |
| `--kfold` | Number of cross-validation folds. | `10` |
| `--n-conformers` | Number of lowest-energy RDKit conformers to optimize with xTB. | `5` |
| `--output, -o` | Output JSON file for fitted epsilon values. | — |
| `--csv` | Output CSV with training data and per-molecule predictions. | — |
| `--use-xtb-wbos` | Use xTB Wiberg bond orders (instead of RDKit) for 7-param classification. | — |
| `--cache-dir` | Directory for caching xTB results (enables restart capability). | — |
| `--verbose, -v` | Print per-molecule details instead of a progress bar. | — |

### Subcommand: `predict`

Predict ΔHf° for new molecules using previously fitted atom equivalent energies.

```bash
python -m deltahf predict -i molecules.csv --epsilon params.json --model 4param --n-conformers 5 -o results.csv
```

| Option | Description | Default |
|--------|-------------|---------|
| `--input, -i` | CSV file with a `smiles` column (required). | — |
| `--epsilon` | JSON file with fitted atom equivalent energies (required). | — |
| `--model` | Which model to use: `4param`, `7param`, `hybrid`, or `extended`. | `4param` |
| `--n-conformers` | Number of conformers to optimize with xTB. | `5` |
| `--output, -o` | Output CSV with predicted ΔHf° values. | — |
| `--cache-dir` | Directory for caching xTB results. | — |
| `--verbose, -v` | Print per-molecule details instead of a progress bar. | — |

## Training Data

`deltahf/data/training_data.csv` contains **314 CHNO molecules** with experimental ΔHf° values from two literature sources:

| Source | Count | Description |
|--------|-------|-------------|
| Cawkwell et al. (2021)<sup>1</sup> | 102 | Energetic CHNO molecules + small reference compounds |
| Yalamanchi et al. (2020)<sup>2</sup> | 212 | Cyclic hydrocarbons (CH only) |

Each molecule is assigned a `category` label:

| Category | Count | Description |
|----------|-------|-------------|
| `energetic` | 45 | Explosives, nitro/azide/nitroso compounds |
| `small_CHNO` | 57 | Small reference molecules (methane, water, ethanol, etc.) |
| `cyclic_HC` | 189 | Cyclic hydrocarbons (cyclopentanes through naphthalenes) |
| `strained_3ring` | 13 | Cyclopropane derivatives and quadricyclane |
| `large_HC` | 10 | Large PAHs (pyrene, tetracene, etc.) and decylbenzene |

Elemental composition: all molecules contain only C, H, N, O. The Yalamanchi subset contains only C and H.

## Examples

### Example 1: Process a Single Molecule

```python
from deltahf.pipeline import process_molecule

result = process_molecule("C[N+](=O)[O-]", n_conformers=1, name="nitromethane")
print(f"xTB energy: {result.xtb_energy:.6f} Eh")
print(f"4-param counts: {result.atom_counts_4param}")
print(f"7-param counts: {result.atom_counts_7param}")
print(f"hybrid counts:  {result.atom_counts_hybrid}")
```

### Example 2: Fit All Models

```bash
python -m deltahf fit \
    -i deltahf/data/training_data.csv \
    --model all \
    --kfold 10 \
    --n-conformers 1 \
    --cache-dir .xtb_cache \
    -o params.json
```

This processes all 314 molecules through the pipeline (SMILES -> RDKit conformers -> xTB optimization), fits atom equivalents by least squares for each model, and reports adjusted R², RMSD, MAD, max deviation, and 10-fold CV RMSD. Standard errors on each epsilon are computed from the CV folds.

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

## Benchmarking

`benchmark.py` measures how the number of conformers affects model accuracy and wall-clock time:

```bash
python benchmark.py
```

This runs the full fitting workflow at `n_conformers` = 1, 3, 5, 10, timing each run and reporting Adj. R², RMSD, MAD, max deviation, and CV RMSD for all four models. Results are saved to `benchmark_results.csv`. xTB results are cached per `n_conformers` in `.benchmark_cache/`, so re-runs skip the expensive optimization step.

## Pipeline

For each molecule, deltahf performs the following steps:

1. **SMILES parsing** — Convert SMILES to an RDKit mol object with explicit hydrogens
2. **Conformer generation** — Generate 3D conformers via ETKDG, rank by MMFF94 energy, and prune near-duplicates by RMSD
3. **xTB optimization** — Optimize the lowest *n* unique conformers with GFN2-xTB
4. **Connectivity check** — Verify the optimized geometry hasn't isomerized (atom connectivity preserved)
5. **ΔHf° prediction** — Apply the atom equivalent formula using the lowest xTB energy

## Dependencies

- [Python](https://www.python.org/) >= 3.10
- [RDKit](https://www.rdkit.org/) — molecular representation, 3D embedding, ETKDG conformer generation
- [xTB](https://github.com/grimme-lab/xtb) — GFN2 semi-empirical optimization (CLI)
- [NumPy](https://numpy.org/) — numerical fitting
- [pandas](https://pandas.pydata.org/) — CSV I/O
- [tqdm](https://tqdm.github.io/) — progress bars

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
2. Yalamanchi, K. K.; Monge-Palacios, M.; van Oudenhoven, V. C. O.; Gao, X.; Sarathy, S. M. *J. Phys. Chem. A* **2020**, *124*, 6270–6283 [**DOI:** 10.1021/acs.jpca.0c02785](https://doi.org/10.1021/acs.jpca.0c02785)

---
License: [MIT](https://opensource.org/licenses/MIT)
