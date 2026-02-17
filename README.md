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

**Optional gxtb energies:** The `--use-gxtb` flag enables gxtb (wB97M-V/def2-TZVPPD) single-point energies. **CRITICAL:** gxtb and xTB energies are on completely different scales (~10-17x difference) and **must never be mixed**. If you fit with `--use-gxtb`, you **must** also predict with `--use-gxtb`. See [GXTB_USAGE.md](GXTB_USAGE.md) for details.

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
python -m deltahf fit -i training.csv --model all --kfold 10 --n-conformers 1 -o params.json
```

| Option | Description | Default |
|--------|-------------|---------|
| `--input, -i` | CSV file with `smiles` and `exp_dhf_kcal_mol` columns (required). | — |
| `--model` | Which model(s) to fit: `4param`, `7param`, `hybrid`, `extended`, `both`, or `all`. | `both` |
| `--kfold` | Number of cross-validation folds. | `10` |
| `--n-conformers` | Number of lowest-energy RDKit conformers to optimize with xTB. | `1` |
| `--output, -o` | Output JSON file for fitted epsilon values. | — |
| `--csv` | Output CSV with training data and per-molecule predictions. | — |
| `--use-xtb-wbos` | Use xTB Wiberg bond orders (instead of RDKit) for 7-param classification. | — |
| `--use-gxtb` | Use gxtb energies (wB97M-V/def2-TZVPPD). **WARNING:** Must use consistently for fit AND predict! See [GXTB_USAGE.md](GXTB_USAGE.md). | — |
| `--cache-dir` | Directory for caching results (automatically suffixed with `_xtb` or `_gxtb`). | — |
| `--verbose, -v` | Print per-molecule details instead of a progress bar. | — |

### Subcommand: `predict`

Predict ΔHf° for new molecules using previously fitted atom equivalent energies.

```bash
python -m deltahf predict -i molecules.csv --epsilon params.json --model 4param --n-conformers 1 -o results.csv
```

| Option | Description | Default |
|--------|-------------|---------|
| `--input, -i` | CSV file with a `smiles` column (required). | — |
| `--epsilon` | JSON file with fitted atom equivalent energies (uses defaults from `params/` if not specified). | Default params |
| `--model` | Which model to use: `4param`, `7param`, `hybrid`, or `extended`. | `4param` |
| `--n-conformers` | Number of conformers to optimize with xTB. | `1` |
| `--output, -o` | Output CSV with predicted ΔHf° values. | — |
| `--use-gxtb` | Use gxtb energies. **Must match the method used in fit!** Automatically validated against parameter file metadata. | — |
| `--cache-dir` | Directory for caching results (automatically suffixed with `_xtb` or `_gxtb`). | — |
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

### Running Benchmarks

`benchmark.py` measures how the number of conformers affects model accuracy and wall-clock time:

```bash
# Benchmark with xTB energies
python benchmark.py

# Benchmark with gxtb energies (requires gxtb binary)
python benchmark.py --use-gxtb
```

This runs the full fitting workflow at `n_conformers` = 1, 3, 5, timing each run and reporting Adj. R², RMSD, MAD, max deviation, and CV RMSD for all four models. xTB results are cached per `n_conformers` in `.benchmark_cache/`, so re-runs skip the expensive optimization step.

### Key Findings

Benchmarks on the 314-molecule training set (n_conformers = 1, 3, 5) reveal three important insights:

1. **Number of conformers has minimal impact on accuracy** — Increasing from 1 to 5 conformers provides essentially no improvement in predictive accuracy, while increasing computational cost ~3–4×.

2. **Higher-level energies significantly improve accuracy** — Both gxtb (wB97M-V/def2-TZVPPD single-points on xTB geometries) and UMA (MLIP optimization on GPU) reduce RMSD by ~50–65% compared to xTB, with UMA achieving the best absolute accuracy.

3. **UMA (MLIP) gives the best accuracy** — UMA optimization yields RMSD ~2.6 kcal/mol vs ~3.4 kcal/mol for gxtb and ~7.4 kcal/mol for xTB (extended model, n=1), at the cost of requiring a GPU.

### Results Summary

Times are shown relative to n_conformers = 1 for each method.

**xTB energies (GFN2-xTB optimization, CPU):**

| n_conformers | Model | n_params | Adj. R² | RMSD (kcal/mol) | MAD (kcal/mol) | Max Dev (kcal/mol) | CV RMSD | Rel. Time |
|:---:|-------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 4param | 4 | 0.912 | 11.09 | 7.10 | 76.10 | 11.46 | 1.0× |
| 1 | 7param | 7 | 0.950 | 8.33 | 4.56 | 71.03 | 8.87 | 1.0× |
| 1 | hybrid | 9 | 0.959 | 7.49 | 4.13 | 49.31 | 8.26 | 1.0× |
| 1 | extended | 14 | **0.960** | **7.36** | **3.95** | **49.11** | **8.22** | 1.0× |
| 3 | extended | 14 | 0.960 | 7.38 | 3.96 | 49.12 | 8.25 | 2.2× |
| 5 | extended | 14 | 0.959 | 7.46 | 3.98 | 52.13 | 8.35 | 3.5× |

**gxtb energies (xTB geometry + wB97M-V/def2-TZVPPD single-point, CPU):**

| n_conformers | Model | n_params | Adj. R² | RMSD (kcal/mol) | MAD (kcal/mol) | Max Dev (kcal/mol) | CV RMSD | Rel. Time |
|:---:|-------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 4param | 4 | 0.989 | 4.01 | 3.02 | 15.00 | 4.16 | 1.0× |
| 1 | 7param | 7 | **0.992** | **3.40** | **2.43** | **13.59** | **3.61** | 1.0× |
| 1 | hybrid | 9 | 0.991 | 3.54 | 2.47 | 13.78 | 3.84 | 1.0× |
| 1 | extended | 14 | 0.991 | 3.43 | 2.43 | 14.14 | 3.81 | 1.0× |
| 3 | extended | 14 | 0.991 | 3.43 | 2.42 | 14.13 | 3.81 | 2.1× |
| 5 | extended | 14 | 0.991 | 3.43 | 2.42 | 14.34 | 3.81 | 3.1× |

**UMA energies (MLIP optimization, GPU):**

| n_conformers | Model | n_params | Adj. R² | RMSD (kcal/mol) | MAD (kcal/mol) | Max Dev (kcal/mol) | CV RMSD | Rel. Time |
|:---:|-------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 4param | 4 | 0.993 | 3.11 | 2.24 | 12.04 | 3.26 | 1.0× |
| 1 | 7param | 7 | 0.995 | 2.77 | 1.88 | 11.10 | 2.93 | 1.0× |
| 1 | hybrid | 9 | 0.994 | 2.78 | 1.85 | 11.66 | 2.95 | 1.0× |
| 1 | extended | 14 | **0.995** | **2.62** | **1.80** | **11.54** | **2.85** | 1.0× |
| 3 | extended | 14 | 0.995 | 2.62 | 1.80 | 11.57 | 2.86 | 2.6× |
| 5 | extended | 14 | **0.995** | **2.58** | **1.77** | **11.38** | **2.79** | 4.1× |

**Recommendations:**
- Use `--n-conformers 1` for all methods (now the default) — additional conformers give negligible accuracy gains
- Use `--optimizer uma` for best accuracy if a GPU is available
- Use `--use-gxtb` for improved accuracy over xTB without requiring a GPU
- The extended model gives the best accuracy, but 7param is nearly as good with fewer parameters

Complete benchmark results: [xtb_benchmark_results.md](xtb_benchmark_results.md) | [gxtb_benchmark_results.md](gxtb_benchmark_results.md) | [uma_benchmark_results.md](uma_benchmark_results.md)

## Pipeline

For each molecule, deltahf performs the following steps:

1. **SMILES parsing** — Convert SMILES to an RDKit mol object with explicit hydrogens
2. **Conformer generation** — Generate 3D conformers via ETKDG, rank by MMFF94 energy, and prune near-duplicates by RMSD
3. **xTB optimization** — Optimize the lowest *n* unique conformers with GFN2-xTB
4. **Connectivity check** — Verify the optimized geometry hasn't isomerized (atom connectivity preserved)
5. **gxtb refinement** (optional) — Compute single-point energy with gxtb on the lowest-energy xTB conformer
6. **ΔHf° prediction** — Apply the atom equivalent formula using the lowest xTB (or gxtb) energy

## Dependencies

**Required:**
- [Python](https://www.python.org/) >= 3.10
- [RDKit](https://www.rdkit.org/) — molecular representation, 3D embedding, ETKDG conformer generation
- [xTB](https://github.com/grimme-lab/xtb) — GFN2 semi-empirical optimization (CLI, available via conda-forge)
- [NumPy](https://numpy.org/) — numerical fitting
- [pandas](https://pandas.pydata.org/) — CSV I/O
- [tqdm](https://tqdm.github.io/) — progress bars

**Optional:**
- [gxtb](https://github.com/grimme-lab/gxtb) — Higher-level DFT energies (wB97M-V/def2-TZVPPD). **Must be installed manually from source** (not available via pip/conda). **WARNING:** gxtb and xTB energies are on completely different scales and cannot be mixed - see [GXTB_USAGE.md](GXTB_USAGE.md) for proper usage.

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
