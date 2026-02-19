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

Eight atom classification schemes are available, in order of increasing complexity:

| Model | Max params | Classification |
|-------|:---:|----------------|
| `4param` | 4 | Elemental stoichiometry: C, H, N, O |
| `7param` | 7 | Adds multiply-bonded variants: C′, N′, O′ (bond order > 1.25) |
| `hybrid` | 10 | RDKit hybridization: C/N/O split by sp3, sp2, sp |
| `bondorder` | 10 | Max bond order in Kekulized structure: `_1`, `_2`, `_3` suffix |
| `bondorder_ar` | 13 | As `bondorder` but without Kekulization, preserving aromatic bonds as `_ar` |
| `extended` | 15 | Hybridization + H-count for carbon (e.g. `C_sp3_3H`, `C_sp2_1H`) |
| `bondorder_ext` | 16 | Bond order + H-count for carbon — bond-order analogue of `extended` (**best performer**) |
| `neighbour` | 27 | N/O split by hybridization × highest-priority heavy-atom neighbour (O > N > C); e.g. `N_sp2_O` for nitro N |

Parameters with zero training examples are automatically excluded from fitting, so the actual parameter count is often lower than the maximum.

### Bond-Increment Model (Not Recommended)

A bond-increment scheme was investigated as an alternative, using Kekulized explicit-hydrogen bond counts (19 bond types: C–C, C–H, C–N, C–O, N–N, N–O, etc.) as descriptors instead of atom types. It performed very poorly (RMSD ~85 kcal/mol, CV RMSD ~389 kcal/mol, R² = −4.4) and was not retained.

The failure is fundamental, not incidental: bond counts are approximately linearly dependent on atom counts (e.g. `#H = #(C–H) + #(H–N) + #(H–O) + 2·#(H–H)`), making the design matrix near-singular. Least squares finds wildly large, cancelling coefficients that overfit the training set and generalise poorly. More deeply, the xTB energy decomposes naturally into per-atom (not per-bond) contributions, so atom equivalents are the correct functional form for this approach.

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
| `--model` | Which model(s) to fit: `4param`, `7param`, `hybrid`, `bondorder`, `bondorder_ext`, `bondorder_ar`, `extended`, `neighbour`, `both`, or `all`. | `both` |
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
| `--model` | Which model to use: `4param`, `7param`, `hybrid`, `bondorder`, `bondorder_ext`, `bondorder_ar`, `extended`, or `neighbour`. | `4param` |
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

`benchmark.py` measures accuracy across all eight models, three model chemistries (xTB, gXTB, UMA), and three conformer counts (1, 3, 5). Results are reported for both the full 314-molecule training set and the 102-molecule Cawkwell2021 subset (enabling comparison with the published DFT-B baseline<sup>1</sup>).

```bash
# xTB (CPU)
python benchmark.py --methods xtb

# gXTB (CPU, requires gxtb binary) — append to existing CSV
python benchmark.py --methods gxtb --append

# UMA (GPU, requires fairchem-core and MODEL_DIR) — append to existing CSV
python benchmark.py --methods uma --append

# Pure cheminformatics RF baseline (Morgan FP + Random Forest, no quantum chemistry)
python benchmark.py --methods rf --append  # requires: pip install molpipeline
```

Results are cached per method × n_conformers in `.benchmark_cache/`, so re-runs skip expensive optimizations. See `benchmark_results.csv` for the full output.

### Key Findings

A comprehensive benchmark (314 molecules, n_conformers = 1, 3, 5) reveals four main findings:

1. **Bond-order classification outperforms hybridization** — Using maximum bond order (1/2/3) from the Kekulized structure instead of RDKit hybridization labels improves accuracy at both the coarse level (`bondorder` vs `hybrid`, same 10 params) and the fine-grained level (`bondorder_ext` vs `extended`). The best model overall is `bondorder_ext` (16 params), combining bond-order labels with per-carbon H-counts.

2. **Model chemistry matters far more than parameterisation** — Upgrading from xTB to gXTB (wB97M-V/def2-TZVPPD single-points on xTB geometries) reduces RMSD by ~33% (4.71 → 3.15 kcal/mol for `bondorder_ext`). Upgrading to UMA (MLIP, GPU) reduces it by a further ~21% (3.15 → 2.47 kcal/mol).

3. **Number of conformers has minimal impact** — Increasing from n=1 to n=5 gives negligible accuracy gains (<1% RMSD change) at a cost of 3–4× more computation. Use `--n-conformers 1` (the default).

4. **xTB + bondorder_ext matches or exceeds the published DFT-B baseline** — On the 102-molecule Cawkwell2021 subset, xTB + `bondorder_ext` (RMSD = 6.91 kcal/mol) approaches the DFT-B + 7param result (RMSD = 6.08 kcal/mol) from Cawkwell et al.<sup>1</sup>, while gXTB and UMA substantially surpass it.

5. **The physics prior from xTB is essential** — A pure cheminformatics baseline (Morgan fingerprint + Random Forest, no quantum chemistry) achieved an in-sample RMSD of 8.89 kcal/mol but a cross-validated RMSD of 22.9 kcal/mol on the full training set — dramatically worse than even the simplest xTB model (4param, CV RMSD = 11.5 kcal/mol). The severe overfitting reflects the small dataset size (314 molecules, 2048-dimensional fingerprint). The atom equivalent approach sidesteps this by encoding the quantum mechanical energy decomposition directly, requiring only 4–16 scalar parameters.

> **Note on the `neighbour` model:** The `neighbour` model (27 params) shows competitive training-set RMSD but produces an extremely large cross-validation RMSD (hundreds of kcal/mol), indicating instability with near-singular design matrix in some CV folds. It is not recommended for practical use.

### Results: Effect of Model Parameterisation

Full dataset (311 molecules), n_conformers = 1:

| Model | Params | xTB RMSD | xTB MAD | gXTB RMSD | gXTB MAD | UMA RMSD | UMA MAD |
|-------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `4param` | 4 | 11.09 | 7.10 | 4.01 | 3.02 | 3.11 | 2.24 |
| `7param` | 7 | 8.33 | 4.56 | 3.40 | 2.43 | 2.77 | 1.88 |
| `hybrid` | 9 | 7.49 | 4.13 | 3.54 | 2.47 | 2.78 | 1.85 |
| `bondorder` | 9 | 5.77 | 3.65 | 3.26 | 2.36 | 2.69 | 1.80 |
| `bondorder_ar` | 12 | 5.75 | 3.61 | 3.24 | 2.36 | 2.69 | 1.79 |
| `extended` | 14 | 7.36 | 3.95 | 3.43 | 2.43 | 2.62 | 1.80 |
| **`bondorder_ext`** | **15** | **4.71** | **3.01** | **3.15** | **2.31** | **2.47** | **1.73** |

All values in kcal/mol. Adj. R² and CV RMSD available in `benchmark_results.csv`.

### Results: Comparison with DFT-B Literature Baseline

Cawkwell2021 subset (102 molecules), n_conformers = 1. The DFT-B baseline<sup>1</sup> used DFTB+ geometries and energies; our methods use xTB/gXTB/UMA geometries with re-fitted atom equivalents.

| Method | Model | Params | RMSD (kcal/mol) | Max Dev (kcal/mol) |
|--------|-------|:---:|:---:|:---:|
| DFT-B (lit.)<sup>1</sup> | `4param` | 4 | 7.59 | 25.48 |
| DFT-B (lit.)<sup>1</sup> | `7param` | 7 | 6.08 | 15.01 |
| xTB | `bondorder_ext` | 15 | 6.91 | 23.66 |
| gXTB | `bondorder_ext` | 15 | 3.74 | 10.52 |
| **UMA** | **`bondorder_ext`** | **15** | **2.70** | **8.95** |

### Results: Effect of n_conformers

`bondorder_ext` model, full dataset. gXTB timing is from fresh runs (uncached); xTB and UMA timings reflect cache I/O only and are not directly comparable.

| n_conformers | xTB RMSD | gXTB RMSD | UMA RMSD | gXTB wall time |
|:---:|:---:|:---:|:---:|:---:|
| 1 | 4.71 | 3.15 | 2.47 | 78 s |
| 3 | 4.73 | 3.15 | 2.47 | 163 s |
| 5 | 4.77 | 3.15 | 2.42 | 245 s |

### Recommendations

- **Use `--n-conformers 1`** (the default) — additional conformers give <1% RMSD improvement at 3–4× cost
- **Use `--model bondorder_ext`** for best accuracy
- **Use `--optimizer uma`** for best accuracy if a GPU is available
- **Use `--use-gxtb`** for substantially improved accuracy over xTB on CPU alone

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
