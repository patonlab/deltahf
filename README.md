# deltahf

[![Python versions](https://img.shields.io/badge/python-≥3.10-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

## Introduction

**deltahf** estimates gas-phase standard heats of formation (ΔHf°) from atom equivalent energies and semi-empirical quantum chemistry. This approach was described by Cawkwell et al.<sup>1</sup> using DFTB; this implementation replaces DFTB with [GFN2-xTB](https://github.com/grimme-lab/xtb), gxTB, or other MLIPs, requiring the atom equivalent energies to be re-parameterized.

The core equation is:

```
ΔHf° = u_optimizer − Σ(nl × εl)
```

where `u_optimizer` is the total energy from the chosen geometry optimizer (xTB, gxtb, or MLIP) in kcal/mol, `nl` is the count of atom type `l`, and `εl` is the fitted atom equivalent energy.

**Optional gxtb energies:** The `--use-gxtb` flag enables gxtb single-point energies. **CRITICAL:** gxtb and xTB energies are on completely different scales and **must never be mixed**. If you fit with `--use-gxtb`, you **must** also predict with `--use-gxtb`. See [GXTB_USAGE.md](GXTB_USAGE.md) for details. The same note of caution applies to UMA energies

### Atom Equivalent Models

Eight atom classification schemes are available, in order of increasing complexity:

| Model | Max params | Classification |
|-------|:---:|----------------|
| `element` | 7 | Elemental stoichiometry: one parameter per element type |
| `element_bo` | 11 | Adds multiply-bonded variants: C′, N′, O′, S′ (bond order > 1.25) |
| `hybrid` | 15 | RDKit hybridization: each element split by sp3, sp2, sp |
| `bondorder` | 15 | Max bond order in Kekulized structure: `_1`, `_2`, `_3` suffix |
| `bondorder_ar` | 19 | As `bondorder` but without Kekulization, preserving aromatic bonds as `_ar` |
| `extended` | 19 | Hybridization + H-count for carbon (e.g. `C_sp3_3H`, `C_sp2_1H`) |
| `bondorder_ext` | 21 | Bond order + H-count for carbon — bond-order analogue of `extended` (**best performer**) |
| `neighbour` | 31 | N/O split by hybridization × highest-priority heavy-atom neighbour (O > N > C); e.g. `N_sp2_O` for nitro N |

*Max params assumes the full C,H,N,O,F,S,Cl element set. On CHNO-only data: element=4, element\_bo=7, hybrid=10, bondorder=10, bondorder\_ar=13, extended=15, bondorder\_ext=16, neighbour=27. Model names reflect the original CHNO parameter counts.* Parameters with zero training examples are automatically excluded from fitting, so the actual fitted count is often lower than the maximum.

### Bond-Increment Model (Not Recommended)

A bond-increment scheme was investigated as an alternative, using Kekulized explicit-hydrogen bond counts (19 bond types: C–C, C–H, C–N, C–O, N–N, N–O, etc.) as descriptors instead of atom types. It performed very poorly (RMSD ~85 kcal/mol, CV RMSD ~389 kcal/mol, R² = −4.4) and was not retained.

The failure is fundamental, not incidental: bond counts are approximately linearly dependent on atom counts (e.g. `#H = #(C–H) + #(H–N) + #(H–O) + 2·#(H–H)`), making the design matrix near-singular. Least squares finds wildly large, cancelling coefficients that overfit the training set and generalise poorly. More deeply, the xTB energy decomposes naturally into per-atom (not per-bond) contributions, so atom equivalents are the correct functional form for this approach.


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

deltahf requires elemental parameters, specified in json format, to generate its predictions. If no json file is specified then default parameters generated on the training data (described below) will be used. 

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
| `--model` | Which model(s) to fit: `element`, `element_bo`, `hybrid`, `bondorder`, `bondorder_ext`, `bondorder_ar`, `extended`, `neighbour`, `both`, or `all`. | `both` |
| `--kfold` | Number of cross-validation folds. | `10` |
| `--n-conformers` | Number of lowest-energy RDKit conformers to optimize with xTB. | `1` |
| `--output, -o` | Output JSON file for fitted epsilon values. | — |
| `--csv` | Output CSV with training data and per-molecule predictions. | — |
| `--outliers [N]` | Print the top N outliers by \|error\| after fitting each model (default N=10 if flag is given without argument). | — |
| `--optimizer` | Geometry optimizer: `xtb`, `uma`, `esen`, or `aimnet2`. Non-xTB optimizers require optional dependencies. | `xtb` |
| `--xtb-threads N` | Number of OpenMP threads for xTB (4–16 is optimal; beyond 16 overhead dominates). | all CPUs |
| `--use-xtb-wbos` | Use xTB Wiberg bond orders (instead of RDKit) for `element_bo` classification. | — |
| `--use-gxtb` | Use gxtb energies (wB97M-V/def2-TZVPPD). **WARNING:** Must use consistently for fit AND predict! See [GXTB_USAGE.md](GXTB_USAGE.md). | — |
| `--cache-dir` | Directory for caching results (automatically suffixed with `_xtb` or `_gxtb`). | — |
| `--verbose, -v` | Print per-molecule details instead of a progress bar. | — |

### Subcommand: `predict`

Predict ΔHf° for new molecules using previously fitted atom equivalent energies.

```bash
python -m deltahf predict -i molecules.csv --epsilon params.json --model element --n-conformers 1 -o results.csv
```

| Option | Description | Default |
|--------|-------------|---------|
| `--input, -i` | CSV file with a `smiles` column (required). | — |
| `--epsilon` | JSON file with fitted atom equivalent energies (uses defaults from `params/` if not specified). | Default params |
| `--model` | Which model to use: `element`, `element_bo`, `hybrid`, `bondorder`, `bondorder_ext`, `bondorder_ar`, `extended`, or `neighbour`. | `element` |
| `--n-conformers` | Number of conformers to optimize with xTB. | `1` |
| `--output, -o` | Output CSV with predicted ΔHf° values. | — |
| `--optimizer` | Geometry optimizer: `xtb`, `uma`, `esen`, or `aimnet2`. | `xtb` |
| `--xtb-threads N` | Number of OpenMP threads for xTB. | all CPUs |
| `--use-gxtb` | Use gxtb energies. **Must match the method used in fit!** Automatically validated against parameter file metadata. | — |
| `--cache-dir` | Directory for caching results (automatically suffixed with `_xtb` or `_gxtb`). | — |
| `--verbose, -v` | Print per-molecule details instead of a progress bar. | — |

## Training Data

`deltahf/data/training_data.csv` contains **531 molecules** with experimental ΔHf° values covering elements C, H, N, O, F, S, Cl:

| Source | Count | Description |
|--------|-------|-------------|
| Cawkwell et al. (2021)<sup>1</sup> | 102 | Energetic CHNO molecules + small reference compounds |
| Yalamanchi et al. (2020)<sup>2</sup> | 211 | Cyclic hydrocarbons (CH only) |
| ATcT v1.220<sup>3</sup> | 218 | Gas-phase ΔHf° at 298.15 K; neutral closed-shell molecules with elements C, H, N, O, F, S, Cl |

Each molecule is assigned a `category` label:

| Category | Count | Description |
|----------|-------|-------------|
| `cyclic_HC` | 189 | Cyclic hydrocarbons (cyclopentanes through naphthalenes) |
| `small_CHNO` | 136 | Small reference molecules (methane, water, ethanol, etc.) |
| `energetic` | 45 | Explosives, nitro/azide/nitroso compounds |
| `hydrocarbon` | 41 | Acyclic hydrocarbons |
| `chlorinated` | 42 | Molecules containing Cl |
| `fluorinated` | 39 | Molecules containing F |
| `sulfur` | 17 | Molecules containing S |
| `strained_3ring` | 12 | Cyclopropane derivatives and quadricyclane |
| `large_HC` | 10 | Large PAHs (pyrene, tetracene, etc.) |

Category labels (`cyclic_HC`, `small_CHNO`, `energetic`, `hydrocarbon`, `chlorinated`, `fluorinated`, `sulfur`, `strained_3ring`, `large_HC`) are mutually exclusive in `training_data.csv` and sum to the total 531 molecules. Parameter-file metadata (`_metadata.n_molecules`) records the count of molecules successfully processed during fitting, which may be slightly lower than 531 if any molecules fail the pipeline.

The ATcT data was extracted via `atct/extract_atct.py`; see `atct/USING_ATCT_DATA.md` for the extraction workflow and pipeline statistics.

## Examples

### Example 1: Process a Single Molecule

```python
from deltahf.pipeline import process_molecule

result = process_molecule("C[N+](=O)[O-]", n_conformers=1, name="nitromethane")
print(f"xTB energy: {result.xtb_energy:.6f} Eh")
print(f"element counts:    {result.atom_counts_element}")
print(f"element_bo counts: {result.atom_counts_element_bo}")
print(f"hybrid counts:     {result.atom_counts_hybrid}")
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

This processes all 531 molecules through the pipeline (SMILES -> RDKit conformers -> xTB optimization), fits atom equivalents by least squares for each model, and reports adjusted R², RMSD, MAD, max deviation, and 10-fold CV RMSD. Standard errors on each epsilon are computed from the CV folds.

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
    --model element \
    --n-conformers 3 \
    -o predictions.csv
```

See the [examples/](examples/) directory for additional scripts covering batch prediction, model comparison, and atom classification.

## Benchmarking

### Running Benchmarks

`benchmark.py` measures accuracy across all eight models, three model chemistries (xTB, gXTB, UMA), and three conformer counts (1, 3, 5). Results are reported for both the full training set and the 102-molecule Cawkwell2021 subset (enabling comparison with the published DFT-B baseline<sup>1</sup>).

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

Results are cached per method × n_conformers in `.benchmark_cache/`, so re-runs skip expensive optimizations. See [`model_benchmarks/`](model_benchmarks/) for detailed results and the full CSV output.

### Key Findings

A benchmark across all eight models and three methods (n_conformers = 1, 531 molecules) reveals four main findings:

1. **Bond-order classification outperforms hybridization** — Using maximum bond order (1/2/3) from the Kekulized structure instead of RDKit hybridization labels improves accuracy at both the coarse level (`bondorder` vs `hybrid`) and the fine-grained level (`bondorder_ext` vs `extended`). The best model overall is `bondorder_ext`, combining bond-order labels with per-carbon H-counts.

2. **Model chemistry matters far more than parameterisation** — Upgrading from xTB to gXTB (wB97M-V/def2-TZVPPD single-points on xTB geometries) reduces RMSD by ~42% (6.53 → 3.82 kcal/mol for `bondorder_ext`). Upgrading to UMA (MLIP, GPU) reduces it by a further ~32% (3.82 → 2.59 kcal/mol).

3. **Number of conformers has minimal impact for xTB and gXTB** — Increasing from n=1 to n=5 gives negligible accuracy change at 3–4× the cost. For UMA, n=5 actually *degrades* accuracy (RMSD 2.59 → 3.11 for `bondorder_ext`), likely because UMA's more flexible potential finds lower-energy geometries that have isomerized. Use `--n-conformers 1` (the default).

4. **xTB + bondorder_ext matches or exceeds the published DFT-B baseline** — On the 102-molecule Cawkwell2021 subset, xTB + `bondorder_ext` (RMSD = 6.91 kcal/mol) approaches the DFT-B + `element_bo` result (RMSD = 6.08 kcal/mol) from Cawkwell et al.<sup>1</sup>, while gXTB and UMA substantially surpass it.

5. **The physics prior from xTB is essential** — A pure cheminformatics baseline (Morgan fingerprint + Random Forest) achieved an in-sample RMSD of 11.70 kcal/mol but a cross-validated RMSD of 30.9 kcal/mol — dramatically worse than even the simplest xTB model (`element`, CV RMSD = 11.9 kcal/mol). The atom equivalent approach sidesteps this by encoding the quantum mechanical energy decomposition directly.

> **Note on the `neighbour` model:** The `neighbour` model shows competitive training-set RMSD but produces an extremely large cross-validation RMSD in some CV folds (near-singular design matrix). It is not recommended for practical use.

### Results: Effect of Model Parameterisation

Full training set (531 molecules, C,H,N,O,F,S,Cl), n_conformers = 1:

| Model | Params | xTB RMSD | xTB MAD | gXTB RMSD | gXTB MAD | UMA RMSD | UMA MAD |
|-------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| `element` | 7 | 11.54 | 8.16 | 4.37 | 3.18 | 3.27 | 2.32 |
| `element_bo` | 11 | 9.12 | 5.78 | 4.07 | 2.81 | 2.99 | 1.98 |
| `hybrid` | 13 | 8.92 | 5.50 | 4.11 | 2.88 | 2.78 | 1.74 |
| `bondorder` | 13 | 7.63 | 5.00 | 3.95 | 2.75 | 2.73 | 1.71 |
| `bondorder_ar` | 16 | 7.60 | 4.96 | 3.94 | 2.74 | 2.71 | 1.67 |
| `extended` | 18 | 8.44 | 5.37 | 3.98 | 2.76 | 2.69 | 1.69 |
| **`bondorder_ext`** | **19** | **6.53** | **4.26** | **3.82** | **2.67** | **2.59** | **1.66** |

All values in kcal/mol. Params = active parameters after dropping zero-count types. Adj. R² and CV RMSD available in [`model_benchmarks/benchmark_results.csv`](model_benchmarks/benchmark_results.csv).

### Results: Comparison with DFT-B Literature Baseline

Cawkwell2021 subset (102 molecules), n_conformers = 1. The DFT-B baseline<sup>1</sup> used DFTB+ geometries and energies; our methods use xTB/gXTB/UMA geometries with re-fitted atom equivalents.

| Method | Model | Params | RMSD (kcal/mol) | Max Dev (kcal/mol) |
|--------|-------|:---:|:---:|:---:|
| DFT-B (lit.)<sup>1</sup> | `element` | 4 | 7.59 | 25.48 |
| DFT-B (lit.)<sup>1</sup> | `element_bo` | 7 | 6.08 | 15.01 |
| xTB | `bondorder_ext` | 15 | 6.91 | 23.66 |
| gXTB | `bondorder_ext` | 15 | 3.75 | 10.52 |
| **UMA** | **`bondorder_ext`** | **15** | **2.70** | **8.95** |

### Results: Effect of n_conformers

`bondorder_ext` model, full training set (526 molecules). Timings are from fresh runs (uncached).

| n_conformers | xTB RMSD | gXTB RMSD | UMA RMSD | xTB time | gXTB time | UMA time |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 6.53 | 3.82 | 2.59 | 54 s | 69 s | 285 s |
| 3 | 6.51 | 3.81 | 2.58 | 104 s | 119 s | 707 s |
| 5 | 6.52 | 3.83 | 3.11\* | 151 s | 166 s | 1095 s |

\*UMA n=5 degrades due to isomerization of some conformers during UMA geometry optimization (max deviation jumps from 20 to 41 kcal/mol).

### Recommendations

- **Use `--n-conformers 1`** (the default) — additional conformers give negligible improvement for xTB/gXTB and can degrade UMA results
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
- [gxtb](https://github.com/grimme-lab/gxtb) — **Must be installed manually from source** (not available via pip/conda). **WARNING:** gxtb and xTB energies are on completely different scales and cannot be mixed - see [GXTB_USAGE.md](GXTB_USAGE.md) for proper usage.

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
3. Ruscic, B.; Bross, D. H. Active Thermochemical Tables (ATcT) values based on ver. 1.220 of the Thermochemical Network. Argonne National Laboratory, 2023. Available at https://atct.anl.gov/

---
License: [MIT](https://opensource.org/licenses/MIT)
