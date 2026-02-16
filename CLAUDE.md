# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python package estimating gas-phase standard heats of formation (ΔHf°) from atom equivalent energies and GFN2-xTB semi-empirical quantum chemistry, based on the approach by Cawkwell et al. (J. Chem. Inf. Model. 2021, 61, 3337–3347).

Core equation: `ΔHf° = u_xtb − Σ(nl × εl)` where `u_xtb` is xTB total energy (kcal/mol), `nl` is atom count for type `l`, and `εl` is the fitted atom equivalent energy.

## Dependencies

Core runtime dependencies: RDKit, NumPy, pandas. Development: pytest, ruff.

**External CLI dependency:** [xTB](https://github.com/grimme-lab/xtb) must be on PATH for integration tests, `fit`, and `predict` commands. Install via conda: `conda install -c conda-forge xtb`.

**Optional CLI dependency:** gxtb for refined single-point energies after xTB optimization. Requires manual installation from source (see gxtb documentation). Enable with `--use-gxtb` flag.

**IMPORTANT:** gxtb and xTB energies are on completely different scales (~10-17x difference) and must NEVER be mixed. Use `--use-gxtb` consistently for both fitting AND prediction, or use only xTB (default). See `GXTB_USAGE.md` for details.

## Commands

```bash
# Environment setup
conda env create -f environment.yml
conda activate deltahf
pip install -e ".[dev]"

# Run all unit tests (no xTB required)
pytest -m "not integration"

# Run integration tests (requires xTB on PATH)
pytest -m integration

# Run a single test
pytest tests/test_smiles.py::TestCountAtoms::test_methane

# Lint
ruff check deltahf/ tests/

# Benchmarking
python benchmark.py  # Measures accuracy vs n_conformers (1, 3, 5, 10), outputs benchmark_results.csv
python benchmark.py --use-gxtb  # Same benchmark using gxtb energies

# CLI usage
python -m deltahf fit -i deltahf/data/training_data.csv --model all --kfold 10 --n-conformers 1 -o params.json
python -m deltahf predict -i molecules.csv --epsilon params.json --model 4param --n-conformers 1 -o results.csv
```

### Common CLI Flags

- `--model`: `4param`, `7param`, `hybrid`, `extended`, `both` (4+7), or `all` (for `fit`); single model for `predict`
- `--n-conformers N`: Number of ETKDG conformers to optimize with xTB (default: 1)
- `--cache-dir PATH`: Cache results as JSON for restart capability. Automatically suffixed with `_xtb` or `_gxtb` based on method.
- `--use-xtb-wbos`: Use xTB Wiberg bond orders for 7-param classification instead of RDKit (disables caching)
- `--use-gxtb`: Use gxtb (wB97M-V/def2-TZVPPD) energies instead of xTB. Must be used consistently for both fit AND predict! Cache directory automatically separated.
- `--verbose, -v`: Print per-molecule details instead of progress bar
- `--csv FILE` (fit only): Output CSV with training data + predictions + errors for each model
- `-o, --output FILE`: JSON output for fitted parameters (fit) or CSV for predictions (predict). Parameters include `_metadata` tracking the method used.

## Architecture

### Atom Equivalent Models

Four atom classification schemes, each producing a different set of parameters:

| Model | Max Params | Classification |
| ----- | --------- | -------------- |
| `4param` | 4 | Elemental stoichiometry: C, H, N, O |
| `7param` | 7 | Adds multiply-bonded variants: C', N', O' (bond order > 1.25) |
| `hybrid` | 10 | RDKit hybridization: C/N/O split by sp3, sp2, sp |
| `extended` | 15 | Hybridization + H-count for carbon (e.g. C_sp3_3H, C_sp2_1H) |

Parameters with zero training examples are automatically excluded from fitting, so the actual parameter count may be less than the max (e.g., if no O_sp atoms exist in the training data, the hybrid model fits 9 parameters instead of 10).

### Pipeline Flow

Per molecule: SMILES → RDKit mol → ETKDG conformers → MMFF ranking → xTB optimization of lowest n → connectivity check (detect isomerization) → atom equivalent ΔHf prediction.

### Key Modules

- `deltahf/smiles.py` — SMILES→mol, `count_atoms()` (4-param), `classify_atoms_7param()`, `classify_atoms_hybrid()`, `classify_atoms_extended()`
- `deltahf/conformers.py` — ETKDG conformer generation, MMFF energy ranking, XYZ file writing
- `deltahf/xtb.py` — xTB CLI subprocess wrapper, energy parsing (regex on stdout), Wiberg bond order parsing
- `deltahf/atom_equivalents.py` — least-squares fitting (`np.linalg.lstsq`), k-fold CV with RMSD and epsilon standard errors, adjusted R², RMSD/MAD/max deviation statistics
- `deltahf/pipeline.py` — end-to-end orchestration: `process_molecule()`, `process_csv()`
- `deltahf/cache.py` — JSON-based result caching for xTB runs (enables restart capability)
- `deltahf/__main__.py` — CLI with `fit` and `predict` subcommands
- `deltahf/data/` — 314-molecule training dataset (CSV) and `load_training_data()` helper

### Constants

`HARTREE_TO_KCAL = 627.5094740631` defined in both `xtb.py` and `atom_equivalents.py`.

### Bond Order Classification (7-param)

Uses RDKit `GetBondTypeAsDouble()` from the molecular graph. An atom is "primed" (multiply-bonded) if it has any bond > 1.25 OR has nonzero formal charge with a neighbor that has a bond > 1.25 (handles charge-separated resonance forms like `[N+](=O)[O-]`). Optional `--use-xtb-wbos` flag uses xTB Wiberg bond orders instead.

### Fitting Details

- Dynamic parameter filtering: `_print_model_results()` checks total atom counts and drops parameters with zero examples
- Adjusted R²: `R²_adj = 1 - (1 - R²) × (n-1)/(n-p-1)` where p = number of parameters
- CV reports RMSD (not MSD) and standard errors on each epsilon from fold-to-fold variation

## Training Data

`deltahf/data/training_data.csv` contains **314 CHNO molecules** with columns: `id, name, formula, smiles, exp_dhf_kcal_mol, source, category`.

Sources: **Cawkwell2021** (102 CHNO molecules) and **Yalamanchi2020** (212 CH cyclic hydrocarbons, enthalpies converted from kJ/mol).

Categories: `energetic` (45), `small_CHNO` (57), `cyclic_HC` (189), `strained_3ring` (13), `large_HC` (10). All molecules contain only C, H, N, O.

## Benchmarking

`benchmark.py` (in the repo root) measures how `n_conformers` affects accuracy and runtime:

- Tests with `n_conformers` = 1, 3, 5, 10
- Runs full fitting workflow with 10-fold CV for all four models
- Reports Adj. R², RMSD, MAD, max deviation, CV RMSD, and wall time for each configuration
- Outputs `benchmark_results.csv` with comparative results
- Uses separate cache directories per `n_conformers` (`.benchmark_cache/n1/`, etc.) to enable restart capability

## Cache Mechanism

Both `fit` and `predict` commands support `--cache-dir` to cache xTB results:

- Each molecule's xTB results (energy, coordinates, Wiberg bond orders) are saved as JSON
- Enables restarting failed runs without re-running expensive xTB calculations
- Cache key is based on SMILES string and `n_conformers`
- Used extensively in `benchmark.py` to avoid redundant xTB calls
- **Important:** Caching is automatically disabled when `--use-xtb-wbos` is set

## Output Formats

**fit --csv**: Adds columns to input CSV:
- `xtb_energy_eh`, `xtb_energy_kcal_mol`: xTB total energies
- `gxtb_energy_eh`, `gxtb_energy_kcal_mol`: gxtb single-point energies (if `--use-gxtb` is set)
- `pred_dhf_4param`, `pred_dhf_7param`, etc.: Predicted ΔHf° for each model (kcal/mol, uses gxtb energy if available)
- `error_4param`, `error_7param`, etc.: Prediction error (pred - exp) for each model
- `error`: Error message for failed molecules (None if successful)

**predict -o**: Creates CSV with columns:
- All input columns (including `smiles`, `name` if present)
- `xtb_energy_kcal_mol`: xTB energy
- `gxtb_energy_kcal_mol`: gxtb single-point energy (if `--use-gxtb` is set)
- `pred_dhf_kcal_mol`: Predicted ΔHf° using specified model (uses gxtb energy if available)
- `error`: Error message (None if successful)
