# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python package estimating gas-phase standard heats of formation (ΔHf°) from atom equivalent energies and GFN2-xTB semi-empirical quantum chemistry, based on the approach by Cawkwell et al. (J. Chem. Inf. Model. 2021, 61, 3337–3347).

Core equation: `ΔHf° = u_xtb − Σ(nl × εl)` where `u_xtb` is xTB total energy (kcal/mol), `nl` is atom count for type `l`, and `εl` is the fitted atom equivalent energy.

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

# CLI usage
python -m deltahf fit -i deltahf/data/training_data.csv --model all --kfold 10 --n-conformers 1 -o params.json
python -m deltahf predict -i molecules.csv --epsilon params.json --model 4param --n-conformers 5 -o results.csv
```

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
