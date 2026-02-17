# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python package estimating gas-phase standard heats of formation (ΔHf°) from atom equivalent energies and GFN2-xTB semi-empirical quantum chemistry, based on the approach by Cawkwell et al. (J. Chem. Inf. Model. 2021, 61, 3337–3347).

Core equation: `ΔHf° = u_xtb − Σ(nl × εl)` where `u_xtb` is xTB total energy (kcal/mol), `nl` is atom count for type `l`, and `εl` is the fitted atom equivalent energy.

## Dependencies

Core runtime dependencies: RDKit, NumPy, pandas. Development: pytest, ruff.

**External CLI dependency:** [xTB](https://github.com/grimme-lab/xtb) must be on PATH for integration tests, `fit`, and `predict` commands. Install via conda: `conda install -c conda-forge xtb`.

**Optional CLI dependency:** gxtb for refined single-point energies after xTB optimization. Requires manual installation from source (see gxtb documentation). Enable with `--use-gxtb` flag.

**Optional MLIP optimizers:** Alternative geometry optimizers via ASE can replace xTB using `--optimizer uma|esen|aimnet2`:
- `uma` / `esen`: Requires `fairchem-core` (pip) and model `.pt` files. Set `MODEL_DIR` environment variable (or `.env` file) to the directory containing the model files.
- `aimnet2`: Requires `aimnet2calc` (pip).
- MLIP energies are in eV (converted internally to Hartree); they are on a completely different scale from xTB and gxtb energies. Cache directories are automatically isolated per optimizer.
- `--use-xtb-wbos` is incompatible with non-xTB optimizers.

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
python benchmark.py  # Measures accuracy vs n_conformers (1, 3, 5), outputs benchmark_results.csv
python benchmark.py --use-gxtb  # Same benchmark using gxtb energies
python benchmark.py --optimizer uma  # Benchmark with UMA optimizer

# CLI usage
python -m deltahf fit -i deltahf/data/training_data.csv --model all --kfold 10 --n-conformers 1 -o params.json
python -m deltahf predict -i molecules.csv --epsilon params.json --model 4param --n-conformers 1 -o results.csv
```

### Common CLI Flags

- `--model`: `4param`, `7param`, `hybrid`, `extended`, `both` (4+7), or `all` (for `fit`); single model for `predict`
- `--n-conformers N`: Number of ETKDG conformers to optimize (default: 1)
- `--optimizer`: `xtb` (default), `uma`, `esen`, or `aimnet2`. Selects the geometry optimizer. Non-xTB optimizers require optional dependencies (see Dependencies).
- `--cache-dir PATH`: Cache results as JSON for restart capability. Automatically suffixed with `_xtb`, `_uma`, `_esen`, etc. (plus `_gxtb` if applicable).
- `--use-xtb-wbos`: Use xTB Wiberg bond orders for 7-param classification instead of RDKit (disables caching; incompatible with non-xTB optimizers)
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

Per molecule: SMILES → RDKit mol → ETKDG conformers → MMFF ranking → optimizer (xTB or MLIP) of lowest n conformers → connectivity check (detect isomerization) → optional gxtb single-point → atom equivalent ΔHf prediction.

### Key Modules

- `deltahf/constants.py` — single source of truth for `HARTREE_TO_KCAL`, `HARTREE_TO_EV`, `PARAM_NAMES_*`, and `MODEL_DEFS`. No deltahf imports (leaf module).
- `deltahf/smiles.py` — SMILES→mol, `count_atoms()` (4-param), `classify_atoms_7param()`, `classify_atoms_hybrid()`, `classify_atoms_extended()`
- `deltahf/conformers.py` — ETKDG conformer generation, MMFF energy ranking, XYZ file writing
- `deltahf/xtb.py` — xTB CLI subprocess wrapper, energy parsing (regex on stdout), Wiberg bond order parsing
- `deltahf/uma.py` — MLIP optimizer via FAIRChem/AIMNet2 + ASE. Uses L-BFGS with `fmax=0.05` eV/Å and `DEFAULT_MAX_STEPS=250`. Loads `MODEL_DIR` from `.env` via python-dotenv. Sets `charge` and `spin` in `atoms.info` for FAIRChem. Writes plain XYZ (`format="xyz"`) for RDKit compatibility.
- `deltahf/atom_equivalents.py` — least-squares fitting (`np.linalg.lstsq`), k-fold CV with RMSD and epsilon standard errors, adjusted R², RMSD/MAD/max deviation statistics
- `deltahf/pipeline.py` — end-to-end orchestration: `process_molecule()`, `process_csv()`. `MoleculeResult` carries both `xtb_energy` and `mlip_energy` fields (`mlip_energy` covers UMA, eSEN, AIMNet2). Energy priority via `_best_energy_kcal()`: gxtb > mlip > xtb. Charge is auto-detected from SMILES via `Chem.GetFormalCharge()`.
- `deltahf/cache.py` — JSON-based result caching. Cache key is canonical SMILES; lookup validates `n_conformers`, `optimizer`, and `charge` for correctness. Old cache entries without `optimizer` field default to `"xtb"` for backwards compatibility.
- `deltahf/__main__.py` — CLI with `fit` and `predict` subcommands. MLIP model loaded once before the molecule loop (expensive); passed through to `process_molecule()`.
- `deltahf/data/` — 314-molecule training dataset (CSV) and `load_training_data()` helper

### Constants

All physical constants and model definitions live in `deltahf/constants.py` (the only module with no deltahf imports): `HARTREE_TO_KCAL = 627.5094740631`, `HARTREE_TO_EV = 27.211386245988`, `PARAM_NAMES_*`, and `MODEL_DEFS`. Other modules import from there; `atom_equivalents.py` re-exports `PARAM_NAMES_*` for backwards compatibility.

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

- Tests with `n_conformers` = 1, 3, 5
- Runs full fitting workflow with 10-fold CV for all four models
- Reports Adj. R², RMSD, MAD, max deviation, CV RMSD, and wall time for each configuration
- Outputs `benchmark_results.csv` with comparative results
- Uses separate cache directories per `n_conformers` in `.benchmark_cache/` to enable restart capability

## Cache Mechanism

Both `fit` and `predict` commands support `--cache-dir` to cache optimizer results:

- Each molecule's results are saved as JSON keyed by canonical SMILES
- Lookup validates `n_conformers`, `optimizer`, and `charge` — mismatches result in a cache miss
- Enables restarting failed runs without re-running expensive optimizations
- **Important:** Caching is automatically disabled when `--use-xtb-wbos` is set

## Output Formats

**fit --csv**: Adds columns to input CSV:
- `xtb_energy_eh`, `xtb_energy_kcal_mol`: xTB total energies (when `--optimizer xtb`)
- `mlip_energy_eh`, `mlip_energy_kcal_mol`: MLIP total energies (when `--optimizer uma/esen/aimnet2`; field is named `mlip_energy` regardless of specific MLIP used)
- `gxtb_energy_eh`, `gxtb_energy_kcal_mol`: gxtb single-point energies (if `--use-gxtb` is set)
- `pred_dhf_4param`, `pred_dhf_7param`, etc.: Predicted ΔHf° for each model (kcal/mol)
- `error_4param`, `error_7param`, etc.: Prediction error (pred - exp) for each model
- `error`: Error message for failed molecules (None if successful)

**predict -o**: Creates CSV with columns:
- All input columns (including `smiles`, `name` if present)
- `xtb_energy_kcal_mol` or `mlip_energy_kcal_mol`: optimizer energy
- `gxtb_energy_kcal_mol`: gxtb single-point energy (if `--use-gxtb` is set)
- `pred_dhf_kcal_mol`: Predicted ΔHf° using specified model
- `error`: Error message (None if successful)
