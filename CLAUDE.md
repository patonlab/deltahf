# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python package replicating the workflow from Cawkwell et al., "Atom Equivalent Energies for the Rapid Estimation of the Heat of Formation of Explosive Molecules" (J. Chem. Inf. Model. 2021, 61, 3337–3347), using **GFN2-xTB** instead of DFTB/lanl31.

Core equation: `ΔHf° = u_xtb − Σ(nl × εl)` where `u_xtb` is xTB total energy (kcal/mol), `nl` is atom count, and `εl` is the fitted atom equivalent energy.

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
python -m deltahf fit -i deltahf/data/training_data.csv --model both --kfold 10 --n-conformers 5
python -m deltahf predict -i molecules.csv --epsilon params.json --n-conformers 5 -o results.csv
```

## Architecture

**Two models**: 4-parameter (C, H, N, O counts) and 7-parameter (adds C', N', O' for multiply-bonded atoms, bond order threshold > 1.25).

Pipeline flow per molecule: SMILES → RDKit mol → ETKDG conformers → MMFF ranking → xTB optimization of lowest n → atom equivalent ΔHf prediction.

### Key Modules

- `deltahf/smiles.py` — SMILES→mol, `count_atoms()` (4-param), `classify_atoms_7param()` (7-param with resonance-aware bond order detection)
- `deltahf/conformers.py` — ETKDG conformer generation, MMFF energy ranking, XYZ file writing
- `deltahf/xtb.py` — xTB CLI subprocess wrapper, energy parsing (regex on stdout)
- `deltahf/atom_equivalents.py` — least-squares fitting (`np.linalg.lstsq`), k-fold cross-validation, RMSD/max deviation statistics
- `deltahf/pipeline.py` — end-to-end orchestration: `process_molecule()`, `process_csv()`
- `deltahf/__main__.py` — CLI with `fit` and `predict` subcommands
- `deltahf/data/` — 103-molecule training dataset (CSV) and `load_training_data()` helper

### Constants

`HARTREE_TO_KCAL = 627.5094740631` defined in both `xtb.py` and `atom_equivalents.py`.

### Bond Order Classification (7-param)

Uses RDKit `GetBondTypeAsDouble()` from the molecular graph. An atom is "primed" (multiply-bonded) if it has any bond > 1.25 OR has nonzero formal charge with a neighbor that has a bond > 1.25 (handles charge-separated resonance forms like `[N+](=O)[O-]`).

## Training Data

`deltahf/data/training_data.csv` contains 103 molecules from Table 1 of the paper with columns: `id, name, formula, smiles, exp_dhf_kcal_mol`. SMILES were curated from NIST WebBook (https://webbook.nist.gov/chemistry/) and validated against RDKit-computed molecular formulas.
