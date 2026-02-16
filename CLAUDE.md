# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python program replicating the workflow from "Atom Equivalent Energies for the Rapid Estimation of the Heat of Formation of Explosive Molecules from Density Functional Tight Binding Theory". The reference paper PDF is in the repo root.

## Pipeline Architecture

For a given CSV of SMILES strings, the pipeline processes each molecule through these stages (with error handling at each step):

1. **SMILES parsing** — convert SMILES to RDKit mol object
2. **3D embedding** — generate 3D coordinates with RDKit
3. **Conformer search** — ETKDG conformer analysis with RDKit
4. **Geometry optimization** — xTB optimization of the lowest *n* conformers (n is a CLI argument)
5. **Heat of Formation prediction** — apply atom equivalent formula to the optimized geometry

## Parameterization

Atom equivalent energies must first be parameterized using experimental heats of formation. Training data sources:
- Two tables of experimental data in the reference PDF
- NIST Chemistry WebBook: https://webbook.nist.gov/chemistry/

These should be extracted and stored as structured data (CSV) for the fitting step.

## Key Dependencies

- `rdkit` — molecular representation, 3D embedding, ETKDG conformer generation
- `xtb-python` or `xtb` CLI — semi-empirical quantum chemistry optimization
- `pandas` — CSV I/O and data handling
- `numpy` / `scipy` — numerical fitting of atom equivalents

## Commands

No build/test/lint infrastructure exists yet. As it is set up:
- Use `pyproject.toml` for packaging
- Use `pytest` for testing (`pytest` to run all, `pytest tests/test_file.py::test_name` for a single test)
- Use `ruff` for linting and formatting
