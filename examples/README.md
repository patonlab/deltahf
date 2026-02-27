# Examples

This directory contains example scripts demonstrating how to use `deltahf` via the CLI and Python API.

---

## Prerequisites

Install deltahf and its dependencies (see the main [README](../README.md)):

```bash
conda env create -f environment.yml
conda activate deltahf
pip install -e .
```

For CLI and pipeline examples, **xTB must be on your PATH**. The atom classification example (`atom_classification.py`) requires only RDKit.

---

## Example Data

[`molecules.csv`](molecules.csv) — a small set of 10 molecules used by the batch examples:

```csv
smiles,name
CC(=O)O,acetic acid
c1ccccc1,benzene
CCO,ethanol
C(=O)=O,carbon dioxide
C[N+](=O)[O-],nitromethane
CC(C)C,isobutane
C1CCCCC1,cyclohexane
c1ccncc1,pyridine
CC=O,acetaldehyde
O,water
```

Any CSV with a `smiles` column (and optionally a `name` column) can be used as input.

---

## CLI Usage

### Predict ΔHf° for new molecules

The simplest way to use deltahf is through the CLI. Default parameters (fitted on the 531-molecule training set) are used automatically:

```bash
# Predict with the element model (simplest, 7 parameters)
python -m deltahf predict -i molecules.csv --model element -o predictions.csv

# Predict with the bondorder_ext model (best accuracy)
python -m deltahf predict -i molecules.csv --model bondorder_ext -o predictions.csv
```

### Choose a geometry optimizer

Three optimizer backends are available, each requiring its own matched parameter set:

```bash
# xTB (default, CPU) — uses params/params_xtb.json automatically
python -m deltahf predict -i molecules.csv --model bondorder_ext -o predictions.csv

# gxtb (CPU, higher accuracy) — uses params/params_gxtb.json automatically
python -m deltahf predict -i molecules.csv --model bondorder_ext --use-gxtb -o predictions.csv

# UMA (GPU, best accuracy) — uses params/params_uma.json automatically
python -m deltahf predict -i molecules.csv --model bondorder_ext --optimizer uma -o predictions.csv
```

> **Warning:** Each optimizer produces energies on a different scale. Parameters fitted with one optimizer must **never** be used with another. deltahf validates this automatically when using the default parameter files.

### Use custom parameters

If you have fitted your own parameters (see below), pass them with `--epsilon`:

```bash
python -m deltahf predict -i molecules.csv --epsilon my_params.json --model bondorder_ext -o predictions.csv
```

### Fit new parameters on custom training data

To fit atom equivalent energies on your own training set, provide a CSV with `smiles` and `exp_dhf_kcal_mol` columns:

```bash
# Fit all 8 models with 10-fold cross-validation
python -m deltahf fit \
    -i my_training_data.csv \
    --model all \
    --kfold 10 \
    --n-conformers 1 \
    -o my_params.json

# Fit a single model
python -m deltahf fit \
    -i my_training_data.csv \
    --model bondorder_ext \
    --kfold 10 \
    --n-conformers 1 \
    -o my_params.json

# Show the 10 worst outliers after fitting
python -m deltahf fit \
    -i my_training_data.csv \
    --model bondorder_ext \
    --kfold 10 \
    --outliers 10 \
    -o my_params.json
```

### Increase conformer sampling

By default, only the single lowest-energy RDKit conformer is optimized. For more thorough sampling (at the cost of computation time):

```bash
python -m deltahf predict -i molecules.csv --model bondorder_ext --n-conformers 5 -o predictions.csv
```

> **Note:** Benchmarking shows <1% RMSD improvement going from 1 to 5 conformers, so `--n-conformers 1` (the default) is recommended for most use cases.

---

## Python API

### Predict ΔHf° for a single molecule

[`predict_single.py`](predict_single.py) — process one molecule and print results:

```python
from deltahf.pipeline import process_molecule
import json
from pathlib import Path

params_file = Path("params/params_xtb.json")
with open(params_file) as f:
    params = json.load(f)

result = process_molecule(
    "C[N+](=O)[O-]",
    n_conformers=1,
    name="nitromethane",
    epsilon_bondorder_ext=params["bondorder_ext"],
)

print(f"ΔHf° = {result.dhf_bondorder_ext:.2f} kcal/mol")
print(f"xTB energy = {result.xtb_energy:.6f} Eh")
```

Run it:
```bash
python examples/predict_single.py
```

### Predict ΔHf° for a batch of molecules

[`predict_batch.py`](predict_batch.py) — process a CSV and write results:

```python
from deltahf.pipeline import process_csv
from pathlib import Path
import json

with open("params/params_xtb.json") as f:
    params = json.load(f)

df = process_csv(
    Path("examples/molecules.csv"),
    n_conformers=1,
    epsilon_element=params["element"],
    epsilon_bondorder_ext=params["bondorder_ext"],
    output_path=Path("predictions.csv"),
)
```

Run it:
```bash
python examples/predict_batch.py
```

### Explore atom classification schemes

[`atom_classification.py`](atom_classification.py) — see how each model classifies atoms (no xTB required, RDKit only):

```python
from deltahf.smiles import count_atoms, classify_atoms_hybrid, classify_atoms_bondorder_ext

print(count_atoms("c1ccccc1"))
# {'C': 6, 'H': 6}

print(classify_atoms_hybrid("c1ccccc1"))
# {'C_sp2': 6, 'H': 6, ...}

print(classify_atoms_bondorder_ext("c1ccccc1"))
# {'C_2_1H': 6, 'H': 6, ...}
```

Run it:
```bash
python examples/atom_classification.py
```

### Compare all eight models

[`compare_models.py`](compare_models.py) — predict ΔHf° with every model and compare:

```bash
python examples/compare_models.py
```

Output:
```
Molecule: benzene (c1ccccc1)
xTB energy: -8.370133 Eh (...)
Experimental ΔHf°: 19.8 kcal/mol

Model              ΔHf° (kcal/mol)      Error
------------------ ---------------- ----------
element                        ...       +...
element_bo                     ...       +...
hybrid                         ...       +...
bondorder                      ...       +...
bondorder_ext                  ...       +...
bondorder_ar                   ...       +...
extended                       ...       +...
neighbour                      ...       +...
```

---

## Available Atom Classification Functions

All functions take a SMILES string and return a `dict[str, int]` of atom type counts:

| Function | Model | Description |
|----------|-------|-------------|
| `count_atoms()` | `element` | One parameter per element (C, H, N, O, F, S, Cl) |
| `classify_atoms_7param()` | `element_bo` | Adds primed variants for multiply-bonded atoms (C', N', O', S') |
| `classify_atoms_hybrid()` | `hybrid` | RDKit hybridization: sp3, sp2, sp per element |
| `classify_atoms_bondorder()` | `bondorder` | Max Kekulized bond order: `_1`, `_2`, `_3` suffix |
| `classify_atoms_bondorder_ext()` | `bondorder_ext` | Bond order + H-count for carbon (best performer) |
| `classify_atoms_bondorder_ar()` | `bondorder_ar` | Bond order without Kekulization (preserves `_ar`) |
| `classify_atoms_extended()` | `extended` | Hybridization + H-count for carbon |
| `classify_atoms_neighbour()` | `neighbour` | Extended + heavy-atom neighbour type for N/O |

---

## Parameter Files

Default parameter files are in [`params/`](../params/):

| File | Optimizer | Best CV RMSD (kcal/mol) |
|------|-----------|:-----------------------:|
| `params_xtb.json` | xTB (CPU) | 7.09 |
| `params_gxtb.json` | gxtb (CPU) | 4.11 |
| `params_uma.json` | UMA (GPU) | 2.78 |

Each JSON file contains fitted atom equivalent energies for all eight models, plus metadata. Extract a specific model's parameters with:

```python
import json

with open("params/params_xtb.json") as f:
    params = json.load(f)

epsilon = params["bondorder_ext"]  # dict mapping atom types to energies
```
