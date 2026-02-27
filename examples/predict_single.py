"""Predict ΔHf° for a single molecule using the Python API.

Requires xTB on PATH (or use --optimizer uma for GPU).
"""

import json
from pathlib import Path

from deltahf.pipeline import process_molecule

# Load default xTB parameters
params_file = Path(__file__).resolve().parent.parent / "params" / "params_xtb.json"
with open(params_file) as f:
    params = json.load(f)

# Process nitromethane with the bondorder_ext model (best performer)
result = process_molecule(
    "C[N+](=O)[O-]",
    n_conformers=1,
    name="nitromethane",
    epsilon_bondorder_ext=params["bondorder_ext"],
)

if result.error:
    print(f"Error: {result.error}")
else:
    print(f"Molecule:  {result.name} ({result.smiles})")
    print(f"xTB energy: {result.xtb_energy:.6f} Eh")
    print(f"xTB energy: {result.xtb_energy_kcal:.2f} kcal/mol")
    print(f"ΔHf° (bondorder_ext): {result.dhf_bondorder_ext:.2f} kcal/mol")
    print(f"Experimental ΔHf°:    -17.8 kcal/mol")
    print()
    print(f"Atom counts (element):      {result.atom_counts_element}")
    print(f"Atom counts (bondorder_ext): {result.atom_counts_bondorder_ext}")
    print(f"Conformers: {result.n_conformers_generated} generated, "
          f"{result.n_conformers_optimized} optimized")
