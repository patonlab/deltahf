"""Compare ΔHf° predictions across all eight atom equivalent models.

Requires xTB on PATH. Processes a single molecule and shows predictions
from every model side by side.
"""

import json
from pathlib import Path

from deltahf.pipeline import process_molecule

# Load default xTB parameters
params_file = Path(__file__).resolve().parent.parent / "params" / "params_xtb.json"
with open(params_file) as f:
    params = json.load(f)

# Pass all model epsilons so every prediction is computed
result = process_molecule(
    "c1ccccc1",
    n_conformers=1,
    name="benzene",
    epsilon_element=params["element"],
    epsilon_element_bo=params["element_bo"],
    epsilon_hybrid=params["hybrid"],
    epsilon_bondorder=params["bondorder"],
    epsilon_bondorder_ext=params["bondorder_ext"],
    epsilon_bondorder_ar=params["bondorder_ar"],
    epsilon_extended=params["extended"],
    epsilon_neighbour=params["neighbour"],
)

if result.error:
    print(f"Error: {result.error}")
    raise SystemExit(1)

exp_dhf = 19.8  # Experimental ΔHf° for benzene (kcal/mol)

print(f"Molecule: {result.name} ({result.smiles})")
print(f"xTB energy: {result.xtb_energy:.6f} Eh ({result.xtb_energy_kcal:.2f} kcal/mol)")
print(f"Experimental ΔHf°: {exp_dhf:.1f} kcal/mol")
print()

models = [
    ("element",       result.dhf_element),
    ("element_bo",    result.dhf_element_bo),
    ("hybrid",        result.dhf_hybrid),
    ("bondorder",     result.dhf_bondorder),
    ("bondorder_ext", result.dhf_bondorder_ext),
    ("bondorder_ar",  result.dhf_bondorder_ar),
    ("extended",      result.dhf_extended),
    ("neighbour",     result.dhf_neighbour),
]

print(f"{'Model':18s} {'ΔHf° (kcal/mol)':>16s} {'Error':>10s}")
print(f"{'-' * 18} {'-' * 16} {'-' * 10}")
for name, dhf in models:
    if dhf is not None:
        error = dhf - exp_dhf
        print(f"{name:18s} {dhf:16.2f} {error:+10.2f}")
    else:
        print(f"{name:18s} {'N/A':>16s}")
