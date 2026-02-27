"""Predict ΔHf° for a batch of molecules from a CSV file.

Requires xTB on PATH (or use optimizer="uma" for GPU).
This is the Python API equivalent of the CLI command:

    python -m deltahf predict -i molecules.csv --model bondorder_ext -o predictions.csv
"""

import json
from pathlib import Path

from deltahf.pipeline import process_csv

# Load default xTB parameters
params_file = Path(__file__).resolve().parent.parent / "params" / "params_xtb.json"
with open(params_file) as f:
    params = json.load(f)

# Process all molecules in the CSV
input_csv = Path(__file__).resolve().parent / "molecules.csv"
output_csv = Path(__file__).resolve().parent / "predictions.csv"

df = process_csv(
    input_csv,
    n_conformers=1,
    epsilon_element=params["element"],
    epsilon_bondorder_ext=params["bondorder_ext"],
    output_path=output_csv,
)

# Display results
print(f"Processed {len(df)} molecules\n")
for _, row in df.iterrows():
    if row.get("error"):
        print(f"  {row['name']:20s}  ERROR: {row['error']}")
    else:
        print(f"  {row['name']:20s}  "
              f"element={row['dhf_element']:8.2f}  "
              f"bondorder_ext={row['dhf_bondorder_ext']:8.2f} kcal/mol")

print(f"\nResults saved to {output_csv}")
