"""Demonstrate how different atom classification schemes work.

This script requires only RDKit (no xTB needed) and shows how each
model classifies atoms in a set of example molecules.
"""

from deltahf.smiles import (
    count_atoms,
    classify_atoms_7param,
    classify_atoms_hybrid,
    classify_atoms_extended,
    classify_atoms_bondorder,
    classify_atoms_bondorder_ext,
    classify_atoms_bondorder_ar,
    classify_atoms_neighbour,
)

molecules = {
    "ethane":       "CC",
    "ethylene":     "C=C",
    "acetylene":    "C#C",
    "benzene":      "c1ccccc1",
    "nitromethane": "C[N+](=O)[O-]",
    "acetic acid":  "CC(=O)O",
    "pyridine":     "c1ccncc1",
}

classifiers = [
    ("element",        count_atoms),
    ("element_bo",     classify_atoms_7param),
    ("hybrid",         classify_atoms_hybrid),
    ("bondorder",      classify_atoms_bondorder),
    ("bondorder_ext",  classify_atoms_bondorder_ext),
    ("bondorder_ar",   classify_atoms_bondorder_ar),
    ("extended",       classify_atoms_extended),
    ("neighbour",      classify_atoms_neighbour),
]

for mol_name, smiles in molecules.items():
    print(f"{'=' * 60}")
    print(f"{mol_name} ({smiles})")
    print(f"{'=' * 60}")
    for model_name, classify_fn in classifiers:
        counts = classify_fn(smiles)
        # Show only non-zero counts
        nonzero = {k: v for k, v in counts.items() if v > 0}
        print(f"  {model_name:18s} {nonzero}")
    print()
