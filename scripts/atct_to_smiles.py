"""
Convert ATcT data to deltahf training format.

Takes raw ATcT data (name, formula, ΔHf°) and converts to SMILES.
Filters for CHNO-only molecules compatible with deltahf.

Usage:
    python scripts/atct_to_smiles.py -i atct_raw.csv -o atct_training.csv
"""

import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Common molecules with known SMILES
# You'll need to expand this or use automated conversion
KNOWN_SMILES = {
    "methane": "C",
    "ethane": "CC",
    "propane": "CCC",
    "butane": "CCCC",
    "methanol": "CO",
    "ethanol": "CCO",
    "water": "O",
    "hydrogen": "[H][H]",
    "oxygen": "O=O",
    "nitrogen": "N#N",
    "carbon dioxide": "O=C=O",
    "carbon monoxide": "[C-]#[O+]",
    "ammonia": "N",
    "nitromethane": "C[N+](=O)[O-]",
    "nitroethane": "CC[N+](=O)[O-]",
    "benzene": "c1ccccc1",
    "toluene": "Cc1ccccc1",
    "acetylene": "C#C",
    "ethylene": "C=C",
    "formaldehyde": "C=O",
    "acetaldehyde": "CC=O",
    "acetone": "CC(=O)C",
    "acetic acid": "CC(=O)O",
    "formic acid": "C(=O)O",
    # Add more as needed...
}

def name_to_smiles(name):
    """
    Convert chemical name to SMILES.

    This is a simplified version - for production use:
    1. Use PubChem API to look up SMILES by name
    2. Use chemical name parsing libraries
    3. Manually curate a comprehensive mapping
    """
    name_lower = name.lower().strip()

    if name_lower in KNOWN_SMILES:
        return KNOWN_SMILES[name_lower]

    # Try using RDKit's name-to-structure (limited)
    try:
        mol = Chem.MolFromSmiles(name)  # Won't work for names
        if mol:
            return Chem.MolToSmiles(mol)
    except:
        pass

    return None

def is_chno_only(smiles):
    """Check if molecule contains only C, H, N, O."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False

    allowed_atoms = {'C', 'H', 'N', 'O'}
    atoms = {atom.GetSymbol() for atom in mol.GetAtoms()}

    return atoms.issubset(allowed_atoms)

def convert_atct_data(input_csv, output_csv):
    """Convert ATcT data to deltahf training format."""
    df = pd.read_csv(input_csv)

    # Add SMILES column
    df['smiles'] = df['name'].apply(name_to_smiles)

    # Filter for molecules with SMILES
    df = df[df['smiles'].notna()].copy()

    # Filter for CHNO-only
    df['is_chno'] = df['smiles'].apply(is_chno_only)
    df = df[df['is_chno']].copy()

    # Rename columns to match deltahf format
    df = df.rename(columns={
        'dhf_298K_kcal_mol': 'exp_dhf_kcal_mol'
    })

    # Keep only required columns
    output_cols = ['smiles', 'name', 'formula', 'exp_dhf_kcal_mol']
    df = df[output_cols].copy()

    # Add source and category
    df['source'] = 'ATcT_v1.220'
    df['category'] = 'atct'

    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} CHNO molecules to {output_csv}")
    print(f"\nSample:")
    print(df.head())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert ATcT data to deltahf format")
    parser.add_argument("-i", "--input", required=True, help="Input CSV from ATcT scraping")
    parser.add_argument("-o", "--output", default="atct_training.csv", help="Output training CSV")
    args = parser.parse_args()

    convert_atct_data(args.input, args.output)
