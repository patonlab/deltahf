"""Predicted xTB ΔHf° histogram: Cawkwell training set vs random ZINC drug-like sample."""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem

from deltahf.smiles import SUPPORTED_ELEMENTS, heavy_atom_count, total_atom_count

PLOTS_DIR = Path(__file__).parent
REPO_ROOT = PLOTS_DIR.parent

CAWKWELL_CSV = PLOTS_DIR / "cawkwell_si_atom_counts.csv"
ZINC_CSV = PLOTS_DIR / "250k_rndm_zinc_drugs_clean_3.csv"
ZINC_SAMPLE_CSV = PLOTS_DIR / "zinc_sample_1000.csv"
CAWKWELL_INPUT_CSV = PLOTS_DIR / "cawkwell_input.csv"
CAWKWELL_OUT_CSV = PLOTS_DIR / "cawkwell_gxtb_predictions.csv"
ZINC_OUT_CSV = PLOTS_DIR / "zinc_gxtb_predictions.csv"
CACHE_DIR = PLOTS_DIR / ".cache_zinc_cawkwell_gxtb"
OUTPUT_PNG = PLOTS_DIR / "zinc_vs_cawkwell_histogram.png"


def neutralize_smiles(smiles: str) -> str:
    """Remove formal charges (protonated amines, carboxylates, etc.)."""
    from rdkit.Chem.MolStandardize import rdMolStandardize
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    mol_neutral = rdMolStandardize.Uncharger().uncharge(mol)
    return Chem.MolToSmiles(mol_neutral)


def prepare_inputs(seed: int = 42):
    # Cawkwell: extract smiles + name
    cawk = pd.read_csv(CAWKWELL_CSV)[["smiles", "name"]]
    cawk.to_csv(CAWKWELL_INPUT_CSV, index=False)
    print(f"Cawkwell input: {len(cawk)} molecules -> {CAWKWELL_INPUT_CSV}")

    # ZINC: filter to supported elements, sample 1000, then neutralize
    zinc = pd.read_csv(ZINC_CSV)
    def has_supported_elements(smi):
        mol = Chem.MolFromSmiles(smi.strip())
        if mol is None:
            return False
        return all(a.GetSymbol() in SUPPORTED_ELEMENTS for a in mol.GetAtoms())
    zinc["smiles"] = zinc["smiles"].str.strip()
    supported = zinc[zinc["smiles"].apply(has_supported_elements)]
    print(f"ZINC: {len(supported)}/{len(zinc)} molecules have supported elements")
    sample = supported.sample(n=1000, random_state=seed)[["smiles"]]
    sample["smiles"] = sample["smiles"].apply(neutralize_smiles)
    print(f"ZINC sample: {len(sample)} molecules (neutralized) -> {ZINC_SAMPLE_CSV}")
    sample.to_csv(ZINC_SAMPLE_CSV, index=False)


def run_predict(input_csv: Path, output_csv: Path, label: str):
    if output_csv.exists():
        print(f"{label}: predictions already exist at {output_csv}, skipping.")
        return
    cmd = [
        sys.executable, "-m", "deltahf", "predict",
        "-i", str(input_csv),
        "--model", "bondorder_ext",
        "--use-gxtb",
        "--xtb-threads", "16",
        "--cache-dir", str(CACHE_DIR / label),
        "-o", str(output_csv),
    ]
    print(f"\nRunning {label}:\n  " + " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def plot(cawkwell_out: Path, zinc_out: Path, method: str = "gXTB"):
    cawk = pd.read_csv(cawkwell_out)
    zinc = pd.read_csv(zinc_out)

    n_cawk_fail = cawk["error"].notna().sum() if "error" in cawk.columns else 0
    n_zinc_fail = zinc["error"].notna().sum() if "error" in zinc.columns else 0
    if n_cawk_fail:
        print(f"Cawkwell: {n_cawk_fail} failed molecules excluded from plot")
    if n_zinc_fail:
        print(f"ZINC: {n_zinc_fail} failed molecules excluded from plot")

    for df in (cawk, zinc):
        df["n_heavy"] = df["smiles"].apply(heavy_atom_count)
        df["n_atoms"] = df["smiles"].apply(total_atom_count)
        df["dhf_per_heavy"] = df["dhf_bondorder_ext"] / df["n_heavy"]
        df["dhf_per_atom"] = df["dhf_bondorder_ext"] / df["n_atoms"]

    cawk_ok = cawk.dropna(subset=["dhf_bondorder_ext"])
    zinc_ok = zinc.dropna(subset=["dhf_bondorder_ext"])

    panels = [
        ("dhf_bondorder_ext", "Predicted ΔHf° (kcal/mol)", "Heat of formation"),
        ("dhf_per_heavy", "Predicted ΔHf° / heavy atom (kcal/mol)", "Per heavy atom"),
        ("dhf_per_atom", "Predicted ΔHf° / atom (kcal/mol)", "Per atom (incl. H)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for ax, (col, xlabel, title) in zip(axes, panels):
        c_vals = cawk_ok[col].dropna()
        z_vals = zinc_ok[col].dropna()
        all_vals = pd.concat([c_vals, z_vals])
        lo = all_vals.min() - (all_vals.max() - all_vals.min()) * 0.05
        hi = all_vals.max() + (all_vals.max() - all_vals.min()) * 0.05
        bins = np.linspace(lo, hi, 50)
        ax.hist(z_vals, bins=bins, color="#4C72B0", alpha=0.70,
                label=f"ZINC drug-like (n={len(z_vals)})")
        ax.hist(c_vals, bins=bins, color="#DD4949", alpha=0.85,
                label=f"Cawkwell training set (n={len(c_vals)})")
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel("Count", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].legend(framealpha=0.9, fontsize=10)

    fig.suptitle(f"Predicted ΔHf° ({method} + bondorder_ext): training set vs drug-like molecules",
                 fontsize=13, y=1.02)
    output_png = PLOTS_DIR / f"zinc_vs_cawkwell_{method.lower()}_histogram.png"
    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    print(f"\nSaved {output_png}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for ZINC sampling (default: 42)")
    parser.add_argument("--plot-only", action="store_true",
                        help="Skip predict step and go straight to plotting")
    args = parser.parse_args()

    if not args.plot_only:
        prepare_inputs(seed=args.seed)
        run_predict(CAWKWELL_INPUT_CSV, CAWKWELL_OUT_CSV, "cawkwell")
        run_predict(ZINC_SAMPLE_CSV, ZINC_OUT_CSV, "zinc")

    plot(CAWKWELL_OUT_CSV, ZINC_OUT_CSV, method="gXTB")


if __name__ == "__main__":
    main()
