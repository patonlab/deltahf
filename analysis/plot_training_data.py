"""Histogram of training data ΔHf° values, coloured by category."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from deltahf.smiles import heavy_atom_count, total_atom_count

REPO_ROOT = Path(__file__).parent.parent

df = pd.read_csv(REPO_ROOT / "deltahf/data/training_data.csv")

df["n_heavy"] = df["smiles"].apply(heavy_atom_count)
df["n_atoms"] = df["smiles"].apply(total_atom_count)
df["dhf_per_heavy"] = df["exp_dhf_kcal_mol"] / df["n_heavy"]
df["dhf_per_atom"] = df["exp_dhf_kcal_mol"] / df["n_atoms"]

# Category display order and colours
CATEGORY_STYLE = {
    "cyclic_HC":      ("Cyclic HC",       "#4C72B0"),
    "small_CHNO":     ("Small CHNO",      "#55A868"),
    "energetic":      ("Energetic",       "#DD4949"),
    "chlorinated":    ("Chlorinated",     "#8172B2"),
    "hydrocarbon":    ("Hydrocarbon",     "#CCB974"),
    "fluorinated":    ("Fluorinated",     "#64B5CD"),
    "sulfur":         ("Sulfur",          "#C44E52"),
    "strained_3ring": ("Strained 3-ring", "#DA8BC3"),
    "large_HC":       ("Large HC",        "#937860"),
}

panels = [
    ("exp_dhf_kcal_mol", "ΔHf° (kcal/mol)", "Heat of formation"),
    ("dhf_per_heavy", "ΔHf° / heavy atom (kcal/mol)", "Per heavy atom"),
    ("dhf_per_atom", "ΔHf° / atom (kcal/mol)", "Per atom (incl. H)"),
]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

for ax, (col, xlabel, title) in zip(axes, panels):
    vals = df[col]
    bins = np.linspace(vals.min() - (vals.max() - vals.min()) * 0.05,
                       vals.max() + (vals.max() - vals.min()) * 0.05, 50)

    # Stack histograms by category
    for cat, (label, color) in CATEGORY_STYLE.items():
        subset = df.loc[df["category"] == cat, col]
        if len(subset) == 0:
            continue
        ax.hist(subset, bins=bins, color=color, alpha=0.75,
                label=f"{label} (n={len(subset)})", stacked=False)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)

axes[0].legend(framealpha=0.9, fontsize=8, loc="upper left")

fig.suptitle(f"Training data: experimental ΔHf° distributions ({len(df)} molecules)",
             fontsize=13, y=1.02)
plt.tight_layout()
output = Path(__file__).parent / "training_data_histogram.png"
plt.savefig(output, dpi=150, bbox_inches="tight")
print(f"Saved {output}")
