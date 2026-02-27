"""Histogram of training data ΔHf° values, coloured by energetic vs non-energetic."""

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

is_energetic = df["category"] == "energetic"
n_energetic = is_energetic.sum()
n_non = (~is_energetic).sum()

panels = [
    ("exp_dhf_kcal_mol", "ΔHf° (kcal/mol)", "Heat of formation"),
    ("dhf_per_heavy", "ΔHf° / heavy atom (kcal/mol)", "Per heavy atom"),
    ("dhf_per_atom", "ΔHf° / atom (kcal/mol)", "Per atom (incl. H)"),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

for ax, (col, xlabel, title) in zip(axes, panels):
    vals = df[col]
    bins = np.linspace(vals.min() - (vals.max() - vals.min()) * 0.05,
                       vals.max() + (vals.max() - vals.min()) * 0.05, 50)
    ax.hist(df.loc[~is_energetic, col], bins=bins, color="#4C72B0", alpha=0.75,
            label=f"Non-energetic (n={n_non})")
    ax.hist(df.loc[is_energetic, col], bins=bins, color="#DD4949", alpha=0.85,
            label=f"Energetic (n={n_energetic})")
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)

axes[0].legend(framealpha=0.9, fontsize=10)

fig.suptitle("Training data: heat of formation distributions", fontsize=13, y=1.02)
plt.tight_layout()
output = Path(__file__).parent / "training_data_histogram.png"
plt.savefig(output, dpi=150, bbox_inches="tight")
print(f"Saved {output}")
