# Analysis

Scripts and data for visualising the deltahf training set and comparing predicted heats of formation across chemical domains.

---

## Training Data Distribution

**Script:** `plot_training_data.py`

```bash
python plot_training_data.py
```

The training data (`deltahf/data/training_data.csv`) contains 313 molecules split into energetic (n=45) and non-energetic (n=268) categories. The histograms below show absolute ΔHf°, ΔHf° per heavy atom, and ΔHf° per atom (including H). The energetic molecules are notably shifted toward positive ΔHf° values, and the per-atom normalisation reveals a clearer separation between the two groups.

![Training data ΔHf° distributions](training_data_histogram.png)

---

## ZINC vs Cawkwell Comparison

**Script:** `plot_zinc_vs_cawkwell.py`

This script compares predicted ΔHf° distributions for two molecule sets:

- **Cawkwell training set** (531 molecules) — the full deltahf training data
- **ZINC drug-like sample** (1,000 molecules) — randomly sampled from the ZINC 250k drug-like dataset, filtered to supported elements and neutralised

The comparison assesses whether the training set covers the chemical space of typical drug-like molecules and how the predicted ΔHf° distributions differ.

### Usage

```bash
# Full pipeline: prepare inputs, run predictions, plot
python plot_zinc_vs_cawkwell.py

# Plot only (if predictions already exist)
python plot_zinc_vs_cawkwell.py --plot-only
```

By default the script uses gXTB + `bondorder_ext`. The `comparison_workflow` file contains commands for running the xTB variant separately.

### gXTB predictions

Using gXTB + `bondorder_ext`, the ZINC drug-like molecules tend toward more negative predicted ΔHf° than the training set. When normalised per heavy atom, the distributions overlap more substantially, suggesting the atom equivalent model can extrapolate reasonably to drug-like chemical space.

![gXTB: ZINC vs Cawkwell ΔHf° distributions](zinc_vs_cawkwell_gxtb_histogram.png)

### xTB predictions

The same comparison using xTB + `bondorder_ext` shows a similar pattern. The distributions are broader due to the lower accuracy of xTB energies, but the relative shift between training set and drug-like molecules is consistent.

![xTB: ZINC vs Cawkwell ΔHf° distributions](zinc_vs_cawkwell_xtb_histogram.png)

---

## Data Files

| File | Description |
|------|-------------|
| `250k_rndm_zinc_drugs_clean_3.csv` | ZINC 250k drug-like dataset (source data) |
| `zinc_sample_1000.csv` | 1,000-molecule random sample (neutralised, supported elements only) |
| `cawkwell_input.csv` | Cawkwell training set formatted for prediction |
| `cawkwell_gxtb_predictions.csv` | gXTB + bondorder_ext predictions for Cawkwell set |
| `zinc_gxtb_predictions.csv` | gXTB + bondorder_ext predictions for ZINC sample |
| `comparison_workflow` | Shell commands for the xTB comparison workflow |
