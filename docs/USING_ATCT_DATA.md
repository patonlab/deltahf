# Using ATcT Data for Training

## Overview

The [Active Thermochemical Tables (ATcT)](https://atct.anl.gov/) database contains
high-quality experimental heats of formation. This guide explains how to extract ATcT
data and add it to the deltahf training set.

## Workflow

### Step 1: Manual export from ATcT

ATcT does not offer bulk download. You must copy data manually from the website.

1. Visit <https://atct.anl.gov/Thermochemical%20Data/version%201.220/index.php>
2. Select the species of interest
3. Copy the relevant rows into a CSV with the following columns:

```
name,formula,dhf_298K_kJ_mol,uncertainty_kJ_mol
methane,CH4,-74.87,0.32
hydrogen fluoride,HF,-272.68,0.019
sulfur hexafluoride,F6S,-1220.5,0.4
```

- Use the **gas-phase** ΔHf° at **298.15 K** (the values shown in red on the website).
- Values are in **kJ/mol** — the script converts to kcal/mol automatically.

### Step 2: Save the ATcT HTML page (recommended)

The ATcT website page contains SMILES and InChI for every species in its JavaScript
autocomplete. Saving this page gives near-complete SMILES coverage without any API
calls. In your browser:

1. Visit <https://atct.anl.gov/Thermochemical%20Data/version%201.220/index.php>
2. Save the page as HTML (`File → Save Page As… → Webpage, Complete` or just
   `Webpage, HTML Only` — the JS data is embedded in the page source).
3. Save to `docs/` alongside the CSV.

### Step 3: Run the extraction script

```bash
python scripts/extract_atct.py \
    -i docs/ATcT.csv \
    -o docs/atct_out.csv \
    --html "docs/ATcT Thermochemical Values ver. 1.220.html" \
    --existing deltahf/data/training_data.csv
```

The script will:
- Read SMILES directly from the HTML page for any species with a matching CAS number
  (covers ~98% of entries; no API calls for those).
- Fall back to the PubChem REST API for the remaining ~2%.
- Cross-check the RDKit formula against the ATcT formula to catch wrong PubChem matches.
- Reject radicals, ions, condensed-phase entries, and unsupported elements.
- Deduplicate against the existing `training_data.csv` using canonical SMILES.
- Write an output CSV with the same columns as `training_data.csv`.

Use `--skip-pubchem` to run entirely offline (HTML-only, no API calls).
Use `-v` / `--verbose` to see per-molecule decisions.

### Step 3: Review and append

Inspect `atct_out.csv` before appending. When satisfied:

```bash
# Append (skip header line)
tail -n +2 atct_out.csv >> deltahf/data/training_data.csv
```

Update `tests/test_data.py` to reflect the new row count and any new source/category values.

### Step 4: Refit

```bash
python -m deltahf fit \
    -i deltahf/data/training_data.csv \
    --model all \
    --kfold 10 \
    --n-conformers 1 \
    --cache-dir .cache_xtb \
    -o params.json
```

## Pipeline statistics (ATcT v1.220, full database, run 2026-02-19)

The table below shows how the 3442-row ATcT v1.220 export was reduced to the
233 molecules appended to `deltahf/data/training_data.csv`.

| Stage | Count | Reason |
|---|---:|---|
| Raw ATcT rows | 3442 | |
| − Condensed phase | −729 | `(cr`, `(l)`, `(aq)` in FORMULA field |
| − Open-shell gas | −255 | `triplet`, `doublet`, `radical` in FORMULA field |
| After phase/spin pre-filter | 2458 | |
| CAS deduplication | →1740 | Multiple ATcT entries per species (conformers, symmetry variants); kept lowest DHF gas-phase entry per CAS |
| − No SMILES found | −29 | Not present in ATcT HTML autocomplete |
| − Invalid SMILES | −167 | ATcT SMILES not parseable by RDKit |
| − Formula mismatch | −345 | RDKit formula from SMILES ≠ ATcT stated formula |
| − Charged / radical | −819 | Non-zero formal charge or radical electrons |
| − Unsupported elements | −101 | Elements outside {C, H, N, O, F, S, Cl} (Br, P, Si, …) |
| − Duplicates | −47 | Canonical SMILES already in `training_data.csv` |
| **Accepted** | **233** | |
| − Name contains "anion"/"cation" | −3 | Net-neutral zwitterions or mis-named species removed manually |
| **Final** | **230** | |

The HTML autocomplete covered 1711/1740 candidates (98%), eliminating almost all
PubChem API calls. The formula mismatch filter (345 rejections) was essential:
without it, silently wrong structures from the HTML SMILES would enter the training
set. Three net-neutral species with "anion" in their ATcT name (Nitrosyl nitrite
anion N2O3, Oxygen difluoride anion F2O, Isocyanogen anion C2N2) slipped through
the charge filter because RDKit reports them as uncharged; these were removed
manually.

## Notes

- **Neutral, closed-shell molecules only**: ATcT contains ions and radicals; the script
  filters these automatically (`Chem.GetFormalCharge(mol) == 0` and no radical electrons).
- **Category labels**: assigned automatically based on elements present (`fluorinated`,
  `sulfur`, `chlorinated`, `hydrocarbon`, `small_CHNO`).
- **Uncertainty column**: stored in the input CSV but not used by deltahf — the script
  ignores it. You may wish to filter out high-uncertainty entries manually.
