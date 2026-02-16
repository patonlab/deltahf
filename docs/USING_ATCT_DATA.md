# Using ATcT Data for Training

## Overview

The [Active Thermochemical Tables (ATcT)](https://atct.anl.gov/) database contains high-quality experimental heats of formation for **3,444 species**. This guide explains how to use ATcT data to train deltahf models.

## Why Use ATcT Data?

**Advantages:**
- ✅ High-quality experimental values with uncertainties
- ✅ Larger dataset (3,444 species vs current 314)
- ✅ Covers wide range of molecules
- ✅ Regularly updated (version 1.220 as of March 2025)

**Challenges:**
- ❌ No direct CSV export (requires web scraping)
- ❌ Uses chemical names, not SMILES (requires conversion)
- ❌ Includes many non-CHNO molecules (need filtering)
- ❌ Some species may be radicals/ions (not suitable for neutral molecule predictions)

## Workflow

### Step 1: Extract ATcT Data

The ATcT website doesn't provide bulk download, so you need to scrape it:

```bash
# Option A: Manual copy-paste
# 1. Visit https://atct.anl.gov/Thermochemical%20Data/version%201.220/index.php
# 2. Use browser developer tools to extract table data
# 3. Save as CSV with columns: name, formula, dhf_298K_kcal_mol, uncertainty

# Option B: Automated scraping (requires Selenium)
python scripts/scrape_atct.py -o atct_raw.csv
```

**Expected format:**
```csv
name,formula,dhf_298K_kcal_mol,uncertainty
methane,CH4,-17.89,0.08
ethane,C2H6,-20.04,0.12
benzene,C6H6,19.82,0.12
```

### Step 2: Convert Names to SMILES

This is the **hardest step** because ATcT uses chemical names, not SMILES.

**Recommended approaches:**

#### A. Use PubChem API (Best for automation)

```python
import requests

def pubchem_name_to_smiles(name):
    """Look up SMILES from PubChem by chemical name."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url)
        data = response.json()
        return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except:
        return None

# Example
smiles = pubchem_name_to_smiles("benzene")
print(smiles)  # c1ccccc1
```

#### B. Manual Curation (Most reliable)

Create a mapping file for important molecules, especially those not in PubChem.

#### C. Use the Provided Script (Requires expansion)

```bash
python scripts/atct_to_smiles.py -i atct_raw.csv -o atct_training.csv
```

You'll need to expand the `KNOWN_SMILES` dictionary in the script for comprehensive coverage.

### Step 3: Filter for CHNO Molecules

deltahf currently supports only C, H, N, O. The conversion script automatically filters:

```python
def is_chno_only(smiles):
    mol = Chem.MolFromSmiles(smiles)
    allowed_atoms = {'C', 'H', 'N', 'O'}
    atoms = {atom.GetSymbol() for atom in mol.GetAtoms()}
    return atoms.issubset(allowed_atoms)
```

### Step 4: Validate and Merge

```bash
# Check the converted data
head atct_training.csv

# Validate SMILES
python scripts/validate_training_data.py atct_training.csv

# Merge with existing training data
python scripts/merge_training_data.py \
    deltahf/data/training_data.csv \
    atct_training.csv \
    -o deltahf/data/training_data_expanded.csv
```

### Step 5: Train Models

```bash
# Train with expanded dataset
python -m deltahf fit \
    -i deltahf/data/training_data_expanded.csv \
    --model all \
    --kfold 10 \
    --n-conformers 5 \
    --cache-dir .cache/atct_training \
    -o params_atct.json
```

## Important Considerations

### 1. **Neutral Molecules Only**

ATcT includes ions and radicals. For deltahf training:
- ✅ Use: neutral closed-shell molecules
- ❌ Exclude: radicals (e.g., •OH, •CH3)
- ❌ Exclude: ions (e.g., H+, NO2-)

Filter by checking formal charge:
```python
mol = Chem.MolFromSmiles(smiles)
total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
is_neutral = (total_charge == 0)
```

### 2. **Gas-Phase Values**

Ensure ΔHf° values are **gas-phase** at **298.15 K**:
- ATcT provides both 0 K and 298.15 K values
- Use the **298.15 K** values (shown in red on the website)

### 3. **Unit Conversion**

ATcT might provide values in different units:
- If in **kJ/mol**: divide by 4.184 to get kcal/mol
- deltahf expects **kcal/mol**

### 4. **Overlapping Data**

Some molecules may appear in both your current training set and ATcT:
- Check for duplicates by SMILES (canonical form)
- Keep ATcT values if they're more recent or have lower uncertainty
- Or keep both for cross-validation

### 5. **Quality Control**

After conversion:
```python
# Check for issues
import pandas as pd
from rdkit import Chem

df = pd.read_csv('atct_training.csv')

# Validate all SMILES
df['valid_smiles'] = df['smiles'].apply(
    lambda s: Chem.MolFromSmiles(s) is not None
)

print(f"Invalid SMILES: {(~df['valid_smiles']).sum()}")

# Check for reasonable ΔHf° range
print(f"ΔHf° range: {df['exp_dhf_kcal_mol'].min():.1f} to {df['exp_dhf_kcal_mol'].max():.1f} kcal/mol")

# Outliers (very high/low values)
outliers = df[abs(df['exp_dhf_kcal_mol']) > 200]
print(f"Potential outliers: {len(outliers)}")
```

## Example: Complete Workflow

```bash
# 1. Scrape ATcT (manual or automated)
# ... results in atct_raw.csv

# 2. Convert to SMILES and filter
python scripts/atct_to_smiles.py \
    -i atct_raw.csv \
    -o atct_chno.csv

# 3. Validate
python -c "
import pandas as pd
from rdkit import Chem

df = pd.read_csv('atct_chno.csv')
print(f'Total molecules: {len(df)}')
print(f'CHNO only: {len(df)}')

# Check validity
valid = df['smiles'].apply(lambda s: Chem.MolFromSmiles(s) is not None)
print(f'Valid SMILES: {valid.sum()}/{len(df)}')
"

# 4. Merge with existing data
cat deltahf/data/training_data.csv atct_chno.csv > combined_training.csv

# 5. Train
python -m deltahf fit \
    -i combined_training.csv \
    --model all \
    --kfold 10 \
    -o params_atct.json
```

## Expected Results

With a larger, high-quality dataset from ATcT, you might see:
- Improved parameter estimates (lower standard errors)
- Better generalization to new molecules
- More robust cross-validation results
- Potential discovery of systematic errors in atom equivalent approach

## Alternative: PubChem CID Mapping

If ATcT provides PubChem CIDs, conversion is much easier:

```python
def cid_to_smiles(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
```

## Resources

- **ATcT Database**: https://atct.anl.gov/
- **PubChem API**: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- **RDKit Name Parsing**: Limited, best used with PubChem for lookup
- **ChemSpider API**: Alternative to PubChem for name→SMILES conversion

## Next Steps

1. **Expand KNOWN_SMILES dictionary** in `scripts/atct_to_smiles.py`
2. **Implement PubChem API integration** for automatic name→SMILES
3. **Add validation checks** for radicals, ions, and unusual species
4. **Compare ATcT vs current training data** for overlapping molecules
5. **Benchmark models** with ATcT data vs original Cawkwell/Yalamanchi data
