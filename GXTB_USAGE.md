# Using gxtb with deltahf

## Overview

deltahf now supports using **gxtb** (wB97M-V/def2-TZVPPD) energies in addition to xTB (GFN2-xTB). However, **xTB and gxtb energies are on completely different scales** and **must never be mixed**.

## Energy Scale Difference

| Method | CH₄ Energy | Description |
|--------|-----------|-------------|
| **xTB** | ~-4.18 Eh = -2,620 kcal/mol | Minimal basis, valence electrons only |
| **gxtb** | ~-40.49 Eh = -25,410 kcal/mol | Large basis (def2-TZVPPD), all electrons |

**Ratio: ~10-17x difference!**

This means that atom equivalent energies (εₗ) fitted using xTB cannot be applied to gxtb energies, and vice versa.

## Proper Workflow

### Option 1: xTB Workflow (Default, Recommended)

```bash
# Fit parameters using xTB energies
python -m deltahf fit \
  -i deltahf/data/training_data.csv \
  --model all \
  --n-conformers 5 \
  --cache-dir .cache/fit \
  -o params_xtb.json

# Predict using xTB energies
python -m deltahf predict \
  -i molecules.csv \
  --epsilon params_xtb.json \
  --model 4param \
  --n-conformers 5 \
  --cache-dir .cache/predict \
  -o results.csv
```

**Note:** Cache directory will automatically become `.cache/fit_xtb/` and `.cache/predict_xtb/`

### Option 2: gxtb Workflow (Slower, Higher-Level Theory)

```bash
# Fit parameters using gxtb energies
python -m deltahf fit \
  -i deltahf/data/training_data.csv \
  --model all \
  --n-conformers 5 \
  --use-gxtb \
  --cache-dir .cache/fit \
  -o params_gxtb.json

# Predict using gxtb energies (MUST use --use-gxtb!)
python -m deltahf predict \
  -i molecules.csv \
  --epsilon params_gxtb.json \
  --model 4param \
  --n-conformers 5 \
  --use-gxtb \
  --cache-dir .cache/predict \
  -o results.csv
```

**Note:** Cache directory will automatically become `.cache/fit_gxtb/` and `.cache/predict_gxtb/`

## Automatic Safeguards

### 1. Separate Cache Directories

Cache directories are automatically suffixed with `_xtb` or `_gxtb` to prevent mixing:

```
.cache/myrun      →  .cache/myrun_xtb/      (without --use-gxtb)
.cache/myrun      →  .cache/myrun_gxtb/     (with --use-gxtb)
```

### 2. Method Metadata in Parameters

Fitted parameter files include metadata tracking which method was used:

```json
{
  "4param": {
    "C": -2399.123,
    "H": -313.456,
    ...
  },
  "_metadata": {
    "method": "gxtb",
    "n_conformers": 5,
    "n_molecules": 311
  }
}
```

### 3. Method Mismatch Detection

If you try to predict using the wrong method, you'll get a clear error:

```
WARNING: Method mismatch detected!
Parameters were fitted using: gxtb
You are predicting using: xtb
This will produce incorrect results - energies are on different scales!

Either:
  - Use --use-gxtb to match the fitting method
  - Refit parameters using xtb
```

## Benchmark with gxtb

```bash
# Benchmark gxtb performance
python benchmark.py --use-gxtb
```

This will create separate cache directories:
- `.benchmark_cache/n1_gxtb/`
- `.benchmark_cache/n3_gxtb/`
- `.benchmark_cache/n5_gxtb/`

## Performance Considerations

| Aspect | xTB | gxtb |
|--------|-----|------|
| Speed | Fast (~15s for 314 molecules) | Slow (~10x slower) |
| Accuracy | Good (R²~0.99 with n=5) | Potentially better |
| Theory Level | Semi-empirical (GFN2-xTB) | DFT (wB97M-V/def2-TZVPPD) |
| Recommended | ✓ Yes | Only if needed |

**Recommendation:** Stick with xTB unless you have a specific need for higher-level energies.

## Installation

gxtb is optional but required if using `--use-gxtb`.

**Note:** gxtb cannot be installed via conda or pip - it must be installed manually from source.

1. Download gxtb from the official repository
2. Follow the installation instructions in the gxtb documentation
3. Ensure the `gxtb` binary is on your PATH
4. Verify installation: `which gxtb`

## Technical Details

### Why Are the Energies So Different?

1. **Basis set size:**
   - xTB: Minimal basis with Slater-type orbitals
   - gxtb: Triple-zeta + polarization + diffuse (def2-TZVPPD)

2. **Core electrons:**
   - xTB: Valence-only (frozen core approximation)
   - gxtb: All electrons included

3. **Energy reference:**
   - Different zero points due to different treatments

### Why Can't We Just Offset/Scale?

The energy differences are **not constant** - they depend on:
- Number of atoms
- Element types
- Electronic structure

So you can't simply apply a correction factor. You must use consistent energies throughout the entire workflow.

## Common Errors

### Error: "gxtb binary not found"

**Cause:** gxtb is not installed or not on PATH

**Solution:** Install gxtb manually from source (see Installation section above). Verify with `which gxtb`

### Error: "Method mismatch detected"

**Cause:** Trying to use xTB-fitted parameters with gxtb energies (or vice versa)

**Solution:** Use the same method for both fit and predict, or refit parameters

### Results look wrong (R² negative, huge RMSD)

**Cause:** Accidentally mixed xTB and gxtb energies/parameters

**Solution:** Delete caches and parameter files, start fresh with consistent method

## Migration from Old Caches

If you have old benchmark caches without method suffixes:

```bash
# Backup old caches
mv .benchmark_cache .benchmark_cache_backup

# These were xTB runs (default)
mkdir -p .benchmark_cache
mv .benchmark_cache_backup/n1 .benchmark_cache/n1_xtb
mv .benchmark_cache_backup/n3 .benchmark_cache/n3_xtb
mv .benchmark_cache_backup/n5 .benchmark_cache/n5_xtb
```

Or simply delete and regenerate:

```bash
rm -rf .benchmark_cache
python benchmark.py  # Regenerates with _xtb suffix
```
