# Model Training

This directory contains the parameter files used as defaults by `deltahf`, along with the documentation of how they were produced. Parameter files should be regenerated here whenever the training data or pipeline changes, and then copied to `params/`.

---

## Training set

**Input:** `deltahf/data/training_data.csv` — 531 molecules (C/H/N/O/F/S/Cl)

| Source | Count | Description |
|--------|------:|-------------|
| Cawkwell2021 | 102 | Energetic CHNO molecules |
| Yalamanchi2020 | 211 | Cyclic hydrocarbons |
| ATcT v1.220 | 218 | Small molecules across C/H/N/O/F/S/Cl |
| **Total** | **531** | |

**Successfully processed:** 526 / 531 molecules

### Pipeline failures (5 molecules)

| Name | Reason |
|------|--------|
| (1S,4R,5S)-5-Ethenylbicyclo[2.2.2]oct-2-ene | Conformer generation failed |
| (1S,4R,5R)-5-Ethenylbicyclo[2.2.2]oct-2-ene | Conformer generation failed |
| [2.2]Paracyclophane | Conformer generation failed |
| Sulfur oxide peroxide | Conformer generation failed |
| Azete | Isomerization detected during geometry optimisation |

---

## Fitting commands

All fits used 10-fold cross-validation, 1 conformer per molecule, and all 8 atom-equivalent models.

### xTB

```bash
python -m deltahf fit \
    -i ../deltahf/data/training_data.csv \
    --model all \
    --kfold 10 \
    --n-conformers 1 \
    --xtb-threads 16 \
    -o params_xtb.json
```

Approximate run time: ~55 s (9.5 molecules/s, 16 threads).

### gxtb

```bash
python -m deltahf fit \
    -i ../deltahf/data/training_data.csv \
    --model all \
    --kfold 10 \
    --n-conformers 1 \
    --use-gxtb \
    --xtb-threads 16 \
    -o params_gxtb.json
```

Approximate run time: ~2 min 15 s (3.9 molecules/s). gxtb runs a wB97M-V/def2-TZVPPD single-point on the xTB geometry, making it ~2.4× slower than xTB alone.

### UMA (to be added)

```bash
python -m deltahf fit \
    -i ../deltahf/data/training_data.csv \
    --model all \
    --kfold 10 \
    --n-conformers 1 \
    --optimizer uma \
    -o params_uma.json
```

---

## Performance summary

All energies in kcal/mol. **Params** = number of parameters actually fitted after dropping atom types with zero training examples. **CV RMSD** = 10-fold cross-validation RMSD (leave-10%-out). All metrics are on the 526 successfully processed molecules.

### xTB

| Model | Params | RMSD | MAD | Max dev | Adj. R² | CV RMSD |
|-------|-------:|-----:|----:|--------:|--------:|--------:|
| element | 7 | 11.54 | 8.16 | 71.76 | 0.9533 | 11.87 |
| element_bo | 11 | 9.12 | 5.78 | 68.70 | 0.9706 | 9.59 |
| hybrid | 13 | 8.92 | 5.50 | 52.96 | 0.9716 | 9.59 |
| bondorder | 13 | 7.63 | 5.00 | 39.88 | 0.9792 | 8.18 |
| bondorder_ext | 19 | 6.53 | 4.26 | 39.40 | 0.9846 | 7.09 |
| bondorder_ar | 16 | 7.60 | 4.96 | 39.96 | 0.9792 | 8.20 |
| extended | 18 | 8.44 | 5.37 | 52.66 | 0.9744 | 9.04 |
| neighbour | 28 | 7.79 | 5.00 | 46.03 | 0.9777 | 8.79 |

### gxtb

| Model | Params | RMSD | MAD | Max dev | Adj. R² | CV RMSD |
|-------|-------:|-----:|----:|--------:|--------:|--------:|
| element | 7 | 4.37 | 3.18 | 22.56 | 0.9933 | 4.47 |
| element_bo | 11 | 4.07 | 2.81 | 21.25 | 0.9941 | 4.19 |
| hybrid | 13 | 4.11 | 2.88 | 20.31 | 0.9940 | 4.31 |
| bondorder | 13 | 3.95 | 2.75 | 21.41 | 0.9944 | 4.17 |
| bondorder_ext | 19 | 3.82 | 2.67 | 19.99 | 0.9947 | 4.11 |
| bondorder_ar | 16 | 3.94 | 2.74 | 21.45 | 0.9944 | 4.17 |
| extended | 18 | 3.98 | 2.76 | 20.37 | 0.9943 | 4.23 |
| neighbour | 28 | 3.75 | 2.62 | 24.48 | 0.9948 | 4.27 |

### UMA

*To be added.*

---

## Notes

- gxtb gives substantially better fits than xTB across all models (RMSD roughly 2–3× lower).
- The `bondorder_ext` model gives the best xTB fit (CV RMSD 7.09); the `neighbour` model gives the best gxtb fit (CV RMSD 4.27), though most gxtb models are similar.
- Max deviations are much larger for xTB (~40–72 kcal/mol) than gxtb (~20–24 kcal/mol).
- Parameter counts listed are after dynamic filtering (e.g. O_sp, S_sp, O_3, S_3 are dropped for models where no training examples exist).
