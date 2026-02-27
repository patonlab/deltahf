# Model Benchmarks

This directory contains the benchmark script and results CSV from evaluating `deltahf` across all eight atom equivalent models, three model chemistries (xTB, gXTB, UMA), and three conformer counts (1, 3, 5).

---

## Commands

All commands were run from the `model_benchmarks/` directory.

```bash
# xTB (CPU, 8 threads)
python benchmark.py --methods xtb --xtb-threads 8

# gXTB (CPU, 8 threads) — appends to existing CSV
python benchmark.py --methods gxtb --append --xtb-threads 8

# UMA (GPU) — appends to existing CSV
python benchmark.py --methods uma --append

# RF baseline (no quantum chemistry) — appends to existing CSV
python benchmark.py --methods rf --append
```

Approximate wall times (531 molecules, n_conformers = 1/3/5):

| Method | n=1 | n=3 | n=5 | Hardware |
|--------|----:|----:|----:|----------|
| xTB | 54 s | 104 s | 151 s | CPU (8 threads) |
| gXTB | 69 s | 119 s | 166 s | CPU (8 threads) |
| UMA | 285 s | 707 s | 1095 s | GPU (CUDA) |
| RF | <1 s | — | — | CPU |

---

## Dataset

**Input:** `deltahf/data/training_data.csv` — 531 molecules (C, H, N, O, F, S, Cl).

**Successfully processed:** 526 / 531 molecules (5 failed: 3 conformer generation failures, 1 sulfur compound, 1 isomerization).

Two evaluation subsets are reported:

| Subset | Molecules | Description |
|--------|:---------:|-------------|
| `all` | 526 | Full training set (all successfully processed molecules) |
| `cawkwell` | 102 | Cawkwell2021 subset only (enables comparison with published DFT-B baseline) |

**Literature baseline:** Cawkwell et al. (2021) reported `element` RMSD = 7.59, `element_bo` RMSD = 6.08 kcal/mol using DFT-B (DFTB+ geometries and energies).

---

## Results: Full Training Set (526 molecules, n_conformers = 1)

### xTB

| Model | Params | Adj. R² | RMSD | MAD | Max dev | CV RMSD |
|-------|-------:|--------:|-----:|----:|--------:|--------:|
| element | 7 | 0.9533 | 11.54 | 8.16 | 71.76 | 11.87 |
| element_bo | 11 | 0.9706 | 9.12 | 5.78 | 68.70 | 9.59 |
| hybrid | 13 | 0.9717 | 8.92 | 5.50 | 52.96 | 9.59 |
| bondorder | 13 | 0.9793 | 7.63 | 5.00 | 39.88 | 8.18 |
| **bondorder_ext** | **19** | **0.9847** | **6.53** | **4.26** | **39.40** | **7.09** |
| bondorder_ar | 16 | 0.9794 | 7.60 | 4.96 | 39.96 | 8.20 |
| extended | 18 | 0.9744 | 8.44 | 5.37 | 52.66 | 9.04 |
| neighbour | 28 | 0.9778 | 7.79 | 5.00 | 46.03 | 8.79 |

### gXTB

| Model | Params | Adj. R² | RMSD | MAD | Max dev | CV RMSD |
|-------|-------:|--------:|-----:|----:|--------:|--------:|
| element | 7 | 0.9933 | 4.37 | 3.18 | 22.56 | 4.47 |
| element_bo | 11 | 0.9941 | 4.07 | 2.81 | 21.25 | 4.19 |
| hybrid | 13 | 0.9940 | 4.11 | 2.88 | 20.31 | 4.31 |
| bondorder | 13 | 0.9945 | 3.95 | 2.75 | 21.41 | 4.17 |
| **bondorder_ext** | **19** | **0.9948** | **3.82** | **2.67** | **19.99** | **4.11** |
| bondorder_ar | 16 | 0.9945 | 3.94 | 2.74 | 21.45 | 4.17 |
| extended | 18 | 0.9943 | 3.98 | 2.76 | 20.37 | 4.23 |
| neighbour | 28 | 0.9949 | 3.75 | 2.62 | 24.48 | 4.27 |

### UMA

| Model | Params | Adj. R² | RMSD | MAD | Max dev | CV RMSD |
|-------|-------:|--------:|-----:|----:|--------:|--------:|
| element | 7 | 0.9962 | 3.27 | 2.32 | 21.27 | 3.37 |
| element_bo | 11 | 0.9968 | 2.99 | 1.98 | 21.00 | 3.12 |
| hybrid | 13 | 0.9972 | 2.78 | 1.74 | 20.89 | 2.93 |
| bondorder | 13 | 0.9974 | 2.73 | 1.71 | 20.10 | 2.88 |
| **bondorder_ext** | **19** | **0.9976** | **2.59** | **1.66** | **20.06** | **2.78** |
| bondorder_ar | 16 | 0.9974 | 2.71 | 1.67 | 20.28 | 2.87 |
| extended | 18 | 0.9974 | 2.69 | 1.69 | 20.89 | 2.88 |
| neighbour | 28 | 0.9976 | 2.56 | 1.61 | 21.32 | 2.85 |

### Cross-method comparison (bondorder_ext, best overall)

| Method | RMSD | MAD | Max dev | CV RMSD | Adj. R² |
|--------|-----:|----:|--------:|--------:|--------:|
| xTB | 6.53 | 4.26 | 39.40 | 7.09 | 0.9847 |
| gXTB | 3.82 | 2.67 | 19.99 | 4.11 | 0.9948 |
| UMA | 2.59 | 1.66 | 20.06 | 2.78 | 0.9976 |

---

## Results: Cawkwell2021 Subset (102 molecules, n_conformers = 1)

This subset enables direct comparison with the published DFT-B results from Cawkwell et al. (2021).

### Comparison with literature baseline

| Method | Model | Params | RMSD | Max dev |
|--------|-------|-------:|-----:|--------:|
| DFT-B (lit.) | element | 4 | 7.59 | 25.48 |
| DFT-B (lit.) | element_bo | 7 | 6.08 | 15.01 |
| xTB | bondorder_ext | 15 | 6.91 | 23.66 |
| gXTB | bondorder_ext | 15 | 3.75 | 10.52 |
| **UMA** | **bondorder_ext** | **15** | **2.70** | **8.95** |

### Full model breakdown (Cawkwell subset)

| Model | xTB RMSD | gXTB RMSD | UMA RMSD |
|-------|:--------:|:---------:|:--------:|
| element | 16.97 | 5.31 | 3.91 |
| element_bo | 13.64 | 4.30 | 3.16 |
| hybrid | 12.23 | 4.73 | 3.18 |
| bondorder | 8.83 | 3.90 | 2.95 |
| bondorder_ext | 6.91 | 3.75 | 2.70 |
| bondorder_ar | 8.75 | 3.85 | 2.93 |
| extended | 11.84 | 4.59 | 3.15 |
| neighbour | 7.48 | 3.56 | 2.86 |

---

## Results: Effect of n_conformers

`bondorder_ext` model, full training set (526 molecules). All timings from fresh (uncached) runs.

### RMSD (kcal/mol)

| n_conformers | xTB | gXTB | UMA |
|:---:|:---:|:---:|:---:|
| 1 | 6.53 | 3.82 | 2.59 |
| 3 | 6.51 | 3.81 | 2.58 |
| 5 | 6.52 | 3.83 | 3.11* |

### Max deviation (kcal/mol)

| n_conformers | xTB | gXTB | UMA |
|:---:|:---:|:---:|:---:|
| 1 | 39.40 | 19.99 | 20.06 |
| 3 | 39.35 | 19.99 | 20.05 |
| 5 | 39.33 | 20.02 | 40.85* |

\*UMA n=5 shows degraded accuracy. The max deviation doubles from ~20 to ~41 kcal/mol, and RMSD increases from 2.59 to 3.11. This is likely because UMA's more flexible potential energy surface finds lower-energy conformer geometries that have actually isomerized (changed connectivity) for some molecules, but pass the connectivity check. The xTB and gXTB optimizers do not exhibit this behaviour.

**Recommendation:** Use `--n-conformers 1` (the default) for all methods.

---

## RF Baseline (No Quantum Chemistry)

A Morgan fingerprint + Random Forest baseline was evaluated to quantify the value of the physics-based xTB energy. This model uses no quantum chemistry — only 2D molecular structure.

| Subset | Molecules | In-sample RMSD | CV RMSD |
|--------|:---------:|:--------------:|:-------:|
| all | 531 | 11.70 | 30.92 |
| cawkwell | 102 | 11.56 | 31.13 |

The RF baseline's CV RMSD (30.9 kcal/mol) is ~2.6× worse than even the simplest xTB atom equivalent model (`element`, CV RMSD = 11.9 kcal/mol) and ~11× worse than UMA + `bondorder_ext` (CV RMSD = 2.78 kcal/mol). This demonstrates that the xTB/gXTB/UMA energy provides an essential physics prior that simple cheminformatics descriptors cannot replace.

---

## Notes on the `neighbour` model

The `neighbour` model achieves competitive training-set RMSD (often the best or near-best for each method), but its cross-validation RMSD is unreliable on the smaller Cawkwell subset: CV RMSD values of 994 (xTB), 18665 (gXTB), and 18664 (UMA) indicate a near-singular design matrix in some CV folds. On the full 526-molecule training set the CV RMSD is reasonable (8.79 for xTB, 4.27 for gXTB, 2.85 for UMA), but the instability on smaller datasets makes it unsuitable for general use. The `bondorder_ext` model is recommended instead.

---

## Discussion: Best Model Choice

The `bondorder_ext` model consistently achieves the best or near-best accuracy across all three methods:

- **xTB:** bondorder_ext is clearly best (RMSD 6.53 vs next-best bondorder at 7.63)
- **gXTB:** bondorder_ext is best by a narrow margin (RMSD 3.82 vs neighbour at 3.75, but neighbour has CV instability)
- **UMA:** bondorder_ext is best among stable models (RMSD 2.59 vs neighbour at 2.56, but neighbour has CV instability on smaller datasets)

The key advantage of bond-order classification over hybridization is that it distinguishes atoms by their bonding environment more finely. For example, a carbonyl carbon (C=O, bond order 2) is classified differently from an aromatic carbon (also sp2 in RDKit), and the hydrogen-count extension further differentiates methyl (CH3) from methylene (CH2) from methine (CH) carbons.

For production use, **`bondorder_ext` + UMA** (GPU) gives the best accuracy (RMSD 2.59 kcal/mol, CV RMSD 2.78). On CPU, **`bondorder_ext` + gXTB** (RMSD 3.82, CV RMSD 4.11) is the best option, and **`bondorder_ext` + xTB** (RMSD 6.53, CV RMSD 7.09) is the most accessible with no optional dependencies.

---

## Files

| File | Description |
|------|-------------|
| `benchmark.py` | Benchmark script |
| `benchmark_results.csv` | Full results (149 rows: 8 models × 3 methods × 3 n_conformers × 2 subsets + DFT-B + RF) |
