# Analysis of Table 1 from PFAS Paper

## Source
**Paper:** "Accurate Enthalpies of Formation for PFAS from First Principles: Combining Different Levels of Theory in a Generalized Thermochemical Hierarchy"
**Authors:** Kento Abeywardane and C. Franklin Goldsmith
**Journal:** ACS Phys. Chem Au 2024, 4, 247–258

## Summary

Table 1 contains **16 reference species** with formation enthalpies at 0 K from the Active Thermochemical Tables (ATcT).

### Key Findings:

1. **All 8 CHNO-only molecules are already in the training data**
   - H₂, CH₄, C₂H₆ (ethane), C₃H₈ (propane), H₂O, CH₂O, CH₃OH, CH₃OCH₃, CH₃C(O)OH

2. **All 8 fluorinated compounds are currently unsupported**
   - HF, CH₃F, CH₂F₂, CHF₃, CF₄, C₂F₄, CF₂O
   - deltahf only supports C, H, N, O elements
   - Adding fluorine support would require updating SUPPORTED_ELEMENTS and refitting

3. **Temperature difference: 0 K vs 298 K**
   - Table 1 values are at 0 K
   - Training data values are at 298.15 K (standard conditions)
   - The temperature difference causes systematic differences in ΔHf° values
   - 0 K values are typically less negative (or more positive) than 298 K values

## Comparison with Existing Training Data

| Molecule | ATcT (0 K) | Training (298 K) | Difference | Training ID |
|----------|------------|------------------|------------|-------------|
| H₂ | 0.00 | 0.0 | 0.0 | 102 |
| CH₄ | -15.90 | -17.9 | -2.0 | 46 |
| C₂H₆ | -16.35 | -20.1 | -3.8 | 51 |
| C₃H₈ | -19.77 | -25.0 | -5.2 | 66 |
| H₂O | -57.10 | -57.8 | -0.7 | 48 |
| CH₂O | -25.19 | -26.0 | -0.8 | 53 |
| CH₃OH | -45.42 | -48.0 | -2.6 | 54 |
| CH₃OCH₃ | -39.80 | -44.0 | -4.2 | 92 |
| CH₃COOH | -100.03 | -103.4 | -3.4 | 96 |

*All values in kcal/mol*

**Note:** The temperature corrections (0 K → 298 K) are all in the expected direction (more negative at 298 K) and of reasonable magnitude (typically 2-5 kcal/mol for small molecules).

## Recommendations

1. **No new CHNO molecules to add** - All 8 CHNO molecules from Table 1 are already in training_data.csv

2. **Fluorinated compounds cannot be used without expanding element support** - Would require:
   - Updating `SUPPORTED_ELEMENTS` in `deltahf/smiles.py` to include F
   - Refitting all models with fluorine-containing molecules
   - Validating that the atom equivalent method works for fluorinated species
   - Adding fluorine parameters:
     - 4-param: Add F
     - 7-param: Add F and F' (multiply-bonded)
     - hybrid: Add F_sp3, F_sp2, F_sp
     - extended: Similar to hybrid
   - Table 1 provides 7 fluorinated reference species with high-quality ATcT values

3. **Temperature consistency** - Training data should remain at 298 K for consistency with experimental conventions and most thermochemical databases

4. **Potential future work** - The PFAS paper demonstrates that fluorinated compounds can be handled accurately with hierarchical thermochemical methods, suggesting the atom equivalent approach could potentially be extended to fluorine

## Data Quality Note

The ATcT database is extremely high-quality with very small uncertainties (typically ±0.02-0.40 kJ/mol for these species). The 0 K vs 298 K temperature difference is expected and well-understood - the correction involves heat capacity integrations:

ΔH°(298K) = ΔH°(0K) + ∫₀²⁹⁸ Cp(products) dT - ∫₀²⁹⁸ Cp(reactants) dT

For small molecules, this correction is typically 2-5 kcal/mol and makes values more negative (or less positive) at higher temperature.

## Extracted Data

The extracted data has been saved to `table1_extracted.csv` with the following columns:
- `id`, `name`, `formula`, `smiles`: Molecule identifiers
- `exp_dhf_kcal_mol`: Formation enthalpy at 0 K (converted from kJ/mol)
- `source`, `category`: Classification
- `notes`: Temperature info and training data comparison
- `atct_dhf_0K_kJ_mol`, `uncertainty_kJ_mol`: Original ATcT values
- `in_training_data`: Duplicate status with existing training data
- `contains_fluorine`: Element support flag
