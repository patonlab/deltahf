"""Comprehensive benchmark: parameterisation × model chemistry × n_conformers.

Compares all atom-equivalent models (MODEL_DEFS) across:
  - Methods: xTB (GFN2), gXTB (wB97M-V/def2-TZVPPD single-point), UMA (MLIP)
  - n_conformers: 1, 3, 5
  - Subsets: all 314 training molecules AND the 102-molecule Cawkwell2021 subset

Literature baseline rows (DFT-B + xTB geometry, Cawkwell2021, Table 2) are
included in the CSV for direct comparison against the Cawkwell2021 subset.

TIMING NOTE: wall time (time_sec) is only meaningful for fresh runs with an
empty cache. On subsequent runs the cache is hit and elapsed time reflects
IO overhead only, not actual computation.

Usage:
    python benchmark.py                                 # xTB, all conformers
    python benchmark.py --methods gxtb --append         # add gXTB to CSV
    python benchmark.py --methods uma --append          # add UMA to CSV
    python benchmark.py --methods xtb gxtb uma          # all at once
    python benchmark.py --methods xtb --n-conformers 1  # single conformer
"""

import argparse
import time
from pathlib import Path

import pandas as pd
from rdkit import RDLogger
from tqdm import tqdm

from deltahf.atom_equivalents import (
    fit_atom_equivalents,
    kfold_cross_validation,
    max_abs_deviation,
    mean_abs_deviation,
    r_squared,
    rmsd,
)
from deltahf.cache import ResultCache
from deltahf.constants import MODEL_DEFS
from deltahf.data import load_training_data
from deltahf.pipeline import _best_energy_kcal, process_molecule

RDLogger.logger().setLevel(RDLogger.ERROR)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

METHODS = {
    "xtb":  {"optimizer": "xtb", "use_gxtb": False, "label": "xTB",  "hardware": "CPU"},
    "gxtb": {"optimizer": "xtb", "use_gxtb": True,  "label": "gXTB", "hardware": "CPU"},
    "uma":  {"optimizer": "uma", "use_gxtb": False,  "label": "UMA",  "hardware": "GPU"},
}

N_CONFORMERS = [1, 3, 5]
KFOLD = 10
CACHE_BASE = Path(".benchmark_cache")

# Published DFT-B baseline: Cawkwell2021, J. Chem. Inf. Model. 2021, 61, 3337-3347, Table 2
# Results reported on the 103-molecule Cawkwell2021 CHNO set (we have 102 of these).
_NA = float("nan")
LITERATURE_ROWS = [
    {
        "method": "DFT-B (lit)", "n_conformers": _NA, "model": "4param",
        "n_params": 4, "subset": "cawkwell", "n_molecules": 103,
        "adj_r2": _NA, "rmsd": 7.59, "mad": _NA, "max_dev": 25.48,
        "cv_rmsd": _NA, "time_sec": _NA, "hardware": None,
    },
    {
        "method": "DFT-B (lit)", "n_conformers": _NA, "model": "7param",
        "n_params": 7, "subset": "cawkwell", "n_molecules": 103,
        "adj_r2": _NA, "rmsd": 6.08, "mad": _NA, "max_dev": 15.01,
        "cv_rmsd": _NA, "time_sec": _NA, "hardware": None,
    },
]

CSV_COLUMNS = [
    "method", "n_conformers", "model", "n_params", "subset", "n_molecules",
    "adj_r2", "rmsd", "mad", "max_dev", "cv_rmsd", "time_sec", "hardware",
]


# ---------------------------------------------------------------------------
# Core benchmark logic
# ---------------------------------------------------------------------------

def _compute_model_stats(results, df_indices, u_values, exp_dhf, param_names, counts_attr):
    """Fit one model on a subset and return a stats dict."""
    counts = [getattr(results[i], counts_attr) for i in df_indices]
    totals = {k: sum(c.get(k, 0) for c in counts) for k in param_names}
    active_params = [k for k in param_names if totals[k] > 0]

    epsilon = fit_atom_equivalents(counts, u_values, exp_dhf, active_params)
    predicted = [
        u_values[j] - sum(counts[j].get(k, 0) * epsilon[k] for k in epsilon)
        for j in range(len(u_values))
    ]
    cv = kfold_cross_validation(counts, u_values, exp_dhf, active_params, k=KFOLD)

    return {
        "n_params": len(active_params),
        "adj_r2": r_squared(predicted, exp_dhf, p=len(active_params)),
        "rmsd": rmsd(predicted, exp_dhf),
        "mad": mean_abs_deviation(predicted, exp_dhf),
        "max_dev": max_abs_deviation(predicted, exp_dhf),
        "cv_rmsd": cv["cv_rmsd"],
    }


def run_method(method_key: str, n_conformers_list: list[int], df: pd.DataFrame) -> list[dict]:
    """Run all models × all n_conformers × both subsets for one method."""
    cfg = METHODS[method_key]
    optimizer = cfg["optimizer"]
    use_gxtb = cfg["use_gxtb"]
    label = cfg["label"]
    hardware = cfg["hardware"]

    predictor = None
    if optimizer != "xtb":
        from deltahf.uma import load_mlip_model
        print(f"  Loading {optimizer} model...")
        predictor = load_mlip_model(optimizer)
        print()

    cawkwell_mask = (df["source"] == "Cawkwell2021").values

    rows = []

    for n_conf in n_conformers_list:
        print(f"  {'─' * 56}")
        print(f"  {label}  |  n_conformers = {n_conf}")
        print(f"  {'─' * 56}")

        method_parts = [optimizer]
        if use_gxtb:
            method_parts.append("gxtb")
        cache_dir = CACHE_BASE / f"n{n_conf}_{'_'.join(method_parts)}"
        cache = ResultCache(cache_dir)

        t0 = time.time()
        results = []
        for _, row in tqdm(df.iterrows(), total=len(df), desc=f"    {label} n={n_conf}"):
            result = process_molecule(
                row["smiles"],
                n_conformers=n_conf,
                name=row.get("name", None),
                optimizer=optimizer,
                predictor=predictor,
                use_gxtb=use_gxtb,
                cache=cache,
            )
            results.append(result)
        elapsed = time.time() - t0

        subsets = {
            "all":      [i for i, r in enumerate(results) if r.error is None],
            "cawkwell": [i for i, r in enumerate(results) if r.error is None and cawkwell_mask[i]],
        }

        n_ok_all = len(subsets["all"])
        n_ok_caw = len(subsets["cawkwell"])
        print(f"    {n_ok_all}/{len(df)} succeeded ({n_ok_caw} Cawkwell) in {elapsed:.1f}s")

        for subset_label, indices in subsets.items():
            if not indices:
                continue
            u_values = [_best_energy_kcal(results[i]) for i in indices]
            exp_dhf = [df.iloc[i]["exp_dhf_kcal_mol"] for i in indices]

            for model_name, (param_names, counts_attr) in MODEL_DEFS.items():
                stats = _compute_model_stats(
                    results, indices, u_values, exp_dhf, param_names, counts_attr,
                )
                rows.append({
                    "method": label,
                    "n_conformers": n_conf,
                    "model": model_name,
                    "subset": subset_label,
                    "n_molecules": len(indices),
                    "time_sec": elapsed,
                    "hardware": hardware,
                    **stats,
                })

        print()

    return rows


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive benchmark: parameterisation × method × n_conformers",
    )
    parser.add_argument(
        "--methods", nargs="+", choices=list(METHODS.keys()), default=["xtb"],
        help="Which method(s) to run (default: xtb)",
    )
    parser.add_argument(
        "--n-conformers", nargs="+", type=int, default=N_CONFORMERS, metavar="N",
        help="Conformer counts to test (default: 1 3 5)",
    )
    parser.add_argument(
        "--output", default="benchmark_results.csv",
        help="Output CSV filename (default: benchmark_results.csv)",
    )
    parser.add_argument(
        "--append", action="store_true",
        help="Append to existing CSV instead of overwriting",
    )
    args = parser.parse_args()

    df = load_training_data()
    n_caw = (df["source"] == "Cawkwell2021").sum()
    n_yal = (df["source"] == "Yalamanchi2020").sum()
    print(f"Loaded {len(df)} molecules ({n_caw} Cawkwell2021, {n_yal} Yalamanchi2020)")
    print()

    all_rows = []
    for method_key in args.methods:
        print(f"{'=' * 60}")
        print(f"  Method: {METHODS[method_key]['label']}")
        print(f"{'=' * 60}")
        method_rows = run_method(method_key, args.n_conformers, df)
        all_rows.extend(method_rows)

    new_df = pd.DataFrame(all_rows, columns=CSV_COLUMNS)

    lit_df = pd.DataFrame(LITERATURE_ROWS, columns=CSV_COLUMNS)
    if args.append and Path(args.output).exists():
        existing = pd.read_csv(args.output)
        has_lit = (existing["method"] == "DFT-B (lit)").any()
        frames = [existing, new_df] + ([] if has_lit else [lit_df])
        combined = pd.concat(frames, ignore_index=True)
    else:
        combined = pd.concat([new_df, lit_df], ignore_index=True)

    combined.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}  ({len(combined)} rows)")

    # Print summary for n_conformers=1
    n1 = combined[combined["n_conformers"] == 1].copy()
    if not n1.empty:
        print()
        print(f"{'=' * 100}")
        print("  SUMMARY  (n_conformers=1)")
        print(f"{'=' * 100}")
        fmt = "{:<10}  {:<14}  {:<10}  {:>8}  {:>8}  {:>8}  {:>8}  {:>8}"
        print(fmt.format("method", "model", "subset", "adj_r2", "rmsd", "mad", "max_dev", "cv_rmsd"))
        print("  " + "-" * 96)
        for _, row in n1.sort_values(["subset", "method", "rmsd"]).iterrows():
            print(fmt.format(
                str(row["method"]),
                str(row["model"]),
                str(row["subset"]),
                f"{row['adj_r2']:.4f}"  if pd.notna(row["adj_r2"])  else "—",
                f"{row['rmsd']:.2f}"    if pd.notna(row["rmsd"])    else "—",
                f"{row['mad']:.2f}"     if pd.notna(row["mad"])     else "—",
                f"{row['max_dev']:.2f}" if pd.notna(row["max_dev"]) else "—",
                f"{row['cv_rmsd']:.2f}" if pd.notna(row["cv_rmsd"]) else "—",
            ))

    lit = combined[combined["method"] == "DFT-B (lit)"]
    if not lit.empty:
        print()
        print("  Literature baseline (DFT-B, Cawkwell2021 subset):")
        for _, row in lit.iterrows():
            print(f"    {row['model']:10s}  rmsd={row['rmsd']:.2f}  max_dev={row['max_dev']:.2f}")


if __name__ == "__main__":
    main()
