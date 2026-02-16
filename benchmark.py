"""Benchmark: effect of n_conformers on model accuracy and timing."""

import argparse
import time

import pandas as pd
from pathlib import Path
from tqdm import tqdm

from rdkit import RDLogger

from deltahf.atom_equivalents import (
    PARAM_NAMES_4,
    PARAM_NAMES_7,
    PARAM_NAMES_EXTENDED,
    PARAM_NAMES_HYBRID,
    fit_atom_equivalents,
    kfold_cross_validation,
    max_abs_deviation,
    mean_abs_deviation,
    r_squared,
    rmsd,
)
from deltahf.cache import ResultCache
from deltahf.data import load_training_data
from deltahf.pipeline import process_molecule

RDLogger.logger().setLevel(RDLogger.ERROR)

N_CONFORMERS = [1, 3, 5]
KFOLD = 10

MODEL_DEFS = {
    "4param": (PARAM_NAMES_4, "atom_counts_4param"),
    "7param": (PARAM_NAMES_7, "atom_counts_7param"),
    "hybrid": (PARAM_NAMES_HYBRID, "atom_counts_hybrid"),
    "extended": (PARAM_NAMES_EXTENDED, "atom_counts_extended"),
}


def run_benchmark(use_gxtb: bool = False):
    df = load_training_data()
    print(f"Loaded {len(df)} molecules")
    if use_gxtb:
        print("Using gxtb single-point energies (wB97M-V/def2-TZVPPD approximation)\n")
    else:
        print()

    benchmark_rows = []

    for n_conf in N_CONFORMERS:
        print(f"{'=' * 60}")
        print(f"  n_conformers = {n_conf}")
        print(f"{'=' * 60}")

        # Use separate cache directories for xTB vs gxtb
        method_suffix = "_gxtb" if use_gxtb else "_xtb"
        cache = ResultCache(Path(f".benchmark_cache/n{n_conf}{method_suffix}"))

        t0 = time.time()
        results = []
        errors = []
        for _, row in tqdm(df.iterrows(), total=len(df), desc=f"  n_conf={n_conf}"):
            result = process_molecule(
                row["smiles"],
                n_conformers=n_conf,
                name=row.get("name", None),
                use_gxtb=use_gxtb,
                cache=cache,
            )
            results.append(result)
            if result.error:
                errors.append(f"{row['name']}: {result.error}")
        elapsed = time.time() - t0

        successful = [(i, r) for i, r in enumerate(results) if r.error is None]
        n_ok = len(successful)
        print(f"  {n_ok}/{len(df)} succeeded in {elapsed:.1f}s")
        if errors:
            print(f"  {len(errors)} failed: {errors[0]}")

        indices = [i for i, _ in successful]
        # Use gxtb energy if available, otherwise fall back to xtb energy
        u_values = [
            (r.gxtb_energy_kcal if r.gxtb_energy_kcal is not None else r.xtb_energy_kcal)
            for _, r in successful
        ]
        exp_dhf = [df.iloc[i]["exp_dhf_kcal_mol"] for i in indices]

        for model_name, (param_names, counts_attr) in MODEL_DEFS.items():
            counts = [getattr(results[i], counts_attr) for i in indices]

            # Dynamic parameter filtering (drop zero-count params)
            totals = {k: sum(c.get(k, 0) for c in counts) for k in param_names}
            active_params = [k for k in param_names if totals[k] > 0]

            epsilon = fit_atom_equivalents(counts, u_values, exp_dhf, active_params)
            predicted = [
                u_values[j] - sum(counts[j].get(k, 0) * epsilon[k] for k in epsilon)
                for j in range(len(u_values))
            ]
            cv = kfold_cross_validation(counts, u_values, exp_dhf, active_params, k=KFOLD)

            row_data = {
                "n_conformers": n_conf,
                "model": model_name,
                "n_params": len(active_params),
                "adj_r2": r_squared(predicted, exp_dhf, p=len(active_params)),
                "rmsd": rmsd(predicted, exp_dhf),
                "mad": mean_abs_deviation(predicted, exp_dhf),
                "max_dev": max_abs_deviation(predicted, exp_dhf),
                "cv_rmsd": cv["cv_rmsd"],
                "n_molecules": n_ok,
                "time_sec": elapsed,
            }
            benchmark_rows.append(row_data)

        print()

    results_df = pd.DataFrame(benchmark_rows)

    # Print summary table
    print(f"\n{'=' * 90}")
    print("  BENCHMARK RESULTS")
    print(f"{'=' * 90}")
    fmt = "{:>12s}  {:>10s}  {:>8s}  {:>8s}  {:>8s}  {:>8s}  {:>8s}  {:>8s}  {:>8s}"
    print(fmt.format(
        "n_conformers", "model", "n_params", "adj_r2", "rmsd", "mad",
        "max_dev", "cv_rmsd", "time(s)",
    ))
    print("  " + "-" * 86)
    for _, row in results_df.iterrows():
        print(fmt.format(
            str(row["n_conformers"]),
            row["model"],
            str(row["n_params"]),
            f"{row['adj_r2']:.4f}",
            f"{row['rmsd']:.2f}",
            f"{row['mad']:.2f}",
            f"{row['max_dev']:.2f}",
            f"{row['cv_rmsd']:.2f}",
            f"{row['time_sec']:.1f}",
        ))
    print()

    results_df.to_csv("benchmark_results.csv", index=False)
    print(f"Results saved to benchmark_results.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark deltahf models with varying n_conformers")
    parser.add_argument(
        "--use-gxtb", action="store_true",
        help="Use gxtb single-point energies (wB97M-V/def2-TZVPPD approximation)",
    )
    args = parser.parse_args()
    run_benchmark(use_gxtb=args.use_gxtb)
