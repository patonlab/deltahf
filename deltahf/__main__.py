"""CLI entry point for deltahf."""

import argparse
import json
import sys
from pathlib import Path

from deltahf.atom_equivalents import (
    PARAM_NAMES_4,
    PARAM_NAMES_7,
    fit_atom_equivalents,
    kfold_cross_validation,
    max_abs_deviation,
    rmsd,
)
from deltahf.pipeline import process_csv, process_molecule


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="deltahf",
        description="Estimate gas-phase heats of formation using xTB and atom equivalent energies",
    )
    subparsers = parser.add_subparsers(dest="command")

    # --- fit subcommand ---
    fit_parser = subparsers.add_parser("fit", help="Fit atom equivalent energies to training data")
    fit_parser.add_argument("--input", "-i", required=True, help="CSV with smiles, exp_dhf_kcal_mol columns")
    fit_parser.add_argument("--model", choices=["4param", "7param", "both"], default="both")
    fit_parser.add_argument("--kfold", type=int, default=10, help="Number of CV folds")
    fit_parser.add_argument("--n-conformers", type=int, default=5, help="Number of conformers to optimize")
    fit_parser.add_argument("--output", "-o", help="Output JSON file for fitted epsilon values")
    fit_parser.add_argument("--csv", help="Output CSV with training data and predictions")

    # --- predict subcommand ---
    pred_parser = subparsers.add_parser("predict", help="Predict DeltaHf for new molecules")
    pred_parser.add_argument("--input", "-i", required=True, help="CSV with smiles column")
    pred_parser.add_argument("--epsilon", required=True, help="JSON file with fitted epsilon values")
    pred_parser.add_argument("--model", choices=["4param", "7param"], default="4param")
    pred_parser.add_argument("--n-conformers", type=int, default=5, help="Number of conformers to optimize")
    pred_parser.add_argument("--output", "-o", help="Output CSV with results")

    return parser


def cmd_fit(args):
    """Run the fitting workflow."""
    import pandas as pd

    from deltahf.xtb import find_xtb_binary

    try:
        find_xtb_binary()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} molecules from {args.input}")

    # Process all molecules to get xTB energies
    print("Running xTB optimizations...")
    results = []
    for idx, row in df.iterrows():
        smiles = row["smiles"]
        name = row.get("name", f"mol_{idx}")
        print(f"  [{idx + 1}/{len(df)}] {name}: {smiles}")
        result = process_molecule(smiles, n_conformers=args.n_conformers, name=name)
        if result.error:
            print(f"    ERROR: {result.error}")
        else:
            print(f"    xTB energy: {result.xtb_energy:.6f} Eh")
        results.append(result)

    # Filter successful results
    successful = [(i, r) for i, r in enumerate(results) if r.error is None]
    print(f"\n{len(successful)}/{len(results)} molecules completed successfully")

    if not successful:
        failed = [(i, r) for i, r in enumerate(results) if r.error is not None]
        print("\nAll molecules failed. First error:")
        print(f"  {df.iloc[failed[0][0]]['name']}: {failed[0][1].error}")
        print("\nIs xTB on your PATH? Try: conda activate deltahf")
        sys.exit(1)

    indices = [i for i, _ in successful]
    u_values = [r.xtb_energy_kcal for _, r in successful]
    exp_dhf = [df.iloc[i]["exp_dhf_kcal_mol"] for i in indices]

    output = {}

    if args.model in ("4param", "both"):
        counts_4 = [results[i].atom_counts_4param for i in indices]
        epsilon_4 = fit_atom_equivalents(counts_4, u_values, exp_dhf, PARAM_NAMES_4)
        print("\n4-parameter atom equivalents (kcal/mol):")
        for k, v in epsilon_4.items():
            print(f"  {k}: {v:.6f}")

        predicted_4 = [
            u_values[j] - sum(counts_4[j].get(k, 0) * epsilon_4[k] for k in epsilon_4)
            for j in range(len(u_values))
        ]
        print(f"  RMSD: {rmsd(predicted_4, exp_dhf):.2f} kcal/mol")
        print(f"  Max deviation: {max_abs_deviation(predicted_4, exp_dhf):.2f} kcal/mol")

        cv_4 = kfold_cross_validation(counts_4, u_values, exp_dhf, PARAM_NAMES_4, k=args.kfold)
        print(f"  CV error ({args.kfold}-fold): {cv_4['cv_error']:.2f} (kcal/mol)^2")
        output["4param"] = epsilon_4

    if args.model in ("7param", "both"):
        counts_7 = [results[i].atom_counts_7param for i in indices]
        epsilon_7 = fit_atom_equivalents(counts_7, u_values, exp_dhf, PARAM_NAMES_7)
        print("\n7-parameter atom equivalents (kcal/mol):")
        for k, v in epsilon_7.items():
            print(f"  {k}: {v:.6f}")

        predicted_7 = [
            u_values[j] - sum(counts_7[j].get(k, 0) * epsilon_7[k] for k in epsilon_7)
            for j in range(len(u_values))
        ]
        print(f"  RMSD: {rmsd(predicted_7, exp_dhf):.2f} kcal/mol")
        print(f"  Max deviation: {max_abs_deviation(predicted_7, exp_dhf):.2f} kcal/mol")

        cv_7 = kfold_cross_validation(counts_7, u_values, exp_dhf, PARAM_NAMES_7, k=args.kfold)
        print(f"  CV error ({args.kfold}-fold): {cv_7['cv_error']:.2f} (kcal/mol)^2")
        output["7param"] = epsilon_7

    if args.output:
        with open(args.output, "w") as f:
            json.dump(output, f, indent=2)
        print(f"\nFitted parameters saved to {args.output}")

    if args.csv:
        results_df = df.copy()
        results_df["xtb_energy_eh"] = [
            results[i].xtb_energy if results[i].error is None else None for i in range(len(df))
        ]
        results_df["xtb_energy_kcal_mol"] = [
            results[i].xtb_energy_kcal if results[i].error is None else None for i in range(len(df))
        ]

        # Map predictions back to full DataFrame (None for failed molecules)
        pred_4_full = [None] * len(df)
        pred_7_full = [None] * len(df)
        for j, idx in enumerate(indices):
            if args.model in ("4param", "both"):
                pred_4_full[idx] = predicted_4[j]
            if args.model in ("7param", "both"):
                pred_7_full[idx] = predicted_7[j]

        if args.model in ("4param", "both"):
            results_df["pred_dhf_4param"] = pred_4_full
            results_df["error_4param"] = results_df["pred_dhf_4param"] - results_df["exp_dhf_kcal_mol"]
        if args.model in ("7param", "both"):
            results_df["pred_dhf_7param"] = pred_7_full
            results_df["error_7param"] = results_df["pred_dhf_7param"] - results_df["exp_dhf_kcal_mol"]

        results_df["error"] = [results[i].error for i in range(len(df))]
        results_df.to_csv(args.csv, index=False)
        print(f"Results CSV saved to {args.csv}")


def cmd_predict(args):
    """Run the prediction workflow."""
    from deltahf.xtb import find_xtb_binary

    try:
        find_xtb_binary()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    with open(args.epsilon) as f:
        epsilon_data = json.load(f)

    epsilon = epsilon_data.get(args.model, epsilon_data)

    results_df = process_csv(
        Path(args.input),
        n_conformers=args.n_conformers,
        epsilon_4param=epsilon if args.model == "4param" else None,
        epsilon_7param=epsilon if args.model == "7param" else None,
        output_path=Path(args.output) if args.output else None,
    )

    print(results_df.to_string(index=False))
    if args.output:
        print(f"\nResults saved to {args.output}")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "fit":
        cmd_fit(args)
    elif args.command == "predict":
        cmd_predict(args)


if __name__ == "__main__":
    main()
