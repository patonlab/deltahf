"""CLI entry point for deltahf."""

import argparse
import json
import sys
from pathlib import Path

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
from deltahf.pipeline import process_csv, process_molecule

# Map model name -> (PARAM_NAMES, atom_counts field on MoleculeResult)
MODEL_DEFS = {
    "4param": (PARAM_NAMES_4, "atom_counts_4param"),
    "7param": (PARAM_NAMES_7, "atom_counts_7param"),
    "hybrid": (PARAM_NAMES_HYBRID, "atom_counts_hybrid"),
    "extended": (PARAM_NAMES_EXTENDED, "atom_counts_extended"),
}

# Which models each --model choice expands to
MODEL_GROUPS = {
    "4param": ["4param"],
    "7param": ["7param"],
    "hybrid": ["hybrid"],
    "extended": ["extended"],
    "both": ["4param", "7param"],
    "all": ["4param", "7param", "hybrid", "extended"],
}

BANNER = (
    r"    ____  _  _  ____" "\n"
    r"   (    \/ )( \(  __)" "\n"
    r"    ) D () __ ( ) _)" "\n"
    r"   (____/\_)(_/(__)" "\n"
    "     ~~ deltahf ~~\n"
)

CITATIONS = """\
   Data sources:
      [1] Cawkwell, M. J.; Manner, V. W.; Kress, J. D. J. Chem. Inf. Model. 2021, 61, 3337-3347.
          DOI: 10.1021/acs.jcim.1c00312
      [2] Yalamanchi, K. K.; Monge-Palacios, M.; van Oudenhoven, V. C. O.; Gao, X.; Sarathy, S. M.
          J. Phys. Chem. A 2020, 124, 6270-6283. DOI: 10.1021/acs.jpca.0c02785"""

SEP = "   " + "\u2500" * 50


def _print_model_results(counts, u_values, exp_dhf, param_names, kfold):
    """Print fitting results for one model with formatted tables."""
    totals = {k: sum(c.get(k, 0) for c in counts) for k in param_names}

    # Drop parameters with zero counts (no training examples)
    active_params = [k for k in param_names if totals[k] > 0]
    dropped = [k for k in param_names if totals[k] == 0]

    epsilon = fit_atom_equivalents(counts, u_values, exp_dhf, active_params)

    predicted = [
        u_values[j] - sum(counts[j].get(k, 0) * epsilon[k] for k in epsilon)
        for j in range(len(u_values))
    ]
    cv = kfold_cross_validation(counts, u_values, exp_dhf, active_params, k=kfold)

    n_params = len(active_params)
    print()
    print(SEP)
    print(f"   {n_params}-parameter model ({len(counts)} molecules)")
    print(SEP)

    # Atom environment totals (only active)
    env_parts = ", ".join(f"{k}={totals[k]}" for k in active_params)
    print(f"   Atom environments: {env_parts}")
    if dropped:
        print(f"   Unused parameters (no examples): {', '.join(dropped)}")
    print()

    # Atom equivalents table
    std_eps = cv["std_epsilon"]
    pm = "\u00b1"
    w = 50  # table width (matches SEP)
    print("   Fitted atom equivalents:")
    print(f"   {'Param':<14}{'Epsilon (kcal/mol)':>20}{'Std Err':>16}")
    print(f"   {'─' * 14}{'─' * 20} {'─' * 15}")
    for k, v in epsilon.items():
        print(f"   {k:<14}{v:>20.6f}      {pm} {std_eps[k]:>8.4f}")
    print()

    # Statistics table
    r2_label = "Adj. R\u00b2"
    cv_label = f"CV RMSD ({kfold}-fold)"
    print("   Fit statistics:")
    print(f"   {'Metric':<{w - 5}}{'Value':>5}")
    print(f"   {'─' * w}")
    print(f"   {r2_label:<35}{r_squared(predicted, exp_dhf, p=len(param_names)):>15.4f}")
    print(f"   {'RMSD':<35}{rmsd(predicted, exp_dhf):>15.2f} kcal/mol")
    print(f"   {'MAD':<35}{mean_abs_deviation(predicted, exp_dhf):>15.2f} kcal/mol")
    print(f"   {'Max deviation':<35}{max_abs_deviation(predicted, exp_dhf):>15.2f} kcal/mol")
    print(f"   {cv_label:<35}{cv['cv_rmsd']:>15.2f} kcal/mol")

    return epsilon, predicted


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="deltahf",
        description="Estimate gas-phase heats of formation using xTB and atom equivalent energies",
    )
    subparsers = parser.add_subparsers(dest="command")

    # --- fit subcommand ---
    fit_parser = subparsers.add_parser("fit", help="Fit atom equivalent energies to training data")
    fit_parser.add_argument("--input", "-i", required=True, help="CSV with smiles, exp_dhf_kcal_mol columns")
    fit_parser.add_argument(
        "--model", choices=["4param", "7param", "hybrid", "extended", "both", "all"], default="both",
    )
    fit_parser.add_argument("--kfold", type=int, default=10, help="Number of CV folds")
    fit_parser.add_argument("--n-conformers", type=int, default=1, help="Number of conformers to optimize")
    fit_parser.add_argument("--output", "-o", help="Output JSON file for fitted epsilon values")
    fit_parser.add_argument("--csv", help="Output CSV with training data and predictions")
    fit_parser.add_argument(
        "--use-xtb-wbos", action="store_true",
        help="Use xTB Wiberg bond orders (instead of RDKit) for 7-param atom classification",
    )
    fit_parser.add_argument(
        "--use-gxtb", action="store_true",
        help="Use gxtb single-point energies (wB97M-V/def2-TZVPPD). WARNING: Do not mix with xTB workflows!",
    )
    fit_parser.add_argument(
        "--cache-dir", type=str, default=None,
        help="Directory for caching xTB results (enables restart capability)",
    )
    fit_parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print per-molecule details instead of a progress bar",
    )

    # --- predict subcommand ---
    pred_parser = subparsers.add_parser("predict", help="Predict DeltaHf for new molecules")
    pred_parser.add_argument("--input", "-i", required=True, help="CSV with smiles column")
    pred_parser.add_argument("--epsilon", required=False, help="JSON file with fitted epsilon values (uses defaults if not specified)")
    pred_parser.add_argument("--model", choices=["4param", "7param", "hybrid", "extended"], default="4param")
    pred_parser.add_argument("--n-conformers", type=int, default=1, help="Number of conformers to optimize")
    pred_parser.add_argument("--output", "-o", help="Output CSV with results")
    pred_parser.add_argument(
        "--use-xtb-wbos", action="store_true",
        help="Use xTB Wiberg bond orders (instead of RDKit) for 7-param atom classification",
    )
    pred_parser.add_argument(
        "--use-gxtb", action="store_true",
        help="Use gxtb single-point energies (wB97M-V/def2-TZVPPD). Must match method used in fit!",
    )
    pred_parser.add_argument(
        "--cache-dir", type=str, default=None,
        help="Directory for caching xTB results (enables restart capability)",
    )
    pred_parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print per-molecule details instead of a progress bar",
    )

    return parser


def cmd_fit(args):
    """Run the fitting workflow."""
    import pandas as pd

    from deltahf.xtb import find_xtb_binary

    try:
        find_xtb_binary()
    except FileNotFoundError as e:
        print(f"   Error: {e}")
        sys.exit(1)

    if args.use_gxtb:
        from deltahf.xtb import find_gxtb_binary
        if find_gxtb_binary() is None:
            print("   Error: gxtb binary not found on PATH.")
            print("   gxtb must be installed manually from source (not available via conda/pip).")
            sys.exit(1)
        print("   Using gxtb energies (wB97M-V/def2-TZVPPD approximation)")
        print("   WARNING: gxtb and xTB energies are on different scales - do not mix workflows!")
        print()

    print(CITATIONS)
    print()

    df = pd.read_csv(args.input)
    print(f"   Loaded {len(df)} molecules from {args.input}")
    if args.use_xtb_wbos:
        print("   Using xTB Wiberg bond orders for 7-param atom classification")

    cache = None
    if args.cache_dir:
        from deltahf.cache import ResultCache

        # Append method suffix to cache directory to keep xTB and gxtb separate
        cache_path = Path(args.cache_dir)
        method_suffix = "_gxtb" if args.use_gxtb else "_xtb"
        if not cache_path.name.endswith(method_suffix):
            cache_path = cache_path.parent / (cache_path.name + method_suffix)

        cache = ResultCache(cache_path)
        print(f"   Using cache directory: {cache_path}")
        if args.use_xtb_wbos:
            print("      Note: caching disabled when --use-xtb-wbos is set")
            cache = None

    # Process all molecules to get xTB energies
    print("   Running xTB optimizations...")
    results = []
    errors = []
    warnings = []

    if args.verbose:
        for idx, row in df.iterrows():
            smiles = row["smiles"]
            name = row.get("name", f"mol_{idx}")
            print(f"      [{idx + 1}/{len(df)}] {name}: {smiles}")
            result = process_molecule(
                smiles, n_conformers=args.n_conformers, name=name,
                use_xtb_wbos=args.use_xtb_wbos, use_gxtb=args.use_gxtb, cache=cache,
            )
            if result.error:
                print(f"         ERROR: {result.error}")
            else:
                print(f"         xTB energy: {result.xtb_energy:.6f} Eh")
                if result.gxtb_energy is not None:
                    print(f"         gxtb energy: {result.gxtb_energy:.6f} Eh")
                if result.n_conformers_isomerized > 0:
                    print(f"         WARNING: {result.n_conformers_isomerized} conformer(s) isomerized")
            results.append(result)
    else:
        from tqdm import tqdm

        for idx, row in tqdm(df.iterrows(), total=len(df), desc="   Molecules"):
            smiles = row["smiles"]
            name = row.get("name", f"mol_{idx}")
            result = process_molecule(
                smiles, n_conformers=args.n_conformers, name=name,
                use_xtb_wbos=args.use_xtb_wbos, use_gxtb=args.use_gxtb, cache=cache,
            )
            if result.error:
                errors.append(f"{name}: {result.error}")
            elif result.n_conformers_isomerized > 0:
                warnings.append(f"{name}: {result.n_conformers_isomerized} conformer(s) isomerized")
            results.append(result)

        if errors:
            print(f"\n   {len(errors)} molecule(s) failed:")
            for msg in errors:
                print(f"      {msg}")
        if warnings:
            print(f"\n   {len(warnings)} molecule(s) had isomerized conformers:")
            for msg in warnings:
                print(f"      {msg}")

    # Filter successful results
    successful = [(i, r) for i, r in enumerate(results) if r.error is None]
    print(f"\n   {len(successful)}/{len(results)} molecules completed successfully")

    if not successful:
        failed = [(i, r) for i, r in enumerate(results) if r.error is not None]
        print("\n   All molecules failed. First error:")
        print(f"      {df.iloc[failed[0][0]]['name']}: {failed[0][1].error}")
        print("\n   Is xTB on your PATH? Try: conda activate deltahf")
        sys.exit(1)

    indices = [i for i, _ in successful]
    # Use gxtb energy if available, otherwise fall back to xtb energy
    u_values = [
        (r.gxtb_energy_kcal if r.gxtb_energy_kcal is not None else r.xtb_energy_kcal)
        for _, r in successful
    ]
    exp_dhf = [df.iloc[i]["exp_dhf_kcal_mol"] for i in indices]

    output = {}
    predicted = {}
    models_to_run = MODEL_GROUPS[args.model]

    for model_name in models_to_run:
        param_names, counts_attr = MODEL_DEFS[model_name]
        counts = [getattr(results[i], counts_attr) for i in indices]
        eps, pred = _print_model_results(counts, u_values, exp_dhf, param_names, args.kfold)
        output[model_name] = eps
        predicted[model_name] = pred

    # Add metadata to track which method was used for fitting
    output["_metadata"] = {
        "method": "gxtb" if args.use_gxtb else "xtb",
        "n_conformers": args.n_conformers,
        "n_molecules": len(indices),
    }

    print()
    print(SEP)

    if args.output:
        with open(args.output, "w") as f:
            json.dump(output, f, indent=2)
        print(f"   Fitted parameters saved to {args.output}")
        print(f"   Method: {'gxtb' if args.use_gxtb else 'xtb'}")

    if args.csv:
        results_df = df.copy()
        results_df["xtb_energy_eh"] = [
            results[i].xtb_energy if results[i].error is None else None for i in range(len(df))
        ]
        results_df["xtb_energy_kcal_mol"] = [
            results[i].xtb_energy_kcal if results[i].error is None else None for i in range(len(df))
        ]

        # Map predictions back to full DataFrame (None for failed molecules)
        for model_name in models_to_run:
            pred_full = [None] * len(df)
            for j, idx in enumerate(indices):
                pred_full[idx] = predicted[model_name][j]
            results_df[f"pred_dhf_{model_name}"] = pred_full
            results_df[f"error_{model_name}"] = (
                results_df[f"pred_dhf_{model_name}"] - results_df["exp_dhf_kcal_mol"]
            )

        results_df["error"] = [results[i].error for i in range(len(df))]
        results_df.to_csv(args.csv, index=False)
        print(f"   Results CSV saved to {args.csv}")


def cmd_predict(args):
    """Run the prediction workflow."""
    from deltahf.xtb import find_xtb_binary

    try:
        find_xtb_binary()
    except FileNotFoundError as e:
        print(f"   Error: {e}")
        sys.exit(1)

    if args.use_gxtb:
        from deltahf.xtb import find_gxtb_binary
        if find_gxtb_binary() is None:
            print("   Error: gxtb binary not found on PATH.")
            print("   gxtb must be installed manually from source (not available via conda/pip).")
            sys.exit(1)

    # Determine which epsilon file to use
    if args.epsilon:
        epsilon_file = args.epsilon
    else:
        # Use default parameters from params directory
        params_dir = Path(__file__).parent.parent / "params"
        if args.use_gxtb:
            epsilon_file = params_dir / "params_gxtb.json"
        else:
            epsilon_file = params_dir / "params_xtb.json"

        if not epsilon_file.exists():
            print(f"   Error: Default parameter file not found: {epsilon_file}")
            sys.exit(1)

        print(f"   WARNING: Using default parameters from {epsilon_file}")
        print(f"            Specify --epsilon to use custom parameters\n")

    with open(epsilon_file) as f:
        epsilon_data = json.load(f)

    # Check for method mismatch
    metadata = epsilon_data.get("_metadata", {})
    param_method = metadata.get("method", "xtb")  # Default to xtb for old parameter files
    current_method = "gxtb" if args.use_gxtb else "xtb"

    if param_method != current_method:
        print(f"\n   WARNING: Method mismatch detected!")
        print(f"   Parameters were fitted using: {param_method}")
        print(f"   You are predicting using: {current_method}")
        print(f"   This will produce incorrect results - energies are on different scales!")
        print(f"\n   Either:")
        print(f"     - Use {'--use-gxtb' if param_method == 'gxtb' else 'without --use-gxtb'} to match the fitting method")
        print(f"     - Refit parameters using {current_method}\n")
        sys.exit(1)

    epsilon = epsilon_data.get(args.model, epsilon_data)

    cache = None
    if args.cache_dir:
        from deltahf.cache import ResultCache

        # Append method suffix to cache directory to keep xTB and gxtb separate
        cache_path = Path(args.cache_dir)
        method_suffix = "_gxtb" if args.use_gxtb else "_xtb"
        if not cache_path.name.endswith(method_suffix):
            cache_path = cache_path.parent / (cache_path.name + method_suffix)

        cache = ResultCache(cache_path)
        print(f"   Using cache directory: {cache_path}")
        if args.use_xtb_wbos:
            print("      Note: caching disabled when --use-xtb-wbos is set")
            cache = None

    epsilon_kwargs = {f"epsilon_{args.model}": epsilon}
    results_df = process_csv(
        Path(args.input),
        n_conformers=args.n_conformers,
        output_path=Path(args.output) if args.output else None,
        use_xtb_wbos=args.use_xtb_wbos,
        use_gxtb=args.use_gxtb,
        cache=cache,
        verbose=args.verbose,
        **epsilon_kwargs,
    )

    # Create simplified display with just essential columns
    display_cols = ["smiles"]
    if "name" in results_df.columns:
        display_cols.append("name")
    display_cols.append("xtb_energy_kcal")
    if args.use_gxtb:
        display_cols.append("gxtb_energy_kcal")

    # Add the prediction column based on the model
    dhf_col = f"dhf_{args.model}"
    if dhf_col in results_df.columns:
        display_cols.append(dhf_col)

    if "error" in results_df.columns:
        display_cols.append("error")

    # Filter to columns that exist in the dataframe
    display_cols = [col for col in display_cols if col in results_df.columns]
    display_df = results_df[display_cols]

    print(display_df.to_string(index=False))
    if args.output:
        print(f"\n   Results saved to {args.output}")


def main(argv=None):
    from rdkit import RDLogger

    RDLogger.logger().setLevel(RDLogger.ERROR)

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    print(BANNER)

    if args.command == "fit":
        cmd_fit(args)
    elif args.command == "predict":
        cmd_predict(args)


if __name__ == "__main__":
    main()
