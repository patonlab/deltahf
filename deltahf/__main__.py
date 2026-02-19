"""CLI entry point for deltahf."""

import argparse
import json
import sys
from pathlib import Path

from deltahf.atom_equivalents import (
    fit_atom_equivalents,
    kfold_cross_validation,
    max_abs_deviation,
    mean_abs_deviation,
    r_squared,
    rmsd,
)
from deltahf.constants import MODEL_DEFS
from deltahf.pipeline import _best_energy_kcal, process_csv, process_molecule

# Which models each --model choice expands to
MODEL_GROUPS = {
    "4param":           ["4param"],
    "7param":           ["7param"],
    "hybrid":           ["hybrid"],
    "bondorder":        ["bondorder"],
    "bondorder_ext":    ["bondorder_ext"],
    "bondorder_ar":     ["bondorder_ar"],
    "extended":         ["extended"],
    "neighbour":        ["neighbour"],
    "both":             ["4param", "7param"],
    "all":              ["4param", "7param", "hybrid", "bondorder", "bondorder_ext",
                         "bondorder_ar", "extended", "neighbour"],
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


def _functional_group_frequencies(smiles_list):
    """Return functional group counts across `smiles_list`, sorted descending, excluding zero-count groups."""
    from rdkit import Chem
    from rdkit.Chem import Fragments

    nitramine    = Chem.MolFromSmarts("[NX3][N+](=O)[O-]")
    nitroso      = Chem.MolFromSmarts("[NX2]=[OX1]")
    tert_amine   = Chem.MolFromSmarts("[NX3;H0;!$(N=*);!$(N#*);!$([N+])]")
    groups = [
        ("Nitro (-NO2)",      lambda m: Fragments.fr_nitro(m) > 0),
        ("Nitramine (N-NO2)", lambda m: m.HasSubstructMatch(nitramine)),
        ("Nitroso (N=O)",     lambda m: m.HasSubstructMatch(nitroso)),
        ("Primary amine",     lambda m: Fragments.fr_NH2(m) > 0),
        ("Secondary amine",   lambda m: Fragments.fr_NH1(m) > 0),
        ("Tertiary amine",    lambda m: m.HasSubstructMatch(tert_amine)),
        ("Carbonyl (C=O)",    lambda m: Fragments.fr_C_O(m) > 0),
        ("Carboxylic acid",   lambda m: Fragments.fr_COO(m) > 0),
        ("Hydroxyl (-OH)",    lambda m: (Fragments.fr_Al_OH(m) + Fragments.fr_Ar_OH(m)) > 0),
        ("Nitrile (-CN)",     lambda m: Fragments.fr_nitrile(m) > 0),
        ("Aromatic ring",     lambda m: any(a.GetIsAromatic() for a in m.GetAtoms())),
        ("3-membered ring",   lambda m: any(len(r) == 3 for r in m.GetRingInfo().AtomRings())),
    ]
    mols = [m for m in (Chem.MolFromSmiles(s) for s in smiles_list) if m is not None]
    result = [(name, sum(1 for m in mols if detect(m))) for name, detect in groups]
    result = [(name, count) for name, count in result if count > 0]
    result.sort(key=lambda x: -x[1])
    return result, len(mols)


def _print_model_results(counts, u_values, exp_dhf, param_names, kfold, mol_names=None, n_outliers=0, mol_smiles=None):
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

    if n_outliers > 0 and mol_names is not None:
        errors_abs = [(abs(predicted[j] - exp_dhf[j]), j) for j in range(len(predicted))]
        errors_abs.sort(reverse=True)
        top = errors_abs[:n_outliers]
        print()
        print(f"   Top {len(top)} outliers by |error|:")
        print(f"   {'Rank':<6}{'Name':<40}{'Exp ΔHf°':>12}{'Pred ΔHf°':>12}{'Error':>10}")
        print(f"   {'─' * 82}")
        for rank, (_, j) in enumerate(top, 1):
            err = predicted[j] - exp_dhf[j]
            print(f"   {rank:<6}{mol_names[j]:<40}{exp_dhf[j]:>12.2f}{predicted[j]:>12.2f}{err:>10.2f}")

        if mol_smiles is not None:
            top_smiles = [mol_smiles[j] for _, j in top]
            fg_counts, n_fg_mols = _functional_group_frequencies(top_smiles)
            if fg_counts:
                print()
                print(f"   Functional groups in outliers ({n_fg_mols} molecules):")
                print(f"   {'Group':<25}{'Count':>6}{'%':>5}")
                print(f"   {'─' * 38}")
                for name, count in fg_counts:
                    print(f"   {name:<25}{count:>6}{100 * count // n_fg_mols:>4}%")

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
        "--model",
        choices=["4param", "7param", "hybrid", "bondorder", "bondorder_ext",
                 "bondorder_ar", "extended", "neighbour", "both", "all"],
        default="both",
    )
    fit_parser.add_argument("--kfold", type=int, default=10, help="Number of CV folds")
    fit_parser.add_argument("--n-conformers", type=int, default=1, help="Number of conformers to optimize")
    fit_parser.add_argument(
        "--optimizer", choices=["xtb", "uma", "esen", "aimnet2"], default="xtb",
        help=(
            "Geometry optimizer: 'xtb' (default, requires xTB binary), "
            "'uma'/'esen' (requires fairchem-core and MODEL_DIR), 'aimnet2' (requires aimnet2calc)"
        ),
    )
    fit_parser.add_argument("--output", "-o", help="Output JSON file for fitted epsilon values")
    fit_parser.add_argument("--csv", help="Output CSV with training data and predictions")
    fit_parser.add_argument(
        "--use-xtb-wbos", action="store_true",
        help="Use xTB Wiberg bond orders (instead of RDKit) for 7-param atom classification (requires --optimizer xtb)",
    )
    fit_parser.add_argument(
        "--use-gxtb", action="store_true",
        help="Use gxtb single-point energies (wB97M-V/def2-TZVPPD). WARNING: Do not mix with xTB workflows!",
    )
    fit_parser.add_argument(
        "--cache-dir", type=str, default=None,
        help="Directory for caching optimization results (enables restart capability)",
    )
    fit_parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print per-molecule details instead of a progress bar",
    )
    fit_parser.add_argument(
        "--outliers", nargs="?", const=10, default=0, type=int, metavar="N",
        help="Print the top N outliers by |error| after each model's statistics (default N=10 if flag given)",
    )

    # --- predict subcommand ---
    pred_parser = subparsers.add_parser("predict", help="Predict DeltaHf for new molecules")
    pred_parser.add_argument("--input", "-i", required=True, help="CSV with smiles column")
    pred_parser.add_argument(
        "--epsilon", required=False,
        help="JSON file with fitted epsilon values (uses defaults if not specified)",
    )
    pred_parser.add_argument(
        "--model",
        choices=["4param", "7param", "hybrid", "bondorder", "bondorder_ext",
                 "bondorder_ar", "extended", "neighbour"],
        default="4param",
    )
    pred_parser.add_argument("--n-conformers", type=int, default=1, help="Number of conformers to optimize")
    pred_parser.add_argument(
        "--optimizer", choices=["xtb", "uma", "esen", "aimnet2"], default="xtb",
        help="Geometry optimizer. Must match the optimizer used during fit!",
    )
    pred_parser.add_argument("--output", "-o", help="Output CSV with results")
    pred_parser.add_argument(
        "--use-xtb-wbos", action="store_true",
        help="Use xTB Wiberg bond orders (instead of RDKit) for 7-param atom classification (requires --optimizer xtb)",
    )
    pred_parser.add_argument(
        "--use-gxtb", action="store_true",
        help="Use gxtb single-point energies (wB97M-V/def2-TZVPPD). Must match method used in fit!",
    )
    pred_parser.add_argument(
        "--cache-dir", type=str, default=None,
        help="Directory for caching optimization results (enables restart capability)",
    )
    pred_parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print per-molecule details instead of a progress bar",
    )

    return parser


def cmd_fit(args):
    """Run the fitting workflow."""
    import pandas as pd

    if args.use_xtb_wbos and args.optimizer != "xtb":
        print(f"   Error: --use-xtb-wbos requires --optimizer xtb (got '{args.optimizer}')")
        sys.exit(1)

    predictor = None
    if args.optimizer == "xtb":
        from deltahf.xtb import find_xtb_binary
        try:
            find_xtb_binary()
        except FileNotFoundError as e:
            print(f"   Error: {e}")
            sys.exit(1)
    else:
        from deltahf.uma import load_mlip_model
        print(f"   Loading {args.optimizer} model...")
        try:
            predictor = load_mlip_model(args.optimizer)
        except (ImportError, FileNotFoundError) as e:
            print(f"   Error: {e}")
            sys.exit(1)
        print(f"   Using {args.optimizer} geometry optimizer")

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

        # Append method suffix to keep caches isolated per optimizer/energy method
        cache_path = Path(args.cache_dir)
        method_parts = [args.optimizer]
        if args.use_gxtb:
            method_parts.append("gxtb")
        method_suffix = "_" + "_".join(method_parts)
        if not cache_path.name.endswith(method_suffix):
            cache_path = cache_path.parent / (cache_path.name + method_suffix)

        cache = ResultCache(cache_path)
        print(f"   Using cache directory: {cache_path}")
        if args.use_xtb_wbos:
            print("      Note: caching disabled when --use-xtb-wbos is set")
            cache = None

    # Process all molecules
    print(f"   Running {args.optimizer} optimizations...")
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
                optimizer=args.optimizer, predictor=predictor,
                use_xtb_wbos=args.use_xtb_wbos, use_gxtb=args.use_gxtb, cache=cache,
            )
            if result.error:
                print(f"         ERROR: {result.error}")
            else:
                opt_energy = result.xtb_energy if args.optimizer == "xtb" else result.mlip_energy
                print(f"         {args.optimizer} energy: {opt_energy:.6f} Eh")
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
                optimizer=args.optimizer, predictor=predictor,
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
        print(f"\n   Is {args.optimizer} available? Check your installation.")
        sys.exit(1)

    indices = [i for i, _ in successful]
    u_values = [_best_energy_kcal(r) for _, r in successful]
    exp_dhf = [df.iloc[i]["exp_dhf_kcal_mol"] for i in indices]

    output = {}
    predicted = {}
    models_to_run = MODEL_GROUPS[args.model]
    mol_names = [df.iloc[i].get("name", df.iloc[i]["smiles"]) for i in indices]
    mol_smiles = [df.iloc[i]["smiles"] for i in indices]

    for model_name in models_to_run:
        param_names, counts_attr = MODEL_DEFS[model_name]
        counts = [getattr(results[i], counts_attr) for i in indices]
        eps, pred = _print_model_results(
            counts, u_values, exp_dhf, param_names, args.kfold,
            mol_names=mol_names, n_outliers=args.outliers, mol_smiles=mol_smiles,
        )
        output[model_name] = eps
        predicted[model_name] = pred

    # Add metadata to track optimizer and energy method used for fitting
    output["_metadata"] = {
        "optimizer": args.optimizer,
        "method": "gxtb" if args.use_gxtb else args.optimizer,
        "n_conformers": args.n_conformers,
        "n_molecules": len(indices),
    }

    print()
    print(SEP)

    if args.output:
        with open(args.output, "w") as f:
            json.dump(output, f, indent=2)
        print(f"   Fitted parameters saved to {args.output}")
        energy_method = "gxtb" if args.use_gxtb else args.optimizer
        print(f"   Optimizer: {args.optimizer}  |  Energy method: {energy_method}")

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
    if args.use_xtb_wbos and args.optimizer != "xtb":
        print(f"   Error: --use-xtb-wbos requires --optimizer xtb (got '{args.optimizer}')")
        sys.exit(1)

    predictor = None
    if args.optimizer == "xtb":
        from deltahf.xtb import find_xtb_binary
        try:
            find_xtb_binary()
        except FileNotFoundError as e:
            print(f"   Error: {e}")
            sys.exit(1)
    else:
        from deltahf.uma import load_mlip_model
        print(f"   Loading {args.optimizer} model...")
        try:
            predictor = load_mlip_model(args.optimizer)
        except (ImportError, FileNotFoundError) as e:
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
        elif args.optimizer != "xtb":
            epsilon_file = params_dir / f"params_{args.optimizer}.json"
        else:
            epsilon_file = params_dir / "params_xtb.json"

        if not epsilon_file.exists():
            print(f"   Error: Default parameter file not found: {epsilon_file}")
            sys.exit(1)

        print(f"   WARNING: Using default parameters from {epsilon_file}")
        print("            Specify --epsilon to use custom parameters\n")

    with open(epsilon_file) as f:
        epsilon_data = json.load(f)

    # Check for method/optimizer mismatch
    metadata = epsilon_data.get("_metadata", {})
    param_method = metadata.get("method", "xtb")   # energy method used during fit
    param_optimizer = metadata.get("optimizer", metadata.get("method", "xtb"))  # optimizer used during fit
    current_method = "gxtb" if args.use_gxtb else args.optimizer
    current_optimizer = args.optimizer

    if param_method != current_method or param_optimizer != current_optimizer:
        print("\n   WARNING: Method mismatch detected!")
        print(f"   Parameters fitted with optimizer={param_optimizer}, energy={param_method}")
        print(f"   You are predicting with optimizer={current_optimizer}, energy={current_method}")
        print("   This will produce incorrect results - energies are on different scales!")
        print("\n   Either:")
        print(f"     - Match the fit settings (--optimizer {param_optimizer}" +
              (" --use-gxtb" if param_method == "gxtb" else "") + ")")
        print("     - Refit parameters using your current settings\n")
        sys.exit(1)

    epsilon = epsilon_data.get(args.model, epsilon_data)

    cache = None
    if args.cache_dir:
        from deltahf.cache import ResultCache

        # Append method suffix to keep caches isolated per optimizer/energy method
        cache_path = Path(args.cache_dir)
        method_parts = [args.optimizer]
        if args.use_gxtb:
            method_parts.append("gxtb")
        method_suffix = "_" + "_".join(method_parts)
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
        optimizer=args.optimizer,
        predictor=predictor,
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
    if args.optimizer == "xtb":
        display_cols.append("xtb_energy_kcal")
    else:
        display_cols.append("mlip_energy_kcal")
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
