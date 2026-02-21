"""End-to-end pipeline orchestrating SMILES -> conformers -> optimizer -> DHf prediction."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from deltahf.atom_equivalents import predict_dhf
from deltahf.conformers import (
    check_connectivity,
    generate_conformers,
    get_lowest_conformers,
    prune_conformers,
    write_xyz,
)
from deltahf.constants import HARTREE_TO_KCAL
from deltahf.smiles import (
    classify_atoms_7param,
    classify_atoms_7param_from_wbo,
    classify_atoms_bondorder,
    classify_atoms_bondorder_ar,
    classify_atoms_bondorder_ext,
    classify_atoms_extended,
    classify_atoms_hybrid,
    classify_atoms_neighbour,
    count_atoms,
)
from deltahf.xtb import find_gxtb_binary, parse_wbo_file, run_gxtb_single_point, run_xtb_optimization

if TYPE_CHECKING:
    from deltahf.cache import ResultCache


@dataclass
class MoleculeResult:
    """Result of processing a single molecule through the pipeline.

    Energy fields are in Hartree (``_energy``) and kcal/mol (``_energy_kcal``).
    For prediction, energy priority is: gxtb > mlip > xtb.
    ``mlip_energy`` is populated for any non-xTB optimizer (UMA, eSEN, AIMNet2).
    """

    smiles: str
    name: str | None = None
    xtb_energy: float | None = None
    xtb_energy_kcal: float | None = None
    mlip_energy: float | None = None       # Populated when optimizer != "xtb"
    mlip_energy_kcal: float | None = None
    gxtb_energy: float | None = None
    gxtb_energy_kcal: float | None = None
    atom_counts_element: dict | None = None
    atom_counts_element_bo: dict | None = None
    atom_counts_hybrid: dict | None = None
    atom_counts_bondorder: dict | None = None
    atom_counts_bondorder_ext: dict | None = None
    atom_counts_bondorder_ar: dict | None = None
    atom_counts_extended: dict | None = None
    atom_counts_neighbour: dict | None = None
    dhf_element: float | None = None
    dhf_element_bo: float | None = None
    dhf_hybrid: float | None = None
    dhf_bondorder: float | None = None
    dhf_bondorder_ext: float | None = None
    dhf_bondorder_ar: float | None = None
    dhf_extended: float | None = None
    dhf_neighbour: float | None = None
    n_conformers_generated: int = 0
    n_conformers_optimized: int = 0
    n_conformers_isomerized: int = 0
    error: str | None = None


def _best_energy_kcal(result: "MoleculeResult") -> float | None:
    """Return the best available energy (kcal/mol) using priority: gxtb > mlip > xtb."""
    return result.gxtb_energy_kcal or result.mlip_energy_kcal or result.xtb_energy_kcal


def _apply_predictions(
    result: "MoleculeResult",
    energy_kcal: float | None,
    epsilon_element: dict | None,
    epsilon_element_bo: dict | None,
    epsilon_hybrid: dict | None,
    epsilon_bondorder: dict | None,
    epsilon_bondorder_ext: dict | None,
    epsilon_bondorder_ar: dict | None,
    epsilon_extended: dict | None,
    epsilon_neighbour: dict | None,
) -> None:
    """Populate all dhf_* fields on result from the given energy and epsilon dicts."""
    if epsilon_element:
        result.dhf_element = predict_dhf(energy_kcal, result.atom_counts_element, epsilon_element)
    if epsilon_element_bo:
        result.dhf_element_bo = predict_dhf(energy_kcal, result.atom_counts_element_bo, epsilon_element_bo)
    if epsilon_hybrid:
        result.dhf_hybrid = predict_dhf(energy_kcal, result.atom_counts_hybrid, epsilon_hybrid)
    if epsilon_bondorder:
        result.dhf_bondorder = predict_dhf(energy_kcal, result.atom_counts_bondorder, epsilon_bondorder)
    if epsilon_bondorder_ext:
        result.dhf_bondorder_ext = predict_dhf(energy_kcal, result.atom_counts_bondorder_ext, epsilon_bondorder_ext)
    if epsilon_bondorder_ar:
        result.dhf_bondorder_ar = predict_dhf(energy_kcal, result.atom_counts_bondorder_ar, epsilon_bondorder_ar)
    if epsilon_extended:
        result.dhf_extended = predict_dhf(energy_kcal, result.atom_counts_extended, epsilon_extended)
    if epsilon_neighbour:
        result.dhf_neighbour = predict_dhf(energy_kcal, result.atom_counts_neighbour, epsilon_neighbour)


def process_molecule(
    smiles: str,
    n_conformers: int = 1,
    num_initial_confs: int = 50,
    epsilon_element: dict | None = None,
    epsilon_element_bo: dict | None = None,
    epsilon_hybrid: dict | None = None,
    epsilon_bondorder: dict | None = None,
    epsilon_bondorder_ext: dict | None = None,
    epsilon_bondorder_ar: dict | None = None,
    epsilon_extended: dict | None = None,
    epsilon_neighbour: dict | None = None,
    work_dir: Path | None = None,
    name: str | None = None,
    optimizer: str = "xtb",
    predictor=None,
    use_xtb_wbos: bool = False,
    use_gxtb: bool = False,
    cache: ResultCache | None = None,
    xtb_threads: int | None = None,
) -> MoleculeResult:
    """Full pipeline for a single molecule."""
    result = MoleculeResult(smiles=smiles, name=name)

    try:
        # Validate molecule contains only supported elements (CHNO)
        from deltahf.smiles import validate_elements
        validate_elements(smiles)

        result.atom_counts_element = count_atoms(smiles)
        if not use_xtb_wbos:
            result.atom_counts_element_bo = classify_atoms_7param(smiles)
        result.atom_counts_hybrid = classify_atoms_hybrid(smiles)
        result.atom_counts_bondorder = classify_atoms_bondorder(smiles)
        result.atom_counts_bondorder_ext = classify_atoms_bondorder_ext(smiles)
        result.atom_counts_bondorder_ar = classify_atoms_bondorder_ar(smiles)
        result.atom_counts_extended = classify_atoms_extended(smiles)
        result.atom_counts_neighbour = classify_atoms_neighbour(smiles)
        from rdkit import Chem as _Chem
        charge = _Chem.GetFormalCharge(_Chem.MolFromSmiles(smiles))

        # Check cache before running optimizer
        if cache is not None and not use_xtb_wbos:
            cached = cache.lookup(smiles, n_conformers, optimizer=optimizer, charge=charge)
            if cached is not None:
                result.xtb_energy = cached.xtb_energy
                result.xtb_energy_kcal = cached.xtb_energy_kcal
                result.mlip_energy = cached.mlip_energy
                result.mlip_energy_kcal = cached.mlip_energy_kcal
                result.gxtb_energy = cached.gxtb_energy
                result.gxtb_energy_kcal = cached.gxtb_energy_kcal
                result.n_conformers_optimized = cached.n_conformers_optimized
                result.n_conformers_isomerized = cached.n_conformers_isomerized
                _apply_predictions(
                    result, _best_energy_kcal(result),
                    epsilon_element, epsilon_element_bo, epsilon_hybrid,
                    epsilon_bondorder, epsilon_bondorder_ext, epsilon_bondorder_ar,
                    epsilon_extended, epsilon_neighbour,
                )
                return result

        mol, energies = generate_conformers(smiles, num_confs=num_initial_confs)
        energies = prune_conformers(mol, energies)
        result.n_conformers_generated = len(energies)

        lowest_cids = get_lowest_conformers(mol, energies, n=n_conformers)

        best_energy = float("inf")
        best_wbo_path = None
        best_xyz_path = None
        if work_dir is None:
            work_dir = Path(tempfile.mkdtemp())

        for i, cid in enumerate(lowest_cids):
            conf_dir = work_dir / f"conf_{i}"
            conf_dir.mkdir(parents=True, exist_ok=True)
            xyz_path = conf_dir / "input.xyz"
            write_xyz(mol, cid, xyz_path)

            if optimizer == "xtb":
                opt_result = run_xtb_optimization(xyz_path, parallel=xtb_threads)
                opt_wbo_path = opt_result.wbo_path
            else:
                from deltahf.uma import run_mlip_optimization
                opt_result = run_mlip_optimization(xyz_path, model=optimizer, predictor=predictor, charge=charge)
                opt_wbo_path = None

            if opt_result.converged:
                result.n_conformers_optimized += 1
                if opt_result.optimized_xyz_path and not check_connectivity(mol, opt_result.optimized_xyz_path):
                    result.n_conformers_isomerized += 1
                    continue
                if opt_result.energy < best_energy:
                    best_energy = opt_result.energy
                    best_wbo_path = opt_wbo_path
                    best_xyz_path = opt_result.optimized_xyz_path

        if best_energy == float("inf"):
            if result.n_conformers_isomerized > 0:
                result.error = (
                    f"All {result.n_conformers_isomerized} converged conformer(s) "
                    "isomerized during optimization"
                )
            else:
                result.error = f"All {optimizer} optimizations failed"
            return result

        if optimizer == "xtb":
            result.xtb_energy = best_energy
            result.xtb_energy_kcal = best_energy * HARTREE_TO_KCAL
        else:
            result.mlip_energy = best_energy
            result.mlip_energy_kcal = best_energy * HARTREE_TO_KCAL

        # Optionally run gxtb single-point energy calculation on the best conformer
        if use_gxtb and best_xyz_path is not None:
            if find_gxtb_binary() is not None:
                gxtb_result = run_gxtb_single_point(best_xyz_path)
                if gxtb_result.converged:
                    result.gxtb_energy = gxtb_result.energy
                    result.gxtb_energy_kcal = gxtb_result.energy * HARTREE_TO_KCAL
                # If gxtb fails, we silently continue with xtb energy (no error)

        # Store in cache
        if cache is not None and not use_xtb_wbos:
            from rdkit import Chem as _Chem

            from deltahf.cache import CachedResult

            cache.store(
                CachedResult(
                    canonical_smiles=_Chem.MolToSmiles(_Chem.MolFromSmiles(smiles)),
                    optimizer=optimizer,
                    xtb_energy=result.xtb_energy,
                    xtb_energy_kcal=result.xtb_energy_kcal,
                    n_conformers=n_conformers,
                    n_conformers_optimized=result.n_conformers_optimized,
                    n_conformers_isomerized=result.n_conformers_isomerized,
                    gfn_level=2,
                    charge=charge,
                    gxtb_energy=result.gxtb_energy,
                    gxtb_energy_kcal=result.gxtb_energy_kcal,
                    mlip_energy=result.mlip_energy,
                    mlip_energy_kcal=result.mlip_energy_kcal,
                )
            )
            cache.save()

        if use_xtb_wbos and best_wbo_path is not None:
            wbos = parse_wbo_file(best_wbo_path)
            result.atom_counts_element_bo = classify_atoms_7param_from_wbo(smiles, wbos)

        _apply_predictions(
            result, _best_energy_kcal(result),
            epsilon_element, epsilon_element_bo, epsilon_hybrid,
            epsilon_bondorder, epsilon_bondorder_ext, epsilon_bondorder_ar,
            epsilon_extended, epsilon_neighbour,
        )

    except Exception as e:
        result.error = str(e)

    return result


def process_csv(
    csv_path: Path,
    n_conformers: int = 1,
    epsilon_element: dict | None = None,
    epsilon_element_bo: dict | None = None,
    epsilon_hybrid: dict | None = None,
    epsilon_bondorder: dict | None = None,
    epsilon_bondorder_ext: dict | None = None,
    epsilon_bondorder_ar: dict | None = None,
    epsilon_extended: dict | None = None,
    epsilon_neighbour: dict | None = None,
    output_path: Path | None = None,
    optimizer: str = "xtb",
    predictor=None,
    use_xtb_wbos: bool = False,
    use_gxtb: bool = False,
    cache: ResultCache | None = None,
    verbose: bool = False,
    xtb_threads: int | None = None,
) -> pd.DataFrame:
    """Process a CSV of SMILES and return a DataFrame of MoleculeResult fields.

    Reads ``smiles`` (required) and ``name`` (optional) columns from ``csv_path``.
    Calls ``process_molecule()`` for each row, forwarding all optimizer/cache args.
    If ``output_path`` is given, writes the results DataFrame to CSV before returning.
    The ``predictor`` object (pre-loaded MLIP model) should be passed when
    ``optimizer`` is not ``"xtb"``; it is expensive to load per-molecule.
    """
    df = pd.read_csv(csv_path)
    results = []

    if verbose:
        rows = df.iterrows()
    else:
        from tqdm import tqdm

        rows = tqdm(df.iterrows(), total=len(df), desc="  Molecules")

    for _, row in rows:
        smiles = row["smiles"]
        mol_name = row.get("name", None)
        mol_result = process_molecule(
            smiles,
            n_conformers=n_conformers,
            epsilon_element=epsilon_element,
            epsilon_element_bo=epsilon_element_bo,
            epsilon_hybrid=epsilon_hybrid,
            epsilon_bondorder=epsilon_bondorder,
            epsilon_bondorder_ext=epsilon_bondorder_ext,
            epsilon_bondorder_ar=epsilon_bondorder_ar,
            epsilon_extended=epsilon_extended,
            epsilon_neighbour=epsilon_neighbour,
            name=mol_name,
            optimizer=optimizer,
            predictor=predictor,
            use_xtb_wbos=use_xtb_wbos,
            use_gxtb=use_gxtb,
            cache=cache,
            xtb_threads=xtb_threads,
        )
        results.append(mol_result)

    results_df = pd.DataFrame([vars(r) for r in results])
    if output_path:
        results_df.to_csv(output_path, index=False)
    return results_df
