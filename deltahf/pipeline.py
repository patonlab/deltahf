"""End-to-end pipeline orchestrating SMILES -> conformers -> xTB -> DHf prediction."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from deltahf.atom_equivalents import HARTREE_TO_KCAL, predict_dhf
from deltahf.conformers import (
    check_connectivity,
    generate_conformers,
    get_lowest_conformers,
    prune_conformers,
    write_xyz,
)
from deltahf.smiles import (
    classify_atoms_7param,
    classify_atoms_7param_from_wbo,
    classify_atoms_extended,
    classify_atoms_hybrid,
    count_atoms,
)
from deltahf.xtb import find_gxtb_binary, parse_wbo_file, run_gxtb_single_point, run_xtb_optimization

if TYPE_CHECKING:
    from deltahf.cache import ResultCache


@dataclass
class MoleculeResult:
    smiles: str
    name: str | None = None
    xtb_energy: float | None = None
    xtb_energy_kcal: float | None = None
    gxtb_energy: float | None = None
    gxtb_energy_kcal: float | None = None
    atom_counts_4param: dict | None = None
    atom_counts_7param: dict | None = None
    atom_counts_hybrid: dict | None = None
    atom_counts_extended: dict | None = None
    dhf_4param: float | None = None
    dhf_7param: float | None = None
    dhf_hybrid: float | None = None
    dhf_extended: float | None = None
    n_conformers_generated: int = 0
    n_conformers_optimized: int = 0
    n_conformers_isomerized: int = 0
    error: str | None = None


def process_molecule(
    smiles: str,
    n_conformers: int = 5,
    num_initial_confs: int = 50,
    epsilon_4param: dict | None = None,
    epsilon_7param: dict | None = None,
    epsilon_hybrid: dict | None = None,
    epsilon_extended: dict | None = None,
    work_dir: Path | None = None,
    name: str | None = None,
    use_xtb_wbos: bool = False,
    use_gxtb: bool = False,
    cache: ResultCache | None = None,
) -> MoleculeResult:
    """Full pipeline for a single molecule."""
    result = MoleculeResult(smiles=smiles, name=name)

    try:
        result.atom_counts_4param = count_atoms(smiles)
        if not use_xtb_wbos:
            result.atom_counts_7param = classify_atoms_7param(smiles)
        result.atom_counts_hybrid = classify_atoms_hybrid(smiles)
        result.atom_counts_extended = classify_atoms_extended(smiles)

        # Check cache before running xTB
        if cache is not None and not use_xtb_wbos:
            cached = cache.lookup(smiles, n_conformers)
            if cached is not None:
                result.xtb_energy = cached.xtb_energy
                result.xtb_energy_kcal = cached.xtb_energy_kcal
                result.gxtb_energy = cached.gxtb_energy
                result.gxtb_energy_kcal = cached.gxtb_energy_kcal
                result.n_conformers_optimized = cached.n_conformers_optimized
                result.n_conformers_isomerized = cached.n_conformers_isomerized
                # Use gxtb energy for predictions if available, otherwise fall back to xtb
                energy_for_prediction = result.gxtb_energy_kcal if result.gxtb_energy_kcal is not None else result.xtb_energy_kcal
                if epsilon_4param:
                    result.dhf_4param = predict_dhf(
                        energy_for_prediction, result.atom_counts_4param, epsilon_4param
                    )
                if epsilon_7param:
                    result.dhf_7param = predict_dhf(
                        energy_for_prediction, result.atom_counts_7param, epsilon_7param
                    )
                if epsilon_hybrid:
                    result.dhf_hybrid = predict_dhf(
                        energy_for_prediction, result.atom_counts_hybrid, epsilon_hybrid
                    )
                if epsilon_extended:
                    result.dhf_extended = predict_dhf(
                        energy_for_prediction, result.atom_counts_extended, epsilon_extended
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

            xtb_result = run_xtb_optimization(xyz_path)
            if xtb_result.converged:
                result.n_conformers_optimized += 1
                if xtb_result.optimized_xyz_path and not check_connectivity(mol, xtb_result.optimized_xyz_path):
                    result.n_conformers_isomerized += 1
                    continue
                if xtb_result.energy < best_energy:
                    best_energy = xtb_result.energy
                    best_wbo_path = xtb_result.wbo_path
                    best_xyz_path = xtb_result.optimized_xyz_path

        if best_energy == float("inf"):
            if result.n_conformers_isomerized > 0:
                result.error = (
                    f"All {result.n_conformers_isomerized} converged conformer(s) "
                    "isomerized during optimization"
                )
            else:
                result.error = "All xTB optimizations failed"
            return result

        result.xtb_energy = best_energy
        result.xtb_energy_kcal = best_energy * HARTREE_TO_KCAL

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
                    xtb_energy=result.xtb_energy,
                    xtb_energy_kcal=result.xtb_energy_kcal,
                    n_conformers=n_conformers,
                    n_conformers_optimized=result.n_conformers_optimized,
                    n_conformers_isomerized=result.n_conformers_isomerized,
                    gfn_level=2,
                    charge=0,
                    gxtb_energy=result.gxtb_energy,
                    gxtb_energy_kcal=result.gxtb_energy_kcal,
                )
            )
            cache.save()

        if use_xtb_wbos and best_wbo_path is not None:
            wbos = parse_wbo_file(best_wbo_path)
            result.atom_counts_7param = classify_atoms_7param_from_wbo(smiles, wbos)

        # Use gxtb energy for predictions if available, otherwise fall back to xtb
        energy_for_prediction = result.gxtb_energy_kcal if result.gxtb_energy_kcal is not None else result.xtb_energy_kcal

        if epsilon_4param:
            result.dhf_4param = predict_dhf(energy_for_prediction, result.atom_counts_4param, epsilon_4param)
        if epsilon_7param:
            result.dhf_7param = predict_dhf(energy_for_prediction, result.atom_counts_7param, epsilon_7param)
        if epsilon_hybrid:
            result.dhf_hybrid = predict_dhf(energy_for_prediction, result.atom_counts_hybrid, epsilon_hybrid)
        if epsilon_extended:
            result.dhf_extended = predict_dhf(energy_for_prediction, result.atom_counts_extended, epsilon_extended)

    except Exception as e:
        result.error = str(e)

    return result


def process_csv(
    csv_path: Path,
    n_conformers: int = 5,
    epsilon_4param: dict | None = None,
    epsilon_7param: dict | None = None,
    epsilon_hybrid: dict | None = None,
    epsilon_extended: dict | None = None,
    output_path: Path | None = None,
    use_xtb_wbos: bool = False,
    use_gxtb: bool = False,
    cache: ResultCache | None = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Process a CSV of SMILES and return results."""
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
            epsilon_4param=epsilon_4param,
            epsilon_7param=epsilon_7param,
            epsilon_hybrid=epsilon_hybrid,
            epsilon_extended=epsilon_extended,
            name=mol_name,
            use_xtb_wbos=use_xtb_wbos,
            use_gxtb=use_gxtb,
            cache=cache,
        )
        results.append(mol_result)

    results_df = pd.DataFrame([vars(r) for r in results])
    if output_path:
        results_df.to_csv(output_path, index=False)
    return results_df
