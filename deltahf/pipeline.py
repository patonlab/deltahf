"""End-to-end pipeline orchestrating SMILES -> conformers -> xTB -> DHf prediction."""

import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from deltahf.atom_equivalents import HARTREE_TO_KCAL, predict_dhf
from deltahf.conformers import generate_conformers, get_lowest_conformers, write_xyz
from deltahf.smiles import classify_atoms_7param, count_atoms
from deltahf.xtb import run_xtb_optimization


@dataclass
class MoleculeResult:
    smiles: str
    name: str | None = None
    xtb_energy: float | None = None
    xtb_energy_kcal: float | None = None
    atom_counts_4param: dict | None = None
    atom_counts_7param: dict | None = None
    dhf_4param: float | None = None
    dhf_7param: float | None = None
    n_conformers_generated: int = 0
    n_conformers_optimized: int = 0
    error: str | None = None


def process_molecule(
    smiles: str,
    n_conformers: int = 5,
    num_initial_confs: int = 50,
    epsilon_4param: dict | None = None,
    epsilon_7param: dict | None = None,
    work_dir: Path | None = None,
    name: str | None = None,
) -> MoleculeResult:
    """Full pipeline for a single molecule."""
    result = MoleculeResult(smiles=smiles, name=name)

    try:
        result.atom_counts_4param = count_atoms(smiles)
        result.atom_counts_7param = classify_atoms_7param(smiles)

        mol, energies = generate_conformers(smiles, num_confs=num_initial_confs)
        result.n_conformers_generated = len(energies)

        lowest_cids = get_lowest_conformers(mol, energies, n=n_conformers)

        best_energy = float("inf")
        if work_dir is None:
            work_dir = Path(tempfile.mkdtemp())

        for i, cid in enumerate(lowest_cids):
            conf_dir = work_dir / f"conf_{i}"
            conf_dir.mkdir(parents=True, exist_ok=True)
            xyz_path = conf_dir / "input.xyz"
            write_xyz(mol, cid, xyz_path)

            xtb_result = run_xtb_optimization(xyz_path)
            if xtb_result.converged and xtb_result.energy < best_energy:
                best_energy = xtb_result.energy
                result.n_conformers_optimized += 1

        if best_energy == float("inf"):
            result.error = "All xTB optimizations failed"
            return result

        result.xtb_energy = best_energy
        result.xtb_energy_kcal = best_energy * HARTREE_TO_KCAL

        if epsilon_4param:
            result.dhf_4param = predict_dhf(result.xtb_energy_kcal, result.atom_counts_4param, epsilon_4param)
        if epsilon_7param:
            result.dhf_7param = predict_dhf(result.xtb_energy_kcal, result.atom_counts_7param, epsilon_7param)

    except Exception as e:
        result.error = str(e)

    return result


def process_csv(
    csv_path: Path,
    n_conformers: int = 5,
    epsilon_4param: dict | None = None,
    epsilon_7param: dict | None = None,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Process a CSV of SMILES and return results."""
    df = pd.read_csv(csv_path)
    results = []

    for _, row in df.iterrows():
        smiles = row["smiles"]
        mol_name = row.get("name", None)
        mol_result = process_molecule(
            smiles,
            n_conformers=n_conformers,
            epsilon_4param=epsilon_4param,
            epsilon_7param=epsilon_7param,
            name=mol_name,
        )
        results.append(mol_result)

    results_df = pd.DataFrame([vars(r) for r in results])
    if output_path:
        results_df.to_csv(output_path, index=False)
    return results_df
