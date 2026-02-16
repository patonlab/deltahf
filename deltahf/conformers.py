"""3D conformer generation via RDKit ETKDG."""

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

from deltahf.smiles import smiles_to_mol


def generate_conformers(
    smiles: str,
    num_confs: int = 50,
    random_seed: int = 42,
    num_threads: int = 0,
) -> tuple[Chem.Mol, list[tuple[float, int]]]:
    """Generate 3D conformers via ETKDG and rank by MMFF energy.

    Returns (mol_with_conformers, sorted_list_of_(energy, conf_id)).
    """
    mol = smiles_to_mol(smiles)

    params = rdDistGeom.ETKDGv3()
    params.randomSeed = random_seed
    params.numThreads = num_threads

    cids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    if len(cids) == 0:
        raise RuntimeError(f"Failed to generate conformers for {smiles}")

    results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=num_threads)
    energies = []
    for cid, (converged, energy) in zip(cids, results):
        if converged == 0:  # 0 means converged in RDKit
            energies.append((energy, cid))

    if not energies:
        # Fall back to using unconverged results with their energies
        for cid, (_, energy) in zip(cids, results):
            energies.append((energy, cid))

    energies.sort()
    return mol, energies


def get_lowest_conformers(
    mol: Chem.Mol,
    energies: list[tuple[float, int]],
    n: int,
) -> list[int]:
    """Return conformer IDs for the n lowest-energy conformers."""
    return [cid for _, cid in energies[:n]]


def write_xyz(mol: Chem.Mol, conf_id: int, path: Path) -> None:
    """Write a single conformer to XYZ format."""
    conf = mol.GetConformer(conf_id)
    natoms = mol.GetNumAtoms()
    lines = [str(natoms), ""]
    for i in range(natoms):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        lines.append(f"{atom.GetSymbol():2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")
    path.write_text("\n".join(lines) + "\n")
