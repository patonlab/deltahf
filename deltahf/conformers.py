"""3D conformer generation via RDKit ETKDG."""

from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds, rdDistGeom, rdmolops

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


def check_connectivity(mol_ref: Chem.Mol, xyz_path: Path) -> bool:
    """Check whether optimized geometry has the same connectivity as the reference.

    Reads the optimized XYZ file, perceives bonds using RDKit's
    DetermineConnectivity, and compares adjacency matrices (ignoring
    bond orders) with the reference molecule.

    Uses the connect-the-dots method first. If that finds a mismatch
    (e.g. slightly stretched bonds in H2), retries with the van der Waals
    method at a more generous covFactor before declaring isomerization.

    Returns True if connectivity matches, False if isomerized or unreadable.
    """
    mol_opt = Chem.MolFromXYZFile(str(xyz_path))
    if mol_opt is None:
        return False

    adj_ref = rdmolops.GetAdjacencyMatrix(mol_ref, useBO=False)

    if mol_ref.GetNumAtoms() != mol_opt.GetNumAtoms():
        return False

    # Primary: connect-the-dots (default)
    rdDetermineBonds.DetermineConnectivity(mol_opt)
    adj_opt = rdmolops.GetAdjacencyMatrix(mol_opt, useBO=False)

    if np.array_equal(adj_ref, adj_opt):
        return True

    # Fallback: VdW method with generous covFactor for stretched bonds
    mol_opt2 = Chem.MolFromXYZFile(str(xyz_path))
    rdDetermineBonds.DetermineConnectivity(mol_opt2, useVdw=True, covFactor=2.0)
    adj_opt2 = rdmolops.GetAdjacencyMatrix(mol_opt2, useBO=False)

    return bool(np.array_equal(adj_ref, adj_opt2))
