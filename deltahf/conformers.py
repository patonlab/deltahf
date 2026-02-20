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


def prune_conformers(
    mol: Chem.Mol,
    energies: list[tuple[float, int]],
    rmsd_threshold: float = 0.125,
) -> list[tuple[float, int]]:
    """Remove near-duplicate conformers by RMSD.

    Greedy filter over the energy-sorted list: keep a conformer only if its
    RMSD to all previously kept conformers exceeds the threshold.
    Always keeps at least the lowest-energy conformer.
    """
    if len(energies) <= 1:
        return list(energies)

    kept: list[tuple[float, int]] = [energies[0]]
    for energy, cid in energies[1:]:
        is_unique = True
        for _, kept_cid in kept:
            rms = AllChem.GetConformerRMS(mol, cid, kept_cid)
            if rms < rmsd_threshold:
                is_unique = False
                break
        if is_unique:
            kept.append((energy, cid))
    return kept


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


def _heavy_adj_equal(mol_a: Chem.Mol, mol_b: Chem.Mol) -> bool:
    """Compare adjacency matrices for heavy atoms only.

    Ignores differences in H connectivity (tautomers, proton transfers).
    Falls back to full-adjacency comparison for all-hydrogen molecules
    (e.g. H2) where there are no heavy atoms to compare.
    """
    heavy_a = [a.GetIdx() for a in mol_a.GetAtoms() if a.GetAtomicNum() != 1]
    heavy_b = [a.GetIdx() for a in mol_b.GetAtoms() if a.GetAtomicNum() != 1]
    if len(heavy_a) != len(heavy_b):
        return False
    adj_a = rdmolops.GetAdjacencyMatrix(mol_a, useBO=False)
    adj_b = rdmolops.GetAdjacencyMatrix(mol_b, useBO=False)
    # No heavy atoms (e.g. H2): fall back to full adjacency
    if not heavy_a:
        return bool(np.array_equal(adj_a, adj_b))
    return bool(np.array_equal(
        adj_a[np.ix_(heavy_a, heavy_a)],
        adj_b[np.ix_(heavy_b, heavy_b)],
    ))


def check_connectivity(mol_ref: Chem.Mol, xyz_path: Path) -> bool:
    """Check whether optimized geometry has the same heavy-atom connectivity as the reference.

    Reads the optimized XYZ file, perceives bonds using RDKit's
    DetermineConnectivity, and compares heavy-atom adjacency matrices with
    the reference molecule. H-only changes (tautomerism, proton transfers)
    are tolerated; changes to heavy-atom bonding are rejected.

    Uses the connect-the-dots method first. If that finds a mismatch
    (e.g. slightly stretched bonds in H2), retries with the van der Waals
    method at a more generous covFactor before declaring isomerization.

    Returns True if heavy-atom connectivity matches, False if isomerized or unreadable.
    """
    mol_opt = Chem.MolFromXYZFile(str(xyz_path))
    if mol_opt is None:
        return False

    if mol_ref.GetNumAtoms() != mol_opt.GetNumAtoms():
        return False

    # Primary: connect-the-dots (default)
    rdDetermineBonds.DetermineConnectivity(mol_opt)
    if _heavy_adj_equal(mol_ref, mol_opt):
        return True

    # Fallback: VdW method with generous covFactor for stretched bonds
    mol_opt2 = Chem.MolFromXYZFile(str(xyz_path))
    rdDetermineBonds.DetermineConnectivity(mol_opt2, useVdw=True, covFactor=2.0)
    return _heavy_adj_equal(mol_ref, mol_opt2)
