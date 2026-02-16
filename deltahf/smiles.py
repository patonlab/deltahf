"""SMILES parsing, atom counting, and bond order classification."""

from rdkit import Chem


def smiles_to_mol(smiles: str) -> Chem.Mol:
    """Convert SMILES to RDKit mol object with explicit hydrogens.

    Raises ValueError if SMILES cannot be parsed.
    """
    if not smiles or not smiles.strip():
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    return Chem.AddHs(mol)


def count_atoms(smiles: str) -> dict[str, int]:
    """Count C, H, N, O atoms for the 4-parameter model."""
    mol = smiles_to_mol(smiles)
    counts = {"C": 0, "H": 0, "N": 0, "O": 0}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in counts:
            counts[symbol] += 1
    return counts


def classify_atoms_7param(smiles: str) -> dict[str, int]:
    """Classify atoms into 7 types for the 7-parameter model.

    C, H, N, O for single-bonded atoms; C', N', O' for atoms with any
    bond order > 1.25 (double or triple bonds).

    The molecule is kekulized first so that aromatic bonds are resolved
    into explicit single/double bonds. This ensures that atoms like O in
    furazan rings or [nH] in heterocycles are classified based on their
    Kekule bond orders rather than the uniform aromatic bond order of 1.5.

    Uses RDKit GetBondTypeAsDouble(): SINGLE=1.0, DOUBLE=2.0, TRIPLE=3.0.
    """
    mol = smiles_to_mol(smiles)
    Chem.Kekulize(mol, clearAromaticFlags=False)
    counts = {"C": 0, "H": 0, "N": 0, "O": 0, "C_prime": 0, "N_prime": 0, "O_prime": 0}

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        if _is_multiply_bonded(atom):
            counts[f"{symbol}_prime"] += 1
        else:
            counts[symbol] += 1

    return counts


def classify_atoms_7param_from_wbo(
    smiles: str, wbos: dict[tuple[int, int], float], threshold: float = 1.25
) -> dict[str, int]:
    """Classify atoms into 7 types using xTB Wiberg bond orders.

    Same classification as classify_atoms_7param, but uses xTB-calculated
    Wiberg bond orders instead of RDKit bond types. The wbos dict maps
    (atom_i, atom_j) 0-based index pairs to bond orders.
    """
    mol = smiles_to_mol(smiles)
    counts = {"C": 0, "H": 0, "N": 0, "O": 0, "C_prime": 0, "N_prime": 0, "O_prime": 0}

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        idx = atom.GetIdx()
        has_multiple_bond = any(wbos.get((idx, j), 0.0) > threshold for j in range(mol.GetNumAtoms()) if j != idx)

        if has_multiple_bond:
            counts[f"{symbol}_prime"] += 1
        else:
            counts[symbol] += 1

    return counts


_HYBRID_MAP = {
    Chem.rdchem.HybridizationType.SP3: "sp3",
    Chem.rdchem.HybridizationType.SP2: "sp2",
    Chem.rdchem.HybridizationType.SP: "sp",
}


def classify_atoms_hybrid(smiles: str) -> dict[str, int]:
    """Classify atoms by hybridization for the 10-parameter model.

    Heavy atoms are classified as {element}_sp3, {element}_sp2, or {element}_sp
    based on RDKit's GetHybridization(). Hydrogen is always counted as H.
    Falls back to sp3 if hybridization is unrecognised (e.g. SP3D, UNSPECIFIED).
    """
    mol = smiles_to_mol(smiles)
    counts: dict[str, int] = {
        "C_sp3": 0, "C_sp2": 0, "C_sp": 0,
        "H": 0,
        "N_sp3": 0, "N_sp2": 0, "N_sp": 0,
        "O_sp3": 0, "O_sp2": 0, "O_sp": 0,
    }

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        hyb_label = _HYBRID_MAP.get(atom.GetHybridization(), "sp3")
        key = f"{symbol}_{hyb_label}"
        counts[key] = counts.get(key, 0) + 1

    return counts


def classify_atoms_extended(smiles: str) -> dict[str, int]:
    """Classify atoms by hybridization and H-count for the 15-parameter model.

    Carbon atoms are further split by the number of attached hydrogen atoms
    (sp3 capped at 3 so CH4 → C_sp3_3H; sp2 capped at 2). Nitrogen and
    oxygen are classified by hybridization only.
    """
    mol = smiles_to_mol(smiles)
    counts: dict[str, int] = {
        "C_sp3_3H": 0, "C_sp3_2H": 0, "C_sp3_1H": 0, "C_sp3_0H": 0,
        "C_sp2_2H": 0, "C_sp2_1H": 0, "C_sp2_0H": 0,
        "C_sp": 0,
        "H": 0,
        "N_sp3": 0, "N_sp2": 0, "N_sp": 0,
        "O_sp3": 0, "O_sp2": 0, "O_sp": 0,
    }

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        hyb_label = _HYBRID_MAP.get(atom.GetHybridization(), "sp3")

        if symbol == "C":
            if hyb_label == "sp":
                counts["C_sp"] += 1
            else:
                n_h_raw = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == "H")
                max_h = 3 if hyb_label == "sp3" else 2
                n_h = min(n_h_raw, max_h)
                key = f"C_{hyb_label}_{n_h}H"
                counts[key] += 1
        else:
            key = f"{symbol}_{hyb_label}"
            counts[key] = counts.get(key, 0) + 1

    return counts


def _is_multiply_bonded(atom) -> bool:
    """Check if a heavy atom is in a multiply-bonded environment.

    An atom is considered multiply-bonded if:
    1. It has any bond with order > 1.25 (aromatic, double, triple), OR
    2. It has a nonzero formal charge and is bonded to a neighbor that itself
       has a bond with order > 1.25. This handles charge-separated resonance
       forms like [N+](=O)[O-] where the [O-] is a resonance equivalent of
       a double-bonded oxygen.
    """
    if any(bond.GetBondTypeAsDouble() > 1.25 for bond in atom.GetBonds()):
        return True

    if atom.GetFormalCharge() != 0:
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if any(b.GetBondTypeAsDouble() > 1.25 for b in neighbor.GetBonds()):
                return True

    return False
