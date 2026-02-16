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
