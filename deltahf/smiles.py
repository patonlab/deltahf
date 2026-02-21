"""SMILES parsing, atom counting, and bond order classification."""

from rdkit import Chem

from deltahf.constants import (
    PARAM_NAMES_4,
    PARAM_NAMES_7,
    PARAM_NAMES_BONDORDER,
    PARAM_NAMES_BONDORDER_AR,
    PARAM_NAMES_BONDORDER_EXT,
    PARAM_NAMES_EXTENDED,
    PARAM_NAMES_HYBRID,
    PARAM_NAMES_NEIGHBOUR,
)

SUPPORTED_ELEMENTS = {"C", "H", "N", "O", "F", "S", "Cl"}


def _zero_counts(param_names: list[str]) -> dict[str, int]:
    """Return a zeroed counts dict for the given parameter names."""
    return dict.fromkeys(param_names, 0)


def _h_count(atom) -> int:
    """Count explicit hydrogen neighbours of an atom."""
    return sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == "H")


def validate_elements(smiles: str, supported_elements: set[str] | None = None) -> None:
    """Validate that molecule contains only supported elements.

    Args:
        smiles: SMILES string to validate
        supported_elements: Set of allowed element symbols (defaults to CHNO)

    Raises:
        ValueError: If molecule contains unsupported elements
    """
    if supported_elements is None:
        supported_elements = SUPPORTED_ELEMENTS

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    mol = Chem.AddHs(mol)
    present_elements = {atom.GetSymbol() for atom in mol.GetAtoms()}
    unsupported = present_elements - supported_elements

    if unsupported:
        raise ValueError(
            f"Molecule contains unsupported elements: {sorted(unsupported)}. "
            f"deltahf only supports: {sorted(supported_elements)}. "
            f"SMILES: {smiles}"
        )


def heavy_atom_count(smiles: str) -> int | None:
    """Return the number of heavy (non-hydrogen) atoms, or None if unparseable."""
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol is not None else None


def total_atom_count(smiles: str) -> int | None:
    """Return the total atom count including hydrogens, or None if unparseable."""
    mol = Chem.MolFromSmiles(smiles)
    return Chem.AddHs(mol).GetNumAtoms() if mol is not None else None


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
    """Count C, H, N, O, F, S, Cl atoms for the element-count model."""
    mol = smiles_to_mol(smiles)
    counts = _zero_counts(PARAM_NAMES_4)
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
    counts = _zero_counts(PARAM_NAMES_7)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O", "F", "S", "Cl"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        # F and Cl form only single bonds in organic chemistry — never primed
        if symbol in ("F", "Cl"):
            counts[symbol] += 1
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
    counts = _zero_counts(PARAM_NAMES_7)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O", "F", "S", "Cl"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        if symbol in ("F", "Cl"):
            counts[symbol] += 1
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
    """Classify atoms by hybridization.

    Heavy atoms are classified as {element}_sp3, {element}_sp2, or {element}_sp
    based on RDKit's GetHybridization(). Hydrogen is always counted as H.
    Falls back to sp3 if hybridization is unrecognised (e.g. SP3D, UNSPECIFIED).
    """
    mol = smiles_to_mol(smiles)
    counts = _zero_counts(PARAM_NAMES_HYBRID)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O", "F", "S", "Cl"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        # F and Cl are always sp3 in organic contexts
        if symbol in ("F", "Cl"):
            counts[f"{symbol}_sp3"] += 1
            continue

        hyb_label = _HYBRID_MAP.get(atom.GetHybridization(), "sp3")
        key = f"{symbol}_{hyb_label}"
        counts[key] = counts.get(key, 0) + 1

    return counts


def classify_atoms_extended(smiles: str) -> dict[str, int]:
    """Classify atoms by hybridization and H-count (extended model).

    Carbon atoms are further split by the number of attached hydrogen atoms
    (sp3 capped at 3 so CH4 → C_sp3_3H; sp2 capped at 2). Nitrogen and
    oxygen are classified by hybridization only.
    """
    mol = smiles_to_mol(smiles)
    counts = _zero_counts(PARAM_NAMES_EXTENDED)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O", "F", "S", "Cl"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        # F and Cl: single type each (no H-count or hybridization distinction needed)
        if symbol in ("F", "Cl"):
            counts[symbol] += 1
            continue

        hyb_label = _HYBRID_MAP.get(atom.GetHybridization(), "sp3")

        if symbol == "C":
            if hyb_label == "sp":
                counts["C_sp"] += 1
            else:
                n_h_raw = _h_count(atom)
                max_h = 3 if hyb_label == "sp3" else 2
                n_h = min(n_h_raw, max_h)
                key = f"C_{hyb_label}_{n_h}H"
                counts[key] += 1
        elif symbol == "S":
            # S: sp3 (thioether) or sp2 (thiophene, sulfoxide) — sp is rare
            key = "S_sp2" if hyb_label == "sp2" else "S_sp3"
            counts[key] += 1
        else:
            key = f"{symbol}_{hyb_label}"
            counts[key] = counts.get(key, 0) + 1

    return counts


def classify_atoms_bondorder(smiles: str) -> dict[str, int]:
    """Classify atoms by their highest bond order in the Kekulized structure.

    Assigns atoms based
    on the maximum bond order of any bond they participate in (after Kekulization),
    rather than RDKit's hybridization. Suffix _1/_2/_3 corresponds to single/double/triple.

    Key differences from hybrid: amide N and pyrrole-type N/O have no double bond
    after Kekulization so fall into _1 (not _2). Aromatic carbons alternate between
    _1 and _2 depending on which resonance structure Kekulization selects.
    """
    mol = smiles_to_mol(smiles)
    mol = Chem.AddHs(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)

    counts = _zero_counts(PARAM_NAMES_BONDORDER)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O", "F", "S", "Cl"):
            continue
        if symbol == "H":
            counts["H"] += 1
            continue
        max_order = max((int(b.GetBondTypeAsDouble()) for b in atom.GetBonds()), default=1)
        max_order = min(max_order, 3)
        # F and Cl are capped at 1 in standard organic molecules
        if symbol in ("F", "Cl"):
            counts[f"{symbol}_1"] += 1
        else:
            counts[f"{symbol}_{max_order}"] += 1

    return counts


def classify_atoms_bondorder_ext(smiles: str) -> dict[str, int]:
    """Classify atoms by highest bond order + H-count for carbon (bondorder_ext model).

    The bond-order analogue of classify_atoms_extended(): carbon atoms are split by
    their maximum bond order in the Kekulized structure (_1/_2/_3) and then further
    by the number of attached hydrogens (capped at 3/2/1 respectively). Nitrogen and
    oxygen are classified by bond order only (_1/_2/_3). Hydrogen is always H.

    Caps per bond-order tier:
      C_1: up to 3H (sp3-like)   C_2: up to 2H (sp2-like)   C_3: up to 1H (sp-like)
    """
    mol = smiles_to_mol(smiles)
    mol = Chem.AddHs(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)

    counts = _zero_counts(PARAM_NAMES_BONDORDER_EXT)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O", "F", "S", "Cl"):
            continue
        if symbol == "H":
            counts["H"] += 1
            continue

        # F and Cl: always single bond, no H-count needed
        if symbol in ("F", "Cl"):
            counts[f"{symbol}_1"] += 1
            continue

        max_order = min(max((int(b.GetBondTypeAsDouble()) for b in atom.GetBonds()), default=1), 3)

        if symbol == "S":
            counts[f"S_{max_order}"] += 1
            continue

        if symbol != "C":
            counts[f"{symbol}_{max_order}"] += 1
            continue

        n_h = _h_count(atom)
        max_h = {1: 3, 2: 2, 3: 1}[max_order]
        counts[f"C_{max_order}_{min(n_h, max_h)}H"] += 1

    return counts


def classify_atoms_bondorder_ar(smiles: str) -> dict[str, int]:
    """Classify atoms by highest bond order, preserving aromatic bonds as a distinct type.

    Like bondorder but without Kekulization: aromatic bonds remain at order 1.5 and
    atoms whose highest bond is aromatic are labelled _ar rather than being split
    arbitrarily into _1/_2. This gives all equivalent aromatic atoms the same label.

    Suffixes: _1 = single only, _ar = aromatic (highest), _2 = double (non-aromatic), _3 = triple.
    """
    mol = smiles_to_mol(smiles)
    mol = Chem.AddHs(mol)

    counts = _zero_counts(PARAM_NAMES_BONDORDER_AR)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O", "F", "S", "Cl"):
            continue
        if symbol == "H":
            counts["H"] += 1
            continue
        # F and Cl are never aromatic
        if symbol in ("F", "Cl"):
            counts[f"{symbol}_1"] += 1
            continue
        max_order = max((b.GetBondTypeAsDouble() for b in atom.GetBonds()), default=1.0)
        if max_order >= 3.0:
            label = "3"
        elif max_order >= 2.0:
            label = "2"
        elif max_order >= 1.5:
            label = "ar"
        else:
            label = "1"
        counts[f"{symbol}_{label}"] += 1

    return counts


def classify_atoms_neighbour(smiles: str) -> dict[str, int]:
    """Classify atoms by hybridization, H-count (C only), and heavy-atom neighbour type.

    Carbon is classified identically to the extended model (hybridization + H-count).
    Nitrogen and oxygen are further split by the highest-priority heavy-atom neighbour
    directly bonded to them, using priority O > N > C. This separates chemically
    distinct environments such as nitro N (bonded to O) from pyridine N (bonded to C)
    and carbonyl O (bonded to C) from nitro O (bonded to N).
    """
    mol = smiles_to_mol(smiles)
    counts = _zero_counts(PARAM_NAMES_NEIGHBOUR)

    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H", "N", "O", "F", "S", "Cl"):
            continue

        if symbol == "H":
            counts["H"] += 1
            continue

        if symbol in ("F", "Cl"):
            counts[symbol] += 1
            continue

        if symbol == "S":
            hyb_label = _HYBRID_MAP.get(atom.GetHybridization(), "sp3")
            key = "S_sp2" if hyb_label == "sp2" else "S_sp3"
            counts[key] += 1
            continue

        hyb_label = _HYBRID_MAP.get(atom.GetHybridization(), "sp3")

        if symbol == "C":
            if hyb_label == "sp":
                counts["C_sp"] += 1
            else:
                n_h_raw = _h_count(atom)
                max_h = 3 if hyb_label == "sp3" else 2
                n_h = min(n_h_raw, max_h)
                counts[f"C_{hyb_label}_{n_h}H"] += 1
        else:
            heavy = {n.GetSymbol() for n in atom.GetNeighbors() if n.GetSymbol() != "H"}
            if "O" in heavy:
                nbr = "O"
            elif "N" in heavy:
                nbr = "N"
            else:
                nbr = "C"
            key = f"{symbol}_{hyb_label}_{nbr}"
            if key in counts:
                counts[key] += 1

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
