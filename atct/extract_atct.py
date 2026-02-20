"""Extract ATcT thermochemical data and convert to deltahf training format.

The ATcT website (https://atct.anl.gov/) does not provide bulk download. Users
must manually export data from the site into a CSV file, then run this script to:
  1. Look up SMILES for each species via the PubChem REST API.
  2. Filter for neutral, closed-shell, gas-phase molecules.
  3. Filter by allowed elements (default: C, H, N, O, F, S, Cl).
  4. Deduplicate against the existing training_data.csv.
  5. Output a CSV ready to append to deltahf/data/training_data.csv.

Usage:
    python scripts/extract_atct.py -i atct_raw.csv -o atct_out.csv

The script auto-detects the ATcT website export format (NAME/DHF/ID columns,
asterisk-prefixed CAS numbers, multi-variant entries) and handles it transparently.

Generic CSV format (with lowercase column names):
    name,formula,dhf_298K_kJ_mol,uncertainty_kJ_mol
    methane,CH4,-74.87,0.32

Optional 'smiles' column skips PubChem lookup for that row:
    name,formula,smiles,dhf_298K_kJ_mol
    methane,CH4,C,-74.87

Optional 'cas' column uses the CAS-specific PubChem endpoint (more reliable):
    name,formula,cas,dhf_298K_kJ_mol
    methane,CH4,74-82-8,-74.87

ATcT website export (NAME/DHF/ID columns, auto-detected):
    NAME,FORMULA,DHF,ID
    methane,CH4 (g),-74.87,74-82-8*0
    hydrogen fluoride,HF (g),-272.68,7664-39-3*0

Output columns match training_data.csv:
    id,name,formula,smiles,exp_dhf_kcal_mol,source,category

Then append to the training data and retrain:
    python -m deltahf fit -i deltahf/data/training_data.csv --model all --kfold 10 -o params.json
"""

import argparse
import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

_KJ_TO_KCAL = 1.0 / 4.184
_PUBCHEM_URL_NAME = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
_PUBCHEM_URL_CAS  = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RN/{cas}/property/IsomericSMILES/JSON"
_DEFAULT_ELEMENTS = {"C", "H", "N", "O", "F", "S", "Cl"}
_ATCT_SOURCE = "ATcT_v1.220"

# Patterns in ATcT FORMULA field that indicate non-gas or open-shell entries
_CONDENSED_PHASE = re.compile(r'\(\s*(cr|l|aq|s)\b', re.I)
_OPEN_SHELL      = re.compile(r'\b(triplet|doublet|radical|quartet|sextet)\b', re.I)


# ---------------------------------------------------------------------------
# ATcT-format detection and normalisation
# ---------------------------------------------------------------------------

def _is_atct_format(fieldnames: list[str]) -> bool:
    """Return True if the CSV looks like an ATcT website export."""
    upper = {f.upper() for f in fieldnames}
    return {"NAME", "DHF", "ID"}.issubset(upper)


def clean_cas(raw: str) -> str:
    """Strip ATcT asterisk prefix and variant suffix from a CAS field.

    ATcT IDs look like '*103613-08-7*11' or '74-82-8*0'.
    Returns the bare CAS string, e.g. '103613-08-7'.
    """
    s = raw.strip().lstrip("*")
    s = re.sub(r"\*\d+$", "", s)
    return s.strip()


def is_atct_custom_id(raw: str) -> bool:
    """Return True if the ATcT ID has a leading *, indicating a non-registry CAS."""
    return raw.strip().startswith("*")


def normalise_atct_rows(rows: list[dict]) -> list[dict]:
    """Convert ATcT export rows to the generic format expected by the pipeline.

    ATcT columns: NAME, FORMULA (contains phase/state info), DHF, ID (CAS+variant).
    Generic columns: name, dhf_298K_kJ_mol, cas, _atct_state (internal filter field).

    Pre-filters:
      - Condensed-phase entries (FORMULA contains '(cr', '(l)', '(aq)', '(s)')
      - Open-shell entries (FORMULA contains 'triplet', 'doublet', 'radical', ...)
      - Entries whose ID has a leading * and appears non-standard (still attempted but flagged)

    Deduplication: for the same cleaned CAS, keep the single lowest-DHF gas-phase
    closed-shell entry (most stable conformer/symmetry variant).
    """
    out = []
    for row in rows:
        name    = row.get("NAME", row.get("name", "")).strip()
        formula = row.get("FORMULA", row.get("formula", "")).strip()
        raw_id  = row.get("ID", row.get("id", "")).strip()
        try:
            dhf = float(row.get("DHF", row.get("dhf_298K_kJ_mol", 0)))
        except ValueError:
            continue

        # Phase filter: skip condensed-phase entries
        if _CONDENSED_PHASE.search(formula):
            continue

        # Open-shell filter: skip triplets, doublets, radicals (cheap pre-filter)
        if _OPEN_SHELL.search(formula):
            continue

        cas = clean_cas(raw_id) if raw_id else ""
        # Extract the raw formula prefix for later cross-validation (text before first ' (')
        atct_formula = re.split(r"\s*\(", formula)[0].strip()
        out.append({
            "name": name,
            "dhf_298K_kJ_mol": dhf,
            "cas": cas,
            "_custom_id": is_atct_custom_id(raw_id),  # flag for fallback to name lookup
            "_atct_formula": atct_formula,             # for formula cross-check
        })

    # Deduplicate: for same CAS, keep lowest DHF (most stable gas-phase entry)
    best: dict[str, dict] = {}
    no_cas = []
    for row in out:
        cas = row["cas"]
        if not cas:
            no_cas.append(row)
        elif cas not in best or row["dhf_298K_kJ_mol"] < best[cas]["dhf_298K_kJ_mol"]:
            best[cas] = row

    return list(best.values()) + no_cas


# ---------------------------------------------------------------------------
# ATcT HTML SMILES extraction
# ---------------------------------------------------------------------------

def load_html_smiles(html_path: Path) -> dict[str, str]:
    """Parse the ATcT HTML page and return a CAS → SMILES mapping.

    The HTML contains a JS autocomplete array with entries like:
        { label: "SMILES=CC", value: "74-84-0"}
    When multiple SMILES exist for the same CAS (e.g. resonance forms or vdW
    complexes), prefer single-fragment SMILES (no '.') and take the first.
    """
    with html_path.open(encoding="utf-8", errors="replace") as f:
        content = f.read()
    pairs = re.findall(r'\{ label: "SMILES=([^"]+)", value: "([^"]+)"\}', content)
    cas_to_smiles: dict[str, str] = {}
    for smiles, cas in pairs:
        if cas not in cas_to_smiles:
            cas_to_smiles[cas] = smiles
        elif "." in cas_to_smiles[cas] and "." not in smiles:
            # Prefer single-fragment SMILES over complex/mixture SMILES
            cas_to_smiles[cas] = smiles
    return cas_to_smiles


# ---------------------------------------------------------------------------
# PubChem lookup
# ---------------------------------------------------------------------------

def _pubchem_fetch(url: str, retries: int, delay: float) -> str | None:
    """Fetch IsomericSMILES from a PubChem REST URL, with retry/backoff."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=10) as resp:
                data = json.loads(resp.read())
            props = data["PropertyTable"]["Properties"]
            # PubChem returns the value under "SMILES" regardless of the property requested
            return props[0].get("IsomericSMILES") or props[0].get("SMILES")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None  # not found — no point retrying
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
    return None


def pubchem_cas_to_smiles(cas: str, retries: int = 3, delay: float = 1.0) -> str | None:
    """Return IsomericSMILES from PubChem for a CAS registry number."""
    url = _PUBCHEM_URL_CAS.format(cas=urllib.parse.quote(cas))
    return _pubchem_fetch(url, retries, delay)


def pubchem_name_to_smiles(name: str, retries: int = 3, delay: float = 1.0) -> str | None:
    """Return IsomericSMILES from PubChem for a chemical name."""
    url = _PUBCHEM_URL_NAME.format(name=urllib.parse.quote(name))
    return _pubchem_fetch(url, retries, delay)


# ---------------------------------------------------------------------------
# Molecule validation
# ---------------------------------------------------------------------------

def is_neutral_closed_shell(smiles: str) -> bool:
    """Return True if the molecule is neutral and has no radical electrons."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    if Chem.GetFormalCharge(mol) != 0:
        return False
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons() != 0:
            return False
    return True


def elements_are_supported(smiles: str, allowed: set[str]) -> bool:
    """Return True if all heavy atoms in the molecule are in *allowed*."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol_h = Chem.AddHs(mol)
    symbols = {atom.GetSymbol() for atom in mol_h.GetAtoms()}
    return symbols.issubset(allowed)


def canonical_smiles(smiles: str) -> str | None:
    """Return RDKit canonical SMILES, or None if the SMILES is invalid."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def molecular_formula(smiles: str) -> str | None:
    """Return the molecular formula (with explicit H) for *smiles*."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol_h = Chem.AddHs(mol)
    return rdMolDescriptors.CalcMolFormula(mol_h)


def _parse_formula(formula: str) -> dict[str, int]:
    """Parse a molecular formula string into element counts."""
    counts: dict[str, int] = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if element:
            counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    return counts


def formula_matches(smiles: str, expected_formula: str) -> bool:
    """Return True if the RDKit formula from *smiles* matches *expected_formula*.

    Comparison is element-wise (ignores ordering). Returns True if
    *expected_formula* is empty or unparseable (no false rejections).
    """
    if not expected_formula:
        return True
    expected = _parse_formula(expected_formula)
    if not expected:
        return True
    actual_str = molecular_formula(smiles)
    if actual_str is None:
        return True
    actual = _parse_formula(actual_str)
    return actual == expected


# ---------------------------------------------------------------------------
# Category assignment
# ---------------------------------------------------------------------------

def assign_category(smiles: str) -> str:
    """Assign a category label based on elemental composition."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unknown"
    symbols = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if "F" in symbols:
        return "fluorinated"
    if "S" in symbols:
        return "sulfur"
    if "Cl" in symbols:
        return "chlorinated"
    if symbols <= {"C", "H"}:
        return "hydrocarbon"
    if "N" in symbols or "O" in symbols:
        return "small_CHNO"
    return "other"


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def load_existing_canonical_smiles(csv_path: Path) -> set[str]:
    """Read canonical SMILES already present in training_data.csv."""
    existing: set[str] = set()
    if not csv_path.exists():
        return existing
    with csv_path.open(encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            canon = canonical_smiles(row["smiles"])
            if canon:
                existing.add(canon)
    return existing


def next_id(csv_path: Path) -> int:
    """Return the next sequential id for appending to training_data.csv."""
    if not csv_path.exists():
        return 1
    with csv_path.open(encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        ids = [int(row["id"]) for row in reader if row.get("id", "").isdigit()]
    return (max(ids) + 1) if ids else 1


def load_smiles_cache(cache_path: Path) -> dict[str, str | None]:
    """Load previously fetched name/CAS→SMILES results from a JSON cache."""
    if cache_path.exists():
        with cache_path.open() as f:
            return json.load(f)
    return {}


def save_smiles_cache(cache_path: Path, cache: dict[str, str | None]) -> None:
    with cache_path.open("w") as f:
        json.dump(cache, f, indent=2)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert an ATcT manual-export CSV to deltahf training format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("-i", "--input", required=True,
                   help="Input CSV. ATcT website export (NAME/DHF/ID columns) is auto-detected. "
                        "Optional 'smiles' column skips PubChem; optional 'cas' column uses CAS lookup.")
    p.add_argument("-o", "--output", default="atct_out.csv",
                   help="Output CSV to append to training_data.csv (default: atct_out.csv)")
    p.add_argument(
        "--elements",
        default=",".join(sorted(_DEFAULT_ELEMENTS)),
        help="Comma-separated allowed elements (default: C,Cl,F,H,N,O,S)",
    )
    p.add_argument(
        "--existing",
        default="deltahf/data/training_data.csv",
        help="Existing training CSV to check for duplicates (default: deltahf/data/training_data.csv)",
    )
    p.add_argument(
        "--smiles-cache",
        default=".atct_smiles_cache.json",
        help="JSON file caching PubChem lookups for reruns (default: .atct_smiles_cache.json)",
    )
    p.add_argument(
        "--html",
        default=None,
        help="Path to the ATcT HTML page (e.g. 'ATcT Thermochemical Values ver. 1.220.html'). "
             "When provided, SMILES are read directly from the page (no PubChem calls needed).",
    )
    p.add_argument(
        "--skip-pubchem",
        action="store_true",
        help="Skip PubChem lookups; only process rows already in the cache or with inline SMILES",
    )
    p.add_argument(
        "--start-id",
        type=int,
        default=None,
        help="Override the starting id (default: auto-detected from --existing)",
    )
    p.add_argument(
        "--dhf-column",
        default="dhf_298K_kJ_mol",
        help="Column name for ΔHf° in kJ/mol for non-ATcT CSVs (default: dhf_298K_kJ_mol)",
    )
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    allowed_elements = set(args.elements.split(","))
    input_path = Path(args.input)
    output_path = Path(args.output)
    existing_path = Path(args.existing)
    cache_path = Path(args.smiles_cache)

    # Load existing canonical SMILES for deduplication
    existing_canon = load_existing_canonical_smiles(existing_path)
    print(f"Existing training molecules: {len(existing_canon)}")

    # Load ATcT HTML SMILES map (CAS → SMILES) when provided
    html_smiles: dict[str, str] = {}
    if args.html:
        html_path = Path(args.html)
        html_smiles = load_html_smiles(html_path)
        print(f"Loaded {len(html_smiles)} SMILES from {html_path.name}")

    # Load SMILES cache
    smiles_cache = load_smiles_cache(cache_path)

    # Read input (utf-8-sig handles optional BOM)
    with input_path.open(encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        raw_rows = list(reader)
        fieldnames = reader.fieldnames or []

    # Auto-detect and normalise ATcT website export format
    atct_fmt = _is_atct_format(fieldnames)
    if atct_fmt:
        rows = normalise_atct_rows(raw_rows)
        print(f"ATcT format detected. Input rows: {len(raw_rows)} → {len(rows)} after phase/spin pre-filter")
    else:
        rows = raw_rows
        print(f"Input rows: {len(rows)}")

    next_molecule_id = args.start_id if args.start_id is not None else next_id(existing_path)

    output_rows: list[dict] = []
    stats = {"no_smiles": 0, "invalid_smiles": 0, "formula_mismatch": 0,
             "charged_or_radical": 0, "anion_cation_name": 0, "multi_fragment": 0,
             "diatomic": 0, "unsupported_elements": 0, "duplicate": 0, "accepted": 0}

    for i, row in enumerate(rows, 1):
        name = row.get("name", row.get("NAME", "")).strip()
        dhf_kj = float(row.get("dhf_298K_kJ_mol", row.get("DHF", 0)))
        dhf_kcal = dhf_kj * _KJ_TO_KCAL

        # Reject entries whose name contains problematic keywords:
        #   anion/cation  — net-neutral zwitterions that slip through the charge filter
        #   diyl          — biradicals encoded as zwitterions by ATcT
        if re.search(r"\b(anion|cation|diyl)\b", name, re.I):
            stats["anion_cation_name"] += 1
            if args.verbose:
                print(f"    SKIP (anion/cation/diyl name): {name!r}")
            continue

        # SMILES resolution priority:
        #   1. inline 'smiles' column      (no API call)
        #   2. ATcT HTML map (CAS lookup)  (no API call)
        #   3. cache hit                   (no API call)
        #   4. CAS → PubChem xref/RN      (standard CAS only)
        #   5. name → PubChem name         (fallback)
        inline_smiles = row.get("smiles", "").strip()
        cas = row.get("cas", "").strip()
        is_custom = row.get("_custom_id", False)  # ATcT-specific non-registry CAS
        cache_key = cas if cas else name

        if inline_smiles:
            raw_smiles = inline_smiles
        elif cas and cas in html_smiles:
            raw_smiles = html_smiles[cas]
            if args.verbose:
                print(f"  [{i}/{len(rows)}] HTML lookup: {cas} ({name!r}) -> {raw_smiles!r}")
        elif cache_key in smiles_cache:
            raw_smiles = smiles_cache[cache_key]
        elif args.skip_pubchem:
            raw_smiles = None
        elif cas and not is_custom:
            if args.verbose:
                print(f"  [{i}/{len(rows)}] PubChem CAS lookup: {cas} ({name!r})")
            raw_smiles = pubchem_cas_to_smiles(cas)
            if raw_smiles is None and name:
                if args.verbose:
                    print(f"    CAS not found, trying name lookup: {name!r}")
                raw_smiles = pubchem_name_to_smiles(name)
            smiles_cache[cache_key] = raw_smiles
            time.sleep(0.22)
        elif name:
            if args.verbose:
                print(f"  [{i}/{len(rows)}] PubChem name lookup: {name!r}")
            raw_smiles = pubchem_name_to_smiles(name)
            smiles_cache[cache_key] = raw_smiles
            time.sleep(0.22)
        else:
            raw_smiles = None

        if raw_smiles is None:
            stats["no_smiles"] += 1
            if args.verbose:
                print(f"    SKIP (no SMILES): {name!r}")
            continue

        canon = canonical_smiles(raw_smiles)
        if canon is None:
            stats["invalid_smiles"] += 1
            if args.verbose:
                print(f"    SKIP (invalid SMILES): {name!r} -> {raw_smiles!r}")
            continue

        # Cross-check formula when available (catches PubChem name→wrong compound)
        atct_formula = row.get("_atct_formula", "")
        if atct_formula and not formula_matches(canon, atct_formula):
            stats["formula_mismatch"] += 1
            if args.verbose:
                actual = molecular_formula(canon) or "?"
                print(f"    SKIP (formula mismatch): {name!r} ATcT={atct_formula!r} RDKit={actual!r}")
            continue

        if not is_neutral_closed_shell(canon):
            stats["charged_or_radical"] += 1
            if args.verbose:
                print(f"    SKIP (charged/radical): {name!r}")
            continue

        # Reject van der Waals complexes, dimers, and mixtures
        mol = Chem.MolFromSmiles(canon)
        if mol is not None and len(Chem.GetMolFrags(mol)) > 1:
            stats["multi_fragment"] += 1
            if args.verbose:
                print(f"    SKIP (multi-fragment): {name!r}")
            continue

        # Reject diatomics (CO, CS, N2, …) — too few atoms to parameterise usefully
        if mol is not None and mol.GetNumHeavyAtoms() < 2:
            stats["diatomic"] += 1
            if args.verbose:
                print(f"    SKIP (diatomic): {name!r}")
            continue

        if not elements_are_supported(canon, allowed_elements):
            stats["unsupported_elements"] += 1
            if args.verbose:
                print(f"    SKIP (unsupported elements): {name!r}")
            continue

        if canon in existing_canon:
            stats["duplicate"] += 1
            if args.verbose:
                print(f"    SKIP (duplicate): {name!r}")
            continue

        # Accept
        actual_formula = molecular_formula(canon) or ""
        category = assign_category(canon)
        output_rows.append({
            "id": next_molecule_id,
            "name": name,
            "formula": actual_formula,
            "smiles": canon,
            "exp_dhf_kcal_mol": round(dhf_kcal, 6),
            "source": _ATCT_SOURCE,
            "category": category,
        })
        existing_canon.add(canon)  # prevent intra-batch duplicates
        next_molecule_id += 1
        stats["accepted"] += 1

    # Save SMILES cache after every run (even if 0 accepted, caches failures)
    save_smiles_cache(cache_path, smiles_cache)
    print(f"SMILES cache saved to {cache_path} ({len(smiles_cache)} entries)")

    # Write output
    fieldnames_out = ["id", "name", "formula", "smiles", "exp_dhf_kcal_mol", "source", "category"]
    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames_out)
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"\nResults:")
    print(f"  Accepted:              {stats['accepted']}")
    print(f"  No SMILES found:       {stats['no_smiles']}")
    print(f"  Invalid SMILES:        {stats['invalid_smiles']}")
    print(f"  Formula mismatch:      {stats['formula_mismatch']}")
    print(f"  Charged / radical:     {stats['charged_or_radical']}")
    print(f"  Anion/cation name:     {stats['anion_cation_name']}")
    print(f"  Multi-fragment:        {stats['multi_fragment']}")
    print(f"  Diatomic:              {stats['diatomic']}")
    print(f"  Unsupported elements:  {stats['unsupported_elements']}")
    print(f"  Duplicates:            {stats['duplicate']}")
    print(f"\nOutput written to {output_path}")
    if stats["accepted"]:
        print(f"\nTo append to training data:")
        print(f"  tail -n +2 {output_path} >> {existing_path}")
        print(f"\nThen refit:")
        print(f"  python -m deltahf fit -i {existing_path} --model all --kfold 10 -o params.json")


if __name__ == "__main__":
    main()
