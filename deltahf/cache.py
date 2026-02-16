"""Optional JSON-backed cache for xTB optimization results."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem


@dataclass
class CachedResult:
    """Cached xTB result for a single molecule."""

    canonical_smiles: str
    xtb_energy: float
    xtb_energy_kcal: float
    n_conformers: int
    n_conformers_optimized: int
    n_conformers_isomerized: int
    gfn_level: int
    charge: int
    gxtb_energy: float | None = None
    gxtb_energy_kcal: float | None = None


class ResultCache:
    """JSON-backed cache for xTB optimization results."""

    def __init__(self, cache_dir: Path):
        self.cache_dir = cache_dir
        self.cache_file = cache_dir / "cache.json"
        self._data: dict[str, dict] = {}
        self._load()

    def _load(self) -> None:
        if self.cache_file.exists():
            with open(self.cache_file) as f:
                self._data = json.load(f)

    def save(self) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        with open(self.cache_file, "w") as f:
            json.dump(self._data, f, indent=2)

    @staticmethod
    def _make_key(smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol)

    def lookup(
        self,
        smiles: str,
        n_conformers: int,
        gfn: int = 2,
        charge: int = 0,
    ) -> CachedResult | None:
        key = self._make_key(smiles)
        entry = self._data.get(key)
        if entry is None:
            return None
        if (
            entry.get("n_conformers") != n_conformers
            or entry.get("gfn_level") != gfn
            or entry.get("charge") != charge
        ):
            return None
        return CachedResult(**entry)

    def store(self, result: CachedResult) -> None:
        key = self._make_key(result.canonical_smiles)
        self._data[key] = {
            "canonical_smiles": result.canonical_smiles,
            "xtb_energy": result.xtb_energy,
            "xtb_energy_kcal": result.xtb_energy_kcal,
            "n_conformers": result.n_conformers,
            "n_conformers_optimized": result.n_conformers_optimized,
            "n_conformers_isomerized": result.n_conformers_isomerized,
            "gfn_level": result.gfn_level,
            "charge": result.charge,
            "gxtb_energy": result.gxtb_energy,
            "gxtb_energy_kcal": result.gxtb_energy_kcal,
        }
