"""Optional JSON-backed cache for geometry optimization results (xTB or MLIP)."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from rdkit import Chem


@dataclass
class CachedResult:
    """Cached geometry optimization result for a single molecule."""

    canonical_smiles: str
    xtb_energy: float | None
    xtb_energy_kcal: float | None
    n_conformers: int
    n_conformers_optimized: int
    n_conformers_isomerized: int
    gfn_level: int
    charge: int
    gxtb_energy: float | None = None
    gxtb_energy_kcal: float | None = None
    mlip_energy: float | None = None
    mlip_energy_kcal: float | None = None
    optimizer: str = "xtb"  # Optimizer used to generate this cache entry


class ResultCache:
    """JSON-backed cache for geometry optimization results."""

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
        optimizer: str = "xtb",
        charge: int = 0,
    ) -> CachedResult | None:
        key = self._make_key(smiles)
        entry = self._data.get(key)
        if entry is None:
            return None
        # "optimizer" field added later; old xTB-only cache entries default to "xtb"
        cached_optimizer = entry.get("optimizer", "xtb")
        if (
            entry.get("n_conformers") != n_conformers
            or cached_optimizer != optimizer
            or entry.get("charge") != charge
        ):
            return None
        # Use .get() with defaults so old cache entries without new fields still load
        return CachedResult(
            canonical_smiles=entry["canonical_smiles"],
            xtb_energy=entry.get("xtb_energy"),
            xtb_energy_kcal=entry.get("xtb_energy_kcal"),
            n_conformers=entry["n_conformers"],
            n_conformers_optimized=entry["n_conformers_optimized"],
            n_conformers_isomerized=entry["n_conformers_isomerized"],
            gfn_level=entry.get("gfn_level", 2),
            charge=entry["charge"],
            gxtb_energy=entry.get("gxtb_energy"),
            gxtb_energy_kcal=entry.get("gxtb_energy_kcal"),
            # "mlip_energy" renamed from "uma_energy"; fall back for old cache files
            mlip_energy=entry.get("mlip_energy", entry.get("uma_energy")),
            mlip_energy_kcal=entry.get("mlip_energy_kcal", entry.get("uma_energy_kcal")),
            optimizer=cached_optimizer,
        )

    def store(self, result: CachedResult) -> None:
        key = self._make_key(result.canonical_smiles)
        self._data[key] = asdict(result)
