"""Tests for geometry optimization result caching (xTB and MLIP)."""

from unittest.mock import patch

from deltahf.cache import CachedResult, ResultCache


def _make_cached_result(**overrides):
    defaults = {
        "canonical_smiles": "CCO",
        "xtb_energy": -5.0,
        "xtb_energy_kcal": -3137.547,
        "n_conformers": 5,
        "n_conformers_optimized": 5,
        "n_conformers_isomerized": 0,
        "gfn_level": 2,
        "charge": 0,
        "gxtb_energy": None,
        "gxtb_energy_kcal": None,
    }
    defaults.update(overrides)
    return CachedResult(**defaults)


class TestResultCache:
    def test_store_and_lookup(self, tmp_path):
        cache = ResultCache(tmp_path)
        result = _make_cached_result()
        cache.store(result)
        found = cache.lookup("CCO", n_conformers=5)
        assert found is not None
        assert found.xtb_energy == -5.0
        assert found.canonical_smiles == "CCO"

    def test_lookup_missing_returns_none(self, tmp_path):
        cache = ResultCache(tmp_path)
        assert cache.lookup("CCO", n_conformers=5) is None

    def test_lookup_wrong_n_conformers_returns_none(self, tmp_path):
        cache = ResultCache(tmp_path)
        cache.store(_make_cached_result(n_conformers=5))
        assert cache.lookup("CCO", n_conformers=10) is None

    def test_lookup_wrong_optimizer_returns_none(self, tmp_path):
        cache = ResultCache(tmp_path)
        cache.store(_make_cached_result(optimizer="xtb"))
        assert cache.lookup("CCO", n_conformers=5, optimizer="uma") is None

    def test_lookup_wrong_charge_returns_none(self, tmp_path):
        cache = ResultCache(tmp_path)
        cache.store(_make_cached_result(charge=0))
        assert cache.lookup("CCO", n_conformers=5, charge=-1) is None

    def test_canonical_smiles_normalization(self, tmp_path):
        """'OCC' and 'CCO' canonicalize to the same SMILES."""
        cache = ResultCache(tmp_path)
        cache.store(_make_cached_result(canonical_smiles="CCO"))
        found = cache.lookup("OCC", n_conformers=5)
        assert found is not None
        assert found.xtb_energy == -5.0

    def test_persistence_across_instances(self, tmp_path):
        cache1 = ResultCache(tmp_path)
        cache1.store(_make_cached_result())
        cache1.save()

        cache2 = ResultCache(tmp_path)
        found = cache2.lookup("CCO", n_conformers=5)
        assert found is not None
        assert found.xtb_energy == -5.0

    def test_cache_dir_created_if_missing(self, tmp_path):
        new_dir = tmp_path / "subdir" / "cache"
        cache = ResultCache(new_dir)
        cache.store(_make_cached_result())
        cache.save()
        assert (new_dir / "cache.json").exists()

    def test_save_creates_valid_json(self, tmp_path):
        import json

        cache = ResultCache(tmp_path)
        cache.store(_make_cached_result())
        cache.save()
        with open(tmp_path / "cache.json") as f:
            data = json.load(f)
        assert len(data) == 1

    def test_overwrite_existing_entry(self, tmp_path):
        cache = ResultCache(tmp_path)
        cache.store(_make_cached_result(xtb_energy=-5.0))
        cache.store(_make_cached_result(xtb_energy=-6.0))
        found = cache.lookup("CCO", n_conformers=5)
        assert found.xtb_energy == -6.0

    def test_gxtb_fields_stored_and_retrieved(self, tmp_path):
        """Test that gxtb energy fields are properly cached."""
        cache = ResultCache(tmp_path)
        result = _make_cached_result(
            gxtb_energy=-4.5,
            gxtb_energy_kcal=-2823.792,
        )
        cache.store(result)
        cache.save()

        cache2 = ResultCache(tmp_path)
        found = cache2.lookup("CCO", n_conformers=5)
        assert found is not None
        assert found.gxtb_energy == -4.5
        assert found.gxtb_energy_kcal == -2823.792

    def test_backward_compatibility_without_gxtb(self, tmp_path):
        """Cache files without gxtb fields should load with None values."""
        import json

        # Simulate old cache file without gxtb fields
        old_cache_data = {
            "CCO": {
                "canonical_smiles": "CCO",
                "xtb_energy": -5.0,
                "xtb_energy_kcal": -3137.547,
                "n_conformers": 5,
                "n_conformers_optimized": 5,
                "n_conformers_isomerized": 0,
                "gfn_level": 2,
                "charge": 0,
            }
        }
        cache_file = tmp_path / "cache.json"
        with open(cache_file, "w") as f:
            json.dump(old_cache_data, f)

        cache = ResultCache(tmp_path)
        found = cache.lookup("CCO", n_conformers=5)
        assert found is not None
        assert found.xtb_energy == -5.0
        assert found.gxtb_energy is None
        assert found.gxtb_energy_kcal is None


class TestCacheInPipeline:
    def test_cache_prevents_xtb_call(self, tmp_path):
        """When cache has a hit, xTB should not be called."""
        cache = ResultCache(tmp_path)
        cache.store(_make_cached_result(canonical_smiles="C"))

        with patch("deltahf.pipeline.run_xtb_optimization") as mock_xtb:
            from deltahf.pipeline import process_molecule

            result = process_molecule("C", n_conformers=5, cache=cache)

        mock_xtb.assert_not_called()
        assert result.xtb_energy == -5.0
        assert result.error is None

    def test_cache_populated_after_xtb(self, tmp_path):
        """After xTB runs, result should be stored in cache."""
        from deltahf.pipeline import process_molecule
        from deltahf.xtb import XtbResult

        cache = ResultCache(tmp_path)

        def mock_xtb(xyz_path, **kwargs):
            return XtbResult(
                energy=-10.0,
                optimized_xyz_path=None,
                converged=True,
                stdout="",
                wbo_path=None,
            )

        with (
            patch("deltahf.pipeline.run_xtb_optimization", side_effect=mock_xtb),
            patch("deltahf.pipeline.check_connectivity", return_value=True),
        ):
            process_molecule("C", n_conformers=2, work_dir=tmp_path / "work", cache=cache)

        found = cache.lookup("C", n_conformers=2)
        assert found is not None
        assert found.xtb_energy == -10.0

    def test_cache_skipped_with_xtb_wbos(self, tmp_path):
        """When use_xtb_wbos=True, cache should not be used."""
        cache = ResultCache(tmp_path)
        cache.store(_make_cached_result(canonical_smiles="C"))

        from deltahf.pipeline import process_molecule
        from deltahf.xtb import XtbResult

        def mock_xtb(xyz_path, **kwargs):
            return XtbResult(
                energy=-10.0,
                optimized_xyz_path=None,
                converged=True,
                stdout="",
                wbo_path=None,
            )

        with (
            patch("deltahf.pipeline.run_xtb_optimization", side_effect=mock_xtb),
            patch("deltahf.pipeline.check_connectivity", return_value=True),
        ):
            result = process_molecule(
                "C", n_conformers=5, work_dir=tmp_path / "work",
                use_xtb_wbos=True, cache=cache,
            )

        # Should have run xTB (not used cache) since use_xtb_wbos is set
        assert result.xtb_energy == -10.0


class TestCacheCLI:
    def test_fit_accepts_cache_dir(self):
        from deltahf.__main__ import build_parser

        args = build_parser().parse_args(["fit", "-i", "data.csv", "--cache-dir", "/tmp/cache"])
        assert args.cache_dir == "/tmp/cache"

    def test_predict_accepts_cache_dir(self):
        from deltahf.__main__ import build_parser

        args = build_parser().parse_args(
            ["predict", "-i", "mol.csv", "--epsilon", "e.json", "--cache-dir", "/tmp/cache"]
        )
        assert args.cache_dir == "/tmp/cache"

    def test_cache_dir_default_is_none(self):
        from deltahf.__main__ import build_parser

        args = build_parser().parse_args(["fit", "-i", "data.csv"])
        assert args.cache_dir is None
