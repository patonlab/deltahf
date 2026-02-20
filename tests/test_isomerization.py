"""Tests for isomerization detection via connectivity checking."""

from pathlib import Path
from unittest.mock import patch

from rdkit import Chem
from rdkit.Chem import rdDistGeom

from deltahf.conformers import check_connectivity, write_xyz
from deltahf.smiles import smiles_to_mol


class TestCheckConnectivity:
    def _make_xyz(self, smiles: str, tmp_path: Path) -> tuple[Chem.Mol, Path]:
        """Helper: embed a molecule and write its XYZ."""
        mol = smiles_to_mol(smiles)
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 42
        rdDistGeom.EmbedMolecule(mol, params)
        xyz_path = tmp_path / "mol.xyz"
        write_xyz(mol, 0, xyz_path)
        return mol, xyz_path

    def test_same_connectivity_returns_true(self, tmp_path):
        """Round-trip: write XYZ, read back, connectivity should match."""
        mol, xyz_path = self._make_xyz("CCO", tmp_path)
        assert check_connectivity(mol, xyz_path) is True

    def test_methane_round_trip(self, tmp_path):
        mol, xyz_path = self._make_xyz("C", tmp_path)
        assert check_connectivity(mol, xyz_path) is True

    def test_nitromethane_round_trip(self, tmp_path):
        mol, xyz_path = self._make_xyz("C[N+](=O)[O-]", tmp_path)
        assert check_connectivity(mol, xyz_path) is True

    def test_different_molecule_returns_false(self, tmp_path):
        """Pass ethanol mol but propane XYZ — different connectivity."""
        mol_ethanol = smiles_to_mol("CCO")
        mol_propane = smiles_to_mol("CCC")
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 42
        rdDistGeom.EmbedMolecule(mol_propane, params)
        xyz_path = tmp_path / "propane.xyz"
        write_xyz(mol_propane, 0, xyz_path)
        # mol_ethanol has different connectivity than propane XYZ
        assert check_connectivity(mol_ethanol, xyz_path) is False

    def test_broken_bond_detected(self, tmp_path):
        """Move a heavy atom far away to simulate a bond break."""
        mol, xyz_path = self._make_xyz("CCO", tmp_path)
        # XYZ atom order: C, C, O, then Hs. Move the O (line index 4) far away.
        lines = xyz_path.read_text().splitlines()
        parts = lines[4].split()
        lines[4] = f"{parts[0]}  100.000000   100.000000   100.000000"
        xyz_path.write_text("\n".join(lines) + "\n")
        assert check_connectivity(mol, xyz_path) is False

    def test_proton_transfer_allowed(self, tmp_path):
        """H-only connectivity change (proton transfer / tautomer) should pass."""
        mol, xyz_path = self._make_xyz("CCO", tmp_path)
        # Move the last H far away — heavy-atom skeleton is unchanged
        lines = xyz_path.read_text().splitlines()
        parts = lines[-1].split()
        lines[-1] = f"{parts[0]}  100.000000   100.000000   100.000000"
        xyz_path.write_text("\n".join(lines) + "\n")
        assert check_connectivity(mol, xyz_path) is True

    def test_nonexistent_file_returns_false(self, tmp_path):
        mol = smiles_to_mol("C")
        assert check_connectivity(mol, tmp_path / "nonexistent.xyz") is False

    def test_empty_file_returns_false(self, tmp_path):
        mol = smiles_to_mol("C")
        empty_file = tmp_path / "empty.xyz"
        empty_file.write_text("")
        assert check_connectivity(mol, empty_file) is False

    def test_atom_count_mismatch_returns_false(self, tmp_path):
        """Mol has 9 atoms (ethanol+Hs), XYZ has 5 atoms (methane+Hs)."""
        mol_ethanol = smiles_to_mol("CCO")
        mol_methane = smiles_to_mol("C")
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 42
        rdDistGeom.EmbedMolecule(mol_methane, params)
        xyz_path = tmp_path / "methane.xyz"
        write_xyz(mol_methane, 0, xyz_path)
        assert check_connectivity(mol_ethanol, xyz_path) is False

    def test_h2_round_trip(self, tmp_path):
        """H2 round-trip should pass connectivity check."""
        mol, xyz_path = self._make_xyz("[HH]", tmp_path)
        assert check_connectivity(mol, xyz_path) is True

    def test_h2_stretched_bond_passes(self, tmp_path):
        """H2 with a slightly stretched bond should still pass (VdW fallback)."""
        mol = smiles_to_mol("[HH]")
        # Write XYZ with H-H distance of 1.1 A (beyond CTD threshold, within VdW fallback)
        xyz_path = tmp_path / "h2_stretched.xyz"
        xyz_path.write_text("2\n\nH  0.000000  0.000000  0.000000\nH  0.000000  0.000000  1.100000\n")
        assert check_connectivity(mol, xyz_path) is True

    def test_h2_genuinely_dissociated_fails(self, tmp_path):
        """H2 with atoms 5 A apart should fail (genuine dissociation)."""
        mol = smiles_to_mol("[HH]")
        xyz_path = tmp_path / "h2_dissociated.xyz"
        xyz_path.write_text("2\n\nH  0.000000  0.000000  0.000000\nH  0.000000  0.000000  5.000000\n")
        assert check_connectivity(mol, xyz_path) is False


class TestIsomerizationInPipeline:
    def test_isomerized_conformer_skipped(self, tmp_path):
        """Mock check_connectivity to return False for first conformer, True for second."""
        from deltahf.pipeline import process_molecule
        from deltahf.xtb import XtbResult

        call_count = 0

        def mock_check(mol, xyz_path):
            nonlocal call_count
            call_count += 1
            return call_count > 1  # First conformer isomerizes, second is fine

        def mock_xtb(xyz_path, **kwargs):
            return XtbResult(
                energy=-10.0,
                optimized_xyz_path=xyz_path.parent / "xtbopt.xyz",
                converged=True,
                stdout="",
                wbo_path=None,
            )

        with (
            patch("deltahf.pipeline.run_xtb_optimization", side_effect=mock_xtb),
            patch("deltahf.pipeline.check_connectivity", side_effect=mock_check),
        ):
            result = process_molecule("C", n_conformers=2, work_dir=tmp_path)

        assert result.n_conformers_isomerized == 1
        assert result.n_conformers_optimized == 2
        assert result.xtb_energy is not None
        assert result.error is None

    def test_all_conformers_isomerized_error(self, tmp_path):
        """When all conformers isomerize, error should mention isomerization."""
        from deltahf.pipeline import process_molecule
        from deltahf.xtb import XtbResult

        def mock_xtb(xyz_path, **kwargs):
            return XtbResult(
                energy=-10.0,
                optimized_xyz_path=xyz_path.parent / "xtbopt.xyz",
                converged=True,
                stdout="",
                wbo_path=None,
            )

        with (
            patch("deltahf.pipeline.run_xtb_optimization", side_effect=mock_xtb),
            patch("deltahf.pipeline.check_connectivity", return_value=False),
        ):
            result = process_molecule("C", n_conformers=3, work_dir=tmp_path)

        assert result.error is not None
        assert "isomerized" in result.error
        assert result.n_conformers_isomerized >= 1
