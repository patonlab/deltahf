"""Tests for the end-to-end pipeline."""

from unittest import mock

import pytest

from deltahf.pipeline import MoleculeResult, process_molecule
from deltahf.xtb import XtbResult


class TestProcessMolecule:
    def test_returns_molecule_result(self):
        """Pipeline returns a MoleculeResult dataclass."""
        mock_result = XtbResult(energy=-5.07, optimized_xyz_path=None, converged=True, stdout="")
        with mock.patch("deltahf.pipeline.run_xtb_optimization", return_value=mock_result):
            result = process_molecule("C", n_conformers=1)
        assert isinstance(result, MoleculeResult)
        assert result.smiles == "C"

    def test_records_atom_counts(self):
        mock_result = XtbResult(energy=-5.07, optimized_xyz_path=None, converged=True, stdout="")
        with mock.patch("deltahf.pipeline.run_xtb_optimization", return_value=mock_result):
            result = process_molecule("C", n_conformers=1)
        assert {k: result.atom_counts_element[k] for k in ("C", "H", "N", "O")} == {"C": 1, "H": 4, "N": 0, "O": 0}
        assert result.atom_counts_element_bo is not None

    def test_records_xtb_energy(self):
        mock_result = XtbResult(energy=-5.07, optimized_xyz_path=None, converged=True, stdout="")
        with mock.patch("deltahf.pipeline.run_xtb_optimization", return_value=mock_result):
            result = process_molecule("C", n_conformers=1)
        assert result.xtb_energy == pytest.approx(-5.07)

    def test_invalid_smiles_records_error(self):
        result = process_molecule("not_valid", n_conformers=1)
        assert result.error is not None
        assert result.xtb_energy is None

    def test_failed_xtb_records_error(self):
        mock_result = XtbResult(energy=float("nan"), optimized_xyz_path=None, converged=False, stdout="failed")
        with mock.patch("deltahf.pipeline.run_xtb_optimization", return_value=mock_result):
            result = process_molecule("C", n_conformers=1)
        assert result.error is not None

    def test_predicts_dhf_with_epsilon(self):
        epsilon_4 = {"C": -200.0, "H": -80.0, "N": -190.0, "O": -95.0}
        mock_result = XtbResult(energy=-5.07, optimized_xyz_path=None, converged=True, stdout="")
        with mock.patch("deltahf.pipeline.run_xtb_optimization", return_value=mock_result):
            result = process_molecule("C", n_conformers=1, epsilon_element=epsilon_4)
        assert result.dhf_element is not None

    def test_name_is_stored(self):
        mock_result = XtbResult(energy=-5.07, optimized_xyz_path=None, converged=True, stdout="")
        with mock.patch("deltahf.pipeline.run_xtb_optimization", return_value=mock_result):
            result = process_molecule("C", n_conformers=1, name="methane")
        assert result.name == "methane"


class TestProcessMoleculeIntegration:
    @pytest.mark.integration
    def test_full_pipeline_methane(self, tmp_path, xtb_available):
        result = process_molecule("C", n_conformers=1, work_dir=tmp_path)
        assert result.error is None
        assert result.xtb_energy is not None
        assert result.xtb_energy < 0

    @pytest.mark.integration
    def test_full_pipeline_water(self, tmp_path, xtb_available):
        result = process_molecule("O", n_conformers=1, work_dir=tmp_path)
        assert result.error is None
        assert result.xtb_energy is not None
