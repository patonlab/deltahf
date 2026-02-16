"""Tests for 3D conformer generation."""

import pytest

from deltahf.conformers import generate_conformers, get_lowest_conformers, write_xyz


class TestGenerateConformers:
    def test_single_conformer(self):
        """Can generate at least 1 conformer for ethanol."""
        mol, energies = generate_conformers("CCO", num_confs=1)
        assert mol.GetNumConformers() >= 1
        assert len(energies) >= 1

    def test_multiple_conformers(self):
        """Request multiple conformers for butane."""
        mol, energies = generate_conformers("CCCC", num_confs=20)
        assert mol.GetNumConformers() >= 1

    def test_energies_sorted(self):
        """Returned energies are sorted ascending by energy."""
        mol, energies = generate_conformers("CCCC", num_confs=20)
        energy_vals = [e for e, _ in energies]
        assert energy_vals == sorted(energy_vals)

    def test_energies_are_tuples(self):
        """Each entry is (energy, conf_id) tuple."""
        mol, energies = generate_conformers("CCO", num_confs=5)
        for item in energies:
            assert len(item) == 2
            energy, cid = item
            assert isinstance(energy, float)
            assert isinstance(cid, int)

    def test_conformer_has_3d_coords(self):
        """Generated conformer has 3D coordinates."""
        mol, energies = generate_conformers("CCO", num_confs=1)
        conf = mol.GetConformer(energies[0][1])
        assert conf.Is3D()

    def test_invalid_smiles_raises(self):
        with pytest.raises((ValueError, RuntimeError)):
            generate_conformers("not_valid", num_confs=1)

    def test_simple_molecule(self):
        """Methane should embed successfully."""
        mol, energies = generate_conformers("C", num_confs=1)
        assert len(energies) >= 1


class TestGetLowestConformers:
    def test_returns_n_lowest(self):
        mol, energies = generate_conformers("CCCC", num_confs=20)
        lowest = get_lowest_conformers(mol, energies, n=3)
        assert len(lowest) <= 3

    def test_returns_fewer_if_not_enough(self):
        mol, energies = generate_conformers("C", num_confs=1)
        lowest = get_lowest_conformers(mol, energies, n=10)
        assert len(lowest) <= 10
        assert len(lowest) >= 1

    def test_returns_conf_ids(self):
        mol, energies = generate_conformers("CCO", num_confs=5)
        lowest = get_lowest_conformers(mol, energies, n=2)
        for cid in lowest:
            assert isinstance(cid, int)
            mol.GetConformer(cid)  # should not raise


class TestWriteXyz:
    def test_write_xyz_file(self, tmp_path):
        mol, energies = generate_conformers("CCO", num_confs=1)
        cid = energies[0][1]
        xyz_path = tmp_path / "test.xyz"
        write_xyz(mol, cid, xyz_path)
        assert xyz_path.exists()

    def test_xyz_format(self, tmp_path):
        mol, energies = generate_conformers("CCO", num_confs=1)
        cid = energies[0][1]
        xyz_path = tmp_path / "test.xyz"
        write_xyz(mol, cid, xyz_path)
        content = xyz_path.read_text()
        lines = content.strip().split("\n")
        natoms = int(lines[0])
        assert natoms == mol.GetNumAtoms()
        # Line 2 is comment, lines 3+ are atom coordinates
        assert len(lines) == natoms + 2

    def test_xyz_has_atom_symbols(self, tmp_path):
        mol, energies = generate_conformers("O", num_confs=1)
        cid = energies[0][1]
        xyz_path = tmp_path / "test.xyz"
        write_xyz(mol, cid, xyz_path)
        content = xyz_path.read_text()
        lines = content.strip().split("\n")
        symbols = [line.split()[0] for line in lines[2:]]
        assert "O" in symbols
        assert symbols.count("H") == 2
