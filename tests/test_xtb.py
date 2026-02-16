"""Tests for xTB subprocess wrapper."""

import pytest

from deltahf.xtb import (
    build_xtb_command,
    find_gxtb_binary,
    parse_gxtb_energy_file,
    parse_total_energy,
    parse_wbo_file,
)


class TestParseTotalEnergy:
    def test_parse_from_summary_block(self):
        output = """
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         ::                     SUMMARY                    ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         :: total energy             -5.070322476938 Eh    ::
         :: gradient norm              0.019484395925 Eh/a  ::
         :: HOMO-LUMO gap             14.652302902752 eV    ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
        """
        energy = parse_total_energy(output)
        assert energy == pytest.approx(-5.070322476938, abs=1e-10)

    def test_parse_from_optimization_output(self):
        output = """
         :: total energy            -22.618809235400 Eh    ::
        """
        energy = parse_total_energy(output)
        assert energy == pytest.approx(-22.618809235400, abs=1e-10)

    def test_takes_last_match(self):
        """If multiple energy lines, take the last (post-optimization) value."""
        output = """
         :: total energy            -10.000000000000 Eh    ::
         :: total energy            -22.618809235400 Eh    ::
        """
        energy = parse_total_energy(output)
        assert energy == pytest.approx(-22.618809235400, abs=1e-10)

    def test_no_match_raises(self):
        with pytest.raises(RuntimeError, match="Could not parse"):
            parse_total_energy("no energy here")

    def test_empty_string_raises(self):
        with pytest.raises(RuntimeError):
            parse_total_energy("")


class TestBuildXtbCommand:
    def test_optimization_command(self):
        cmd = build_xtb_command("/tmp/mol.xyz", opt=True, gfn=2, charge=0)
        assert cmd == ["xtb", "/tmp/mol.xyz", "--opt", "--gfn", "2", "--chrg", "0"]

    def test_singlepoint_command(self):
        cmd = build_xtb_command("/tmp/mol.xyz", opt=False, gfn=2, charge=0)
        assert "xtb" == cmd[0]
        assert "/tmp/mol.xyz" in cmd
        assert "--opt" not in cmd

    def test_custom_gfn(self):
        cmd = build_xtb_command("/tmp/mol.xyz", gfn=1)
        assert "--gfn" in cmd
        idx = cmd.index("--gfn")
        assert cmd[idx + 1] == "1"

    def test_nonzero_charge(self):
        cmd = build_xtb_command("/tmp/mol.xyz", charge=-1)
        idx = cmd.index("--chrg")
        assert cmd[idx + 1] == "-1"

    def test_uhf(self):
        cmd = build_xtb_command("/tmp/mol.xyz", uhf=2)
        assert "--uhf" in cmd
        idx = cmd.index("--uhf")
        assert cmd[idx + 1] == "2"

    def test_uhf_zero_omitted(self):
        cmd = build_xtb_command("/tmp/mol.xyz", uhf=0)
        assert "--uhf" not in cmd


class TestParseWboFile:
    def test_parse_nitromethane_wbo(self, tmp_path):
        """Parse WBO file for nitromethane (xTB format: 1-based indices)."""
        wbo_content = (
            "           1           2  0.94413409148182448\n"
            "           2           3   1.4839822100747204\n"
            "           2           4   1.4839826063331538\n"
            "           3           4  0.37705420781274634\n"
            "           1           5  0.97358453582634685\n"
            "           1           6  0.97358614085771944\n"
            "           1           7  0.95072904780885048\n"
        )
        wbo_file = tmp_path / "wbo"
        wbo_file.write_text(wbo_content)

        wbos = parse_wbo_file(wbo_file)

        # Converted to 0-based indices
        assert wbos[(0, 1)] == pytest.approx(0.944, abs=0.001)
        assert wbos[(1, 0)] == pytest.approx(0.944, abs=0.001)  # symmetric
        assert wbos[(1, 2)] == pytest.approx(1.484, abs=0.001)  # N=O double bond
        assert wbos[(1, 3)] == pytest.approx(1.484, abs=0.001)

    def test_symmetric_lookup(self, tmp_path):
        wbo_file = tmp_path / "wbo"
        wbo_file.write_text("           1           2  1.5\n")
        wbos = parse_wbo_file(wbo_file)
        assert wbos[(0, 1)] == wbos[(1, 0)]


class TestFindGxtbBinary:
    def test_returns_none_if_not_found(self, monkeypatch):
        """gxtb is optional, so return None if not found."""
        monkeypatch.setattr("shutil.which", lambda x: None)
        assert find_gxtb_binary() is None

    def test_returns_path_if_found(self, monkeypatch):
        monkeypatch.setattr("shutil.which", lambda x: "/usr/bin/gxtb" if x == "gxtb" else None)
        assert find_gxtb_binary() == "/usr/bin/gxtb"


class TestParseGxtbEnergyFile:
    def test_parse_valid_energy_file(self, tmp_path):
        """Parse gxtb energy file format."""
        energy_file = tmp_path / "energy"
        energy_file.write_text("energy\n0.0 -5.123456789012\n$end\n")
        energy = parse_gxtb_energy_file(energy_file)
        assert energy == pytest.approx(-5.123456789012, abs=1e-10)

    def test_parse_with_different_first_column(self, tmp_path):
        """First column value doesn't matter, only second column."""
        energy_file = tmp_path / "energy"
        energy_file.write_text("energy\n42.0 -10.5\n$end\n")
        energy = parse_gxtb_energy_file(energy_file)
        assert energy == pytest.approx(-10.5, abs=1e-10)

    def test_insufficient_lines_raises(self, tmp_path):
        energy_file = tmp_path / "energy"
        energy_file.write_text("energy\n")
        with pytest.raises(RuntimeError, match="insufficient lines"):
            parse_gxtb_energy_file(energy_file)

    def test_insufficient_columns_raises(self, tmp_path):
        energy_file = tmp_path / "energy"
        energy_file.write_text("energy\n-5.0\n")
        with pytest.raises(RuntimeError, match="insufficient columns"):
            parse_gxtb_energy_file(energy_file)

    def test_empty_file_raises(self, tmp_path):
        energy_file = tmp_path / "energy"
        energy_file.write_text("")
        with pytest.raises(RuntimeError):
            parse_gxtb_energy_file(energy_file)
