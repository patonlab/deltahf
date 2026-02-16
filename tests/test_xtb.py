"""Tests for xTB subprocess wrapper."""

import pytest

from deltahf.xtb import build_xtb_command, parse_total_energy


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
