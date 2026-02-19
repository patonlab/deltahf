"""Tests for SMILES parsing, atom counting, and bond order classification."""

import pytest

from deltahf.smiles import (
    classify_atoms_7param,
    classify_atoms_7param_from_wbo,
    classify_atoms_extended,
    classify_atoms_hybrid,
    count_atoms,
    smiles_to_mol,
)


class TestSmilesToMol:
    def test_valid_smiles(self):
        mol = smiles_to_mol("C")
        assert mol is not None
        assert mol.GetNumAtoms() == 5  # CH4 with explicit Hs

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            smiles_to_mol("not_a_molecule")

    def test_empty_smiles_raises(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            smiles_to_mol("")

    def test_adds_hydrogens(self):
        mol = smiles_to_mol("O")
        h_count = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")
        assert h_count == 2


_CHNO = ("C", "H", "N", "O")


class TestCountAtoms:
    def test_methane(self):
        counts = count_atoms("C")
        assert {k: counts[k] for k in _CHNO} == {"C": 1, "H": 4, "N": 0, "O": 0}

    def test_water(self):
        counts = count_atoms("O")
        assert {k: counts[k] for k in _CHNO} == {"C": 0, "H": 2, "N": 0, "O": 1}

    def test_ammonia(self):
        counts = count_atoms("N")
        assert {k: counts[k] for k in _CHNO} == {"C": 0, "H": 3, "N": 1, "O": 0}

    def test_nitromethane(self):
        counts = count_atoms("C[N+](=O)[O-]")
        assert {k: counts[k] for k in _CHNO} == {"C": 1, "H": 3, "N": 1, "O": 2}

    def test_rdx(self):
        counts = count_atoms("O=[N+]([O-])N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1")
        assert {k: counts[k] for k in _CHNO} == {"C": 3, "H": 6, "N": 6, "O": 6}

    def test_benzene(self):
        counts = count_atoms("c1ccccc1")
        assert {k: counts[k] for k in _CHNO} == {"C": 6, "H": 6, "N": 0, "O": 0}

    def test_ethanol(self):
        counts = count_atoms("CCO")
        assert {k: counts[k] for k in _CHNO} == {"C": 2, "H": 6, "N": 0, "O": 1}

    def test_hydrogen_cyanide(self):
        counts = count_atoms("C#N")
        assert {k: counts[k] for k in _CHNO} == {"C": 1, "H": 1, "N": 1, "O": 0}

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError):
            count_atoms("not_valid")


class TestClassifyAtoms7Param:
    def test_benzene_all_c_prime(self):
        """All 6 carbons in benzene are C' (aromatic bonds = 1.5 > 1.25)."""
        counts = classify_atoms_7param("c1ccccc1")
        assert counts["C_prime"] == 6
        assert counts["C"] == 0
        assert counts["H"] == 6

    def test_ethane_all_c_single(self):
        """Both carbons in ethane are C (only single bonds)."""
        counts = classify_atoms_7param("CC")
        assert counts["C"] == 2
        assert counts["C_prime"] == 0

    def test_ethene_all_c_prime(self):
        """Both carbons in ethene are C' (double bond > 1.25)."""
        counts = classify_atoms_7param("C=C")
        assert counts["C_prime"] == 2
        assert counts["C"] == 0

    def test_ethyne_all_c_prime(self):
        """Both carbons in ethyne are C' (triple bond > 1.25)."""
        counts = classify_atoms_7param("C#C")
        assert counts["C_prime"] == 2
        assert counts["C"] == 0

    def test_hydrogen_cyanide(self):
        """HCN: C#N -> both C' and N'."""
        counts = classify_atoms_7param("C#N")
        assert counts["C_prime"] == 1
        assert counts["N_prime"] == 1

    def test_methane_no_primed(self):
        """CH4: only single bonds, no primed atoms."""
        counts = classify_atoms_7param("C")
        assert counts["C"] == 1
        assert counts["C_prime"] == 0
        assert counts["H"] == 4

    def test_h_never_primed(self):
        """H atoms are never classified as primed."""
        counts = classify_atoms_7param("C=C")
        assert counts["H"] == 4
        assert "H_prime" not in counts or counts.get("H_prime", 0) == 0

    def test_nitromethane(self):
        """CH3NO2: C single-bonded -> C; N in NO2 has double bond -> N';
        O atoms in NO2 have resonance -> O'."""
        counts = classify_atoms_7param("C[N+](=O)[O-]")
        assert counts["C"] == 1
        assert counts["C_prime"] == 0
        assert counts["N_prime"] == 1
        assert counts["N"] == 0
        assert counts["H"] == 3

    def test_tnt(self):
        """TNT: 6 aromatic C -> C', 1 methyl C -> C, 3 nitro N -> N', 6 nitro O -> O'."""
        counts = classify_atoms_7param("Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]")
        assert counts["C_prime"] == 6  # aromatic ring
        assert counts["C"] == 1  # methyl
        assert counts["N_prime"] == 3  # three nitro groups
        assert counts["O_prime"] == 6  # six nitro oxygens

    def test_water(self):
        """H2O: O has only single bonds -> O (not O')."""
        counts = classify_atoms_7param("O")
        assert counts["O"] == 1
        assert counts["O_prime"] == 0

    def test_formaldehyde(self):
        """CH2O: C=O -> both C' and O'."""
        counts = classify_atoms_7param("C=O")
        assert counts["C_prime"] == 1
        assert counts["O_prime"] == 1
        assert counts["C"] == 0
        assert counts["O"] == 0

    def test_carbon_dioxide(self):
        """CO2: O=C=O -> C' and 2 O'."""
        counts = classify_atoms_7param("O=C=O")
        assert counts["C_prime"] == 1
        assert counts["O_prime"] == 2

    def test_returns_all_seven_keys(self):
        """Result dict always has all keys (CHNO + FSCl variants)."""
        counts = classify_atoms_7param("C")
        expected_keys = {"C", "H", "N", "O", "C_prime", "N_prime", "O_prime", "F", "S", "S_prime", "Cl"}
        assert set(counts.keys()) == expected_keys


class TestClassifyAtoms7ParamFromWbo:
    """Test 7-param classification using xTB Wiberg bond orders."""

    def test_nitromethane_with_xtb_wbos(self):
        """Nitromethane WBOs: N-O bonds ~1.48 (> 1.25 -> N' and O')."""
        # Atom order for C[N+](=O)[O-] with explicit H: C(0), N(1), O(2), O(3), H(4), H(5), H(6)
        wbos = {
            (0, 1): 0.944, (1, 0): 0.944,
            (1, 2): 1.484, (2, 1): 1.484,
            (1, 3): 1.484, (3, 1): 1.484,
            (2, 3): 0.377, (3, 2): 0.377,
            (0, 4): 0.974, (4, 0): 0.974,
            (0, 5): 0.974, (5, 0): 0.974,
            (0, 6): 0.951, (6, 0): 0.951,
        }
        counts = classify_atoms_7param_from_wbo("C[N+](=O)[O-]", wbos)
        assert counts["C"] == 1  # C-N bond is 0.944 < 1.25
        assert counts["N_prime"] == 1  # N-O bonds are 1.484 > 1.25
        assert counts["O_prime"] == 2  # O-N bonds are 1.484 > 1.25
        assert counts["H"] == 3

    def test_methane_all_single(self):
        """Methane: all C-H bonds ~0.99 (< 1.25 -> C not C')."""
        wbos = {
            (0, 1): 0.99, (1, 0): 0.99,
            (0, 2): 0.99, (2, 0): 0.99,
            (0, 3): 0.99, (3, 0): 0.99,
            (0, 4): 0.99, (4, 0): 0.99,
        }
        counts = classify_atoms_7param_from_wbo("C", wbos)
        assert counts["C"] == 1
        assert counts["C_prime"] == 0
        assert counts["H"] == 4

    def test_returns_all_seven_keys(self):
        """Result dict always has all keys (CHNO + FSCl variants)."""
        wbos = {(0, 1): 0.99, (1, 0): 0.99}
        counts = classify_atoms_7param_from_wbo("C", wbos)
        expected_keys = {"C", "H", "N", "O", "C_prime", "N_prime", "O_prime", "F", "S", "S_prime", "Cl"}
        assert set(counts.keys()) == expected_keys


class TestClassifyAtomsHybrid:
    """Tests for 9-parameter hybridization-based classification."""

    def test_methane(self):
        counts = classify_atoms_hybrid("C")
        assert counts["C_sp3"] == 1
        assert counts["H"] == 4

    def test_ethane(self):
        counts = classify_atoms_hybrid("CC")
        assert counts["C_sp3"] == 2
        assert counts["C_sp2"] == 0
        assert counts["H"] == 6

    def test_ethene(self):
        counts = classify_atoms_hybrid("C=C")
        assert counts["C_sp2"] == 2
        assert counts["C_sp3"] == 0
        assert counts["H"] == 4

    def test_ethyne(self):
        counts = classify_atoms_hybrid("C#C")
        assert counts["C_sp"] == 2
        assert counts["H"] == 2

    def test_benzene(self):
        counts = classify_atoms_hybrid("c1ccccc1")
        assert counts["C_sp2"] == 6
        assert counts["H"] == 6

    def test_hydrogen_cyanide(self):
        counts = classify_atoms_hybrid("C#N")
        assert counts["C_sp"] == 1
        assert counts["N_sp"] == 1
        assert counts["H"] == 1

    def test_water(self):
        counts = classify_atoms_hybrid("O")
        assert counts["O_sp3"] == 1
        assert counts["H"] == 2

    def test_formaldehyde(self):
        counts = classify_atoms_hybrid("C=O")
        assert counts["C_sp2"] == 1
        assert counts["O_sp2"] == 1
        assert counts["H"] == 2

    def test_ethanol(self):
        counts = classify_atoms_hybrid("CCO")
        assert counts["C_sp3"] == 2
        assert counts["O_sp3"] == 1
        assert counts["H"] == 6

    def test_ammonia(self):
        counts = classify_atoms_hybrid("N")
        assert counts["N_sp3"] == 1
        assert counts["H"] == 3

    def test_nitromethane(self):
        """CH3-NO2: methyl C is sp3, N is sp2, both O are sp2."""
        counts = classify_atoms_hybrid("C[N+](=O)[O-]")
        assert counts["C_sp3"] == 1
        assert counts["N_sp2"] == 1
        assert counts["O_sp2"] == 2
        assert counts["H"] == 3

    def test_tnt(self):
        """TNT: 6 aromatic sp2 C, 1 methyl sp3 C, 3 sp2 N, 6 sp2 O."""
        counts = classify_atoms_hybrid("Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]")
        assert counts["C_sp2"] == 6
        assert counts["C_sp3"] == 1
        assert counts["N_sp2"] == 3
        assert counts["O_sp2"] == 6

    def test_returns_all_ten_keys(self):
        counts = classify_atoms_hybrid("C")
        expected = {
            "C_sp3", "C_sp2", "C_sp", "H",
            "N_sp3", "N_sp2", "N_sp",
            "O_sp3", "O_sp2", "O_sp",
            "F_sp3", "S_sp3", "S_sp2", "S_sp", "Cl_sp3",
        }
        assert set(counts.keys()) == expected

    def test_total_atoms_match_4param(self):
        """Hybrid counts should sum to the same totals as 4-param for each element."""
        for smiles in ["C", "CCO", "c1ccccc1", "C#N", "C[N+](=O)[O-]"]:
            c4 = count_atoms(smiles)
            ch = classify_atoms_hybrid(smiles)
            assert ch["C_sp3"] + ch["C_sp2"] + ch["C_sp"] == c4["C"]
            assert ch["H"] == c4["H"]
            assert ch["N_sp3"] + ch["N_sp2"] + ch["N_sp"] == c4["N"]
            assert ch["O_sp3"] + ch["O_sp2"] + ch["O_sp"] == c4["O"]


class TestClassifyAtomsExtended:
    """Tests for 13-parameter extended classification (hybridization + H-count for C)."""

    def test_methane(self):
        """CH4: sp3 carbon with 4 H, capped to C_sp3_3H."""
        counts = classify_atoms_extended("C")
        assert counts["C_sp3_3H"] == 1
        assert counts["H"] == 4

    def test_ethane(self):
        """C2H6: two sp3 CH3 groups."""
        counts = classify_atoms_extended("CC")
        assert counts["C_sp3_3H"] == 2
        assert counts["H"] == 6

    def test_propane(self):
        """C3H8: two methyl sp3 + one methylene sp3."""
        counts = classify_atoms_extended("CCC")
        assert counts["C_sp3_3H"] == 2
        assert counts["C_sp3_2H"] == 1
        assert counts["H"] == 8

    def test_isobutane(self):
        """CH(CH3)3: three methyl + one methine sp3."""
        counts = classify_atoms_extended("CC(C)C")
        assert counts["C_sp3_3H"] == 3
        assert counts["C_sp3_1H"] == 1
        assert counts["H"] == 10

    def test_neopentane(self):
        """C(CH3)4: four methyl + one quaternary sp3."""
        counts = classify_atoms_extended("CC(C)(C)C")
        assert counts["C_sp3_3H"] == 4
        assert counts["C_sp3_0H"] == 1
        assert counts["H"] == 12

    def test_benzene(self):
        """C6H6: six sp2 CH."""
        counts = classify_atoms_extended("c1ccccc1")
        assert counts["C_sp2_1H"] == 6
        assert counts["H"] == 6

    def test_toluene(self):
        """Toluene: 5 sp2 CH, 1 sp2 C(no H), 1 sp3 CH3."""
        counts = classify_atoms_extended("Cc1ccccc1")
        assert counts["C_sp3_3H"] == 1
        assert counts["C_sp2_1H"] == 5
        assert counts["C_sp2_0H"] == 1

    def test_formaldehyde(self):
        """CH2O: sp2 C with 2H → C_sp2_2H."""
        counts = classify_atoms_extended("C=O")
        assert counts["C_sp2_2H"] == 1
        assert counts["O_sp2"] == 1
        assert counts["H"] == 2

    def test_acetylene(self):
        """HC≡CH: two sp carbons."""
        counts = classify_atoms_extended("C#C")
        assert counts["C_sp"] == 2
        assert counts["H"] == 2

    def test_water(self):
        counts = classify_atoms_extended("O")
        assert counts["O_sp3"] == 1
        assert counts["H"] == 2

    def test_returns_all_fifteen_keys(self):
        counts = classify_atoms_extended("C")
        expected = {
            "C_sp3_3H", "C_sp3_2H", "C_sp3_1H", "C_sp3_0H",
            "C_sp2_2H", "C_sp2_1H", "C_sp2_0H", "C_sp",
            "H", "N_sp3", "N_sp2", "N_sp", "O_sp3", "O_sp2", "O_sp",
            "F", "S_sp3", "S_sp2", "Cl",
        }
        assert set(counts.keys()) == expected

    def test_total_carbon_matches_4param(self):
        """Extended C counts should sum to the same total as 4-param."""
        for smiles in ["C", "CCC", "CC(C)C", "c1ccccc1", "Cc1ccccc1"]:
            c4 = count_atoms(smiles)
            ce = classify_atoms_extended(smiles)
            c_total = (
                ce["C_sp3_3H"] + ce["C_sp3_2H"] + ce["C_sp3_1H"] + ce["C_sp3_0H"]
                + ce["C_sp2_2H"] + ce["C_sp2_1H"] + ce["C_sp2_0H"] + ce["C_sp"]
            )
            assert c_total == c4["C"]
