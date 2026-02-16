"""Tests for SMILES parsing, atom counting, and bond order classification."""

import pytest

from deltahf.smiles import classify_atoms_7param, count_atoms, smiles_to_mol


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


class TestCountAtoms:
    def test_methane(self):
        counts = count_atoms("C")
        assert counts == {"C": 1, "H": 4, "N": 0, "O": 0}

    def test_water(self):
        counts = count_atoms("O")
        assert counts == {"C": 0, "H": 2, "N": 0, "O": 1}

    def test_ammonia(self):
        counts = count_atoms("N")
        assert counts == {"C": 0, "H": 3, "N": 1, "O": 0}

    def test_nitromethane(self):
        counts = count_atoms("C[N+](=O)[O-]")
        assert counts == {"C": 1, "H": 3, "N": 1, "O": 2}

    def test_rdx(self):
        counts = count_atoms("O=[N+]([O-])N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1")
        assert counts == {"C": 3, "H": 6, "N": 6, "O": 6}

    def test_benzene(self):
        counts = count_atoms("c1ccccc1")
        assert counts == {"C": 6, "H": 6, "N": 0, "O": 0}

    def test_ethanol(self):
        counts = count_atoms("CCO")
        assert counts == {"C": 2, "H": 6, "N": 0, "O": 1}

    def test_hydrogen_cyanide(self):
        counts = count_atoms("C#N")
        assert counts == {"C": 1, "H": 1, "N": 1, "O": 0}

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
        """Result dict always has all 7 keys."""
        counts = classify_atoms_7param("C")
        expected_keys = {"C", "H", "N", "O", "C_prime", "N_prime", "O_prime"}
        assert set(counts.keys()) == expected_keys
