"""Tests for training data integrity."""

import re

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from deltahf.data import load_training_data


def _parse_formula(formula: str) -> dict[str, int]:
    """Parse molecular formula string into element counts."""
    counts: dict[str, int] = {}
    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if element:
            counts[element] = counts.get(element, 0) + (int(count) if count else 1)
    return counts


class TestTrainingData:
    @pytest.fixture
    def df(self):
        return load_training_data()

    def test_has_292_rows(self, df):
        assert len(df) == 292

    def test_required_columns(self, df):
        for col in ["id", "name", "formula", "smiles", "exp_dhf_kcal_mol", "source"]:
            assert col in df.columns

    def test_ids_are_sequential(self, df):
        assert list(df["id"]) == list(range(1, 293))

    def test_all_smiles_parseable(self, df):
        for idx, row in df.iterrows():
            mol = Chem.MolFromSmiles(row["smiles"])
            assert mol is not None, f"ID {row['id']} ({row['name']}): failed to parse SMILES {row['smiles']}"

    def test_nitromethane_entry(self, df):
        row = df[df["id"] == 1].iloc[0]
        assert row["name"] == "nitromethane"
        assert row["exp_dhf_kcal_mol"] == pytest.approx(-19.3)

    def test_methane_entry(self, df):
        row = df[df["id"] == 46].iloc[0]
        assert row["name"] == "methane"
        assert row["exp_dhf_kcal_mol"] == pytest.approx(-17.9)

    def test_water_entry(self, df):
        row = df[df["id"] == 48].iloc[0]
        assert row["name"] == "water"
        assert row["exp_dhf_kcal_mol"] == pytest.approx(-57.8)

    def test_tnt_entry(self, df):
        row = df[df["id"] == 26].iloc[0]
        assert row["name"] == "TNT"
        assert row["exp_dhf_kcal_mol"] == pytest.approx(5.75)

    def test_cyclohexane_entry(self, df):
        row = df[df["name"] == "cyclohexane"].iloc[0]
        assert row["exp_dhf_kcal_mol"] == pytest.approx(-123.4 / 4.184, abs=0.01)
        assert row["source"] == "Yalamanchi2020"

    def test_sources(self, df):
        sources = set(df["source"])
        assert sources == {"Cawkwell2021", "Yalamanchi2020"}
        assert len(df[df["source"] == "Cawkwell2021"]) == 103
        assert len(df[df["source"] == "Yalamanchi2020"]) == 189

    def test_exp_dhf_values_reasonable(self, df):
        """Experimental ΔHf° should be in a reasonable range for organic molecules."""
        assert df["exp_dhf_kcal_mol"].min() > -200
        assert df["exp_dhf_kcal_mol"].max() < 200

    def test_no_duplicate_ids(self, df):
        assert df["id"].is_unique

    def test_formulas_match_smiles(self, df):
        """Atom counts from RDKit should match the stated formula."""
        mismatches = []
        for _, row in df.iterrows():
            mol = Chem.MolFromSmiles(row["smiles"])
            if mol is None:
                continue
            mol = Chem.AddHs(mol)
            rdkit_formula = rdMolDescriptors.CalcMolFormula(mol)
            if _parse_formula(rdkit_formula) != _parse_formula(row["formula"]):
                mismatches.append(
                    f"ID {row['id']} ({row['name']}): stated {row['formula']}, "
                    f"RDKit gives {rdkit_formula} from {row['smiles']}"
                )
        assert not mismatches, "Formula mismatches:\n" + "\n".join(mismatches)
