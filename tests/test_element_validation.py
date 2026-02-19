"""Tests for element validation."""

import pytest

from deltahf.pipeline import process_molecule
from deltahf.smiles import SUPPORTED_ELEMENTS, validate_elements


class TestElementValidation:
    """Test that unsupported elements are properly rejected."""

    def test_validate_elements_chno_only(self):
        """Valid CHNO molecules should pass."""
        validate_elements("C")  # methane
        validate_elements("CCO")  # ethanol
        validate_elements("c1ccccc1")  # benzene
        validate_elements("C[N+](=O)[O-]")  # nitromethane

    def test_validate_elements_fluorine(self):
        """Molecules with fluorine should be rejected."""
        with pytest.raises(ValueError, match="unsupported elements.*\\['F'\\]"):
            validate_elements("CF")

        with pytest.raises(ValueError, match="unsupported elements.*\\['F'\\]"):
            validate_elements("FC(F)(F)C(F)(F)F")  # hexafluoroethane

    def test_validate_elements_chlorine(self):
        """Molecules with chlorine should be rejected."""
        with pytest.raises(ValueError, match="unsupported elements.*\\['Cl'\\]"):
            validate_elements("CCCl")

    def test_validate_elements_sulfur(self):
        """Molecules with sulfur should be rejected."""
        with pytest.raises(ValueError, match="unsupported elements.*\\['S'\\]"):
            validate_elements("CSC")  # dimethyl sulfide

    def test_validate_elements_multiple_unsupported(self):
        """Molecules with multiple unsupported elements should list all."""
        with pytest.raises(ValueError, match="unsupported elements"):
            validate_elements("FS(F)(F)(F)(F)F")  # SF6

    def test_process_molecule_unsupported_elements(self):
        """Pipeline should catch unsupported elements and set error."""
        result = process_molecule("FC(F)(F)F", n_conformers=1)
        assert result.error is not None
        assert "unsupported elements" in result.error.lower()
        assert "F" in result.error

    def test_process_molecule_supported_elements(self):
        """Pipeline should accept valid CHNO molecules."""
        result = process_molecule("CCO", n_conformers=1)
        # Error might be set for other reasons (xTB not available in tests),
        # but should not be about unsupported elements
        if result.error:
            assert "unsupported elements" not in result.error.lower()

    def test_supported_elements_constant(self):
        """Verify SUPPORTED_ELEMENTS constant."""
        assert SUPPORTED_ELEMENTS == {"C", "H", "N", "O"}
