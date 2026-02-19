"""Tests for element validation."""

import pytest

from deltahf.pipeline import process_molecule
from deltahf.smiles import SUPPORTED_ELEMENTS, validate_elements


class TestElementValidation:
    """Test that supported elements are accepted and unsupported elements are rejected."""

    def test_validate_elements_chno(self):
        """Valid CHNO molecules should pass."""
        validate_elements("C")  # methane
        validate_elements("CCO")  # ethanol
        validate_elements("c1ccccc1")  # benzene
        validate_elements("C[N+](=O)[O-]")  # nitromethane

    def test_validate_elements_fluorine(self):
        """Molecules with fluorine should now be accepted."""
        validate_elements("CF")  # fluoromethane
        validate_elements("FC(F)(F)C(F)(F)F")  # hexafluoroethane

    def test_validate_elements_chlorine(self):
        """Molecules with chlorine should now be accepted."""
        validate_elements("CCCl")  # chloropropane

    def test_validate_elements_sulfur(self):
        """Molecules with sulfur should now be accepted."""
        validate_elements("CSC")  # dimethyl sulfide

    def test_validate_elements_fscl_combined(self):
        """Molecules containing F, S, and Cl together should be accepted."""
        validate_elements("ClCS")  # chloromethyl sulfide (simplified)
        validate_elements("FS(F)(F)(F)(F)F")  # SF6

    def test_validate_elements_unsupported(self):
        """Molecules with truly unsupported elements (e.g. Br, P, Si) should be rejected."""
        with pytest.raises(ValueError, match="unsupported elements"):
            validate_elements("CBr")  # bromomethane

        with pytest.raises(ValueError, match="unsupported elements"):
            validate_elements("CP")  # methylphosphine

        with pytest.raises(ValueError, match="unsupported elements"):
            validate_elements("C[Si](C)(C)C")  # tetramethylsilane

    def test_process_molecule_unsupported_elements(self):
        """Pipeline should catch unsupported elements and set error."""
        result = process_molecule("CBr", n_conformers=1)
        assert result.error is not None
        assert "unsupported elements" in result.error.lower()
        assert "Br" in result.error

    def test_process_molecule_supported_elements(self):
        """Pipeline should accept valid CHNO molecules."""
        result = process_molecule("CCO", n_conformers=1)
        # Error might be set for other reasons (xTB not available in tests),
        # but should not be about unsupported elements
        if result.error:
            assert "unsupported elements" not in result.error.lower()

    def test_supported_elements_constant(self):
        """Verify SUPPORTED_ELEMENTS constant includes F, S, Cl."""
        assert SUPPORTED_ELEMENTS == {"C", "H", "N", "O", "F", "S", "Cl"}
