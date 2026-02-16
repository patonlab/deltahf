"""Tests for atom equivalent energy fitting, prediction, and cross-validation."""

import pytest

from deltahf.atom_equivalents import (
    PARAM_NAMES_4,
    fit_atom_equivalents,
    kfold_cross_validation,
    max_abs_deviation,
    predict_dhf,
    rmsd,
)


class TestPredictDhf:
    def test_simple_prediction(self):
        """ΔHf = u - Σ(nl × εl)."""
        epsilon = {"C": -200.0, "H": -80.0, "N": -190.0, "O": -95.0}
        atom_counts = {"C": 1, "H": 4, "N": 0, "O": 0}
        u_kcal = -520.0
        # -520 - (1*-200 + 4*-80 + 0 + 0) = -520 - (-520) = 0
        assert predict_dhf(u_kcal, atom_counts, epsilon) == pytest.approx(0.0)

    def test_with_all_atom_types(self):
        epsilon = {"C": -200.0, "H": -80.0, "N": -190.0, "O": -95.0}
        atom_counts = {"C": 1, "H": 3, "N": 1, "O": 2}
        u_kcal = -1000.0
        expected = -1000.0 - (1 * -200.0 + 3 * -80.0 + 1 * -190.0 + 2 * -95.0)
        assert predict_dhf(u_kcal, atom_counts, epsilon) == pytest.approx(expected)

    def test_7param_prediction(self):
        epsilon = {
            "C": -200.0, "H": -80.0, "N": -190.0, "O": -95.0,
            "C_prime": -205.0, "N_prime": -185.0, "O_prime": -93.0,
        }
        atom_counts = {
            "C": 0, "H": 3, "N": 0, "O": 0,
            "C_prime": 1, "N_prime": 1, "O_prime": 2,
        }
        u_kcal = -818.0
        expected = -818.0 - (3 * -80.0 + 1 * -205.0 + 1 * -185.0 + 2 * -93.0)
        assert predict_dhf(u_kcal, atom_counts, epsilon) == pytest.approx(expected)

    def test_missing_atom_type_treated_as_zero(self):
        epsilon = {"C": -200.0, "H": -80.0, "N": -190.0, "O": -95.0}
        atom_counts = {"C": 1, "H": 4}  # missing N and O
        u_kcal = -520.0
        expected = -520.0 - (1 * -200.0 + 4 * -80.0)
        assert predict_dhf(u_kcal, atom_counts, epsilon) == pytest.approx(expected)


class TestFitAtomEquivalents:
    def test_exact_fit_single_molecule(self):
        """Fitting to a single molecule should reproduce its ΔHf exactly."""
        atom_counts_list = [{"C": 1, "H": 4, "N": 0, "O": 0}]
        u_values = [-2642.0]
        exp_dhf = [-17.9]

        epsilon = fit_atom_equivalents(atom_counts_list, u_values, exp_dhf, PARAM_NAMES_4)
        predicted = predict_dhf(u_values[0], atom_counts_list[0], epsilon)
        assert predicted == pytest.approx(-17.9, abs=0.1)

    def test_recovers_known_epsilon(self):
        """With synthetic data, should recover the true epsilon values."""
        true_epsilon = {"C": -200.0, "H": -80.0, "N": -190.0, "O": -95.0}
        molecules = [
            {"C": 1, "H": 4, "N": 0, "O": 0},
            {"C": 0, "H": 2, "N": 0, "O": 1},
            {"C": 1, "H": 3, "N": 1, "O": 2},
            {"C": 2, "H": 6, "N": 0, "O": 1},
            {"C": 6, "H": 6, "N": 0, "O": 0},
        ]
        exp_dhf = [-17.9, -57.8, -19.3, -56.2, 19.7]
        u_values = [
            exp_dhf[i] + sum(molecules[i].get(k, 0) * true_epsilon[k] for k in true_epsilon)
            for i in range(len(molecules))
        ]

        fitted = fit_atom_equivalents(molecules, u_values, exp_dhf, PARAM_NAMES_4)
        for key in true_epsilon:
            assert fitted[key] == pytest.approx(true_epsilon[key], abs=0.1)

    def test_returns_dict_with_param_names(self):
        atom_counts_list = [
            {"C": 1, "H": 4, "N": 0, "O": 0},
            {"C": 0, "H": 2, "N": 0, "O": 1},
        ]
        u_values = [-2642.0, -450.0]
        exp_dhf = [-17.9, -57.8]
        epsilon = fit_atom_equivalents(atom_counts_list, u_values, exp_dhf, PARAM_NAMES_4)
        assert set(epsilon.keys()) == set(PARAM_NAMES_4)


class TestKfoldCrossValidation:
    @pytest.fixture
    def synthetic_data(self):
        """Generate synthetic training data with known epsilon."""
        true_epsilon = {"C": -200.0, "H": -80.0, "N": -190.0, "O": -95.0}
        molecules = [
            {"C": 1, "H": 4, "N": 0, "O": 0},
            {"C": 0, "H": 2, "N": 0, "O": 1},
            {"C": 1, "H": 3, "N": 1, "O": 2},
            {"C": 2, "H": 6, "N": 0, "O": 1},
            {"C": 6, "H": 6, "N": 0, "O": 0},
            {"C": 0, "H": 3, "N": 1, "O": 0},
            {"C": 2, "H": 6, "N": 0, "O": 0},
            {"C": 3, "H": 8, "N": 0, "O": 0},
            {"C": 1, "H": 1, "N": 1, "O": 0},
            {"C": 0, "H": 0, "N": 2, "O": 0},
        ]
        exp_dhf = [-17.9, -57.8, -19.3, -56.2, 19.7, -5.5, -20.1, -25.0, 31.5, 0.0]
        u_values = [
            exp_dhf[i] + sum(molecules[i].get(k, 0) * true_epsilon[k] for k in true_epsilon)
            for i in range(len(molecules))
        ]
        return molecules, u_values, exp_dhf

    def test_returns_correct_structure(self, synthetic_data):
        molecules, u_values, exp_dhf = synthetic_data
        result = kfold_cross_validation(molecules, u_values, exp_dhf, PARAM_NAMES_4, k=5)
        assert "cv_error" in result
        assert "mean_epsilon" in result
        assert "fold_results" in result
        assert len(result["fold_results"]) == 5

    def test_cv_error_is_non_negative(self, synthetic_data):
        molecules, u_values, exp_dhf = synthetic_data
        result = kfold_cross_validation(molecules, u_values, exp_dhf, PARAM_NAMES_4, k=5)
        assert result["cv_error"] >= 0

    def test_mean_epsilon_has_all_params(self, synthetic_data):
        molecules, u_values, exp_dhf = synthetic_data
        result = kfold_cross_validation(molecules, u_values, exp_dhf, PARAM_NAMES_4, k=5)
        assert set(result["mean_epsilon"].keys()) == set(PARAM_NAMES_4)

    def test_perfect_data_gives_small_cv_error(self, synthetic_data):
        """With perfectly consistent synthetic data, CV error should be very small."""
        molecules, u_values, exp_dhf = synthetic_data
        result = kfold_cross_validation(molecules, u_values, exp_dhf, PARAM_NAMES_4, k=5)
        assert result["cv_error"] < 1.0  # should be near zero for perfect data


class TestStatistics:
    def test_rmsd(self):
        predicted = [1.0, 2.0, 3.0]
        experimental = [1.1, 1.9, 3.2]
        # sqrt(mean([0.01, 0.01, 0.04])) = sqrt(0.02) ≈ 0.1414
        assert rmsd(predicted, experimental) == pytest.approx(0.1414, abs=0.001)

    def test_rmsd_perfect(self):
        assert rmsd([1.0, 2.0], [1.0, 2.0]) == pytest.approx(0.0)

    def test_max_abs_deviation(self):
        predicted = [1.0, 2.0, 3.0]
        experimental = [1.1, 1.9, 3.5]
        assert max_abs_deviation(predicted, experimental) == pytest.approx(0.5)

    def test_max_abs_deviation_symmetric(self):
        assert max_abs_deviation([5.0], [3.0]) == pytest.approx(2.0)
        assert max_abs_deviation([3.0], [5.0]) == pytest.approx(2.0)
