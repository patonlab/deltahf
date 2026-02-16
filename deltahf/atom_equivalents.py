"""Atom equivalent energy fitting, prediction, and cross-validation."""

import numpy as np
from numpy.typing import NDArray

HARTREE_TO_KCAL = 627.5094740631

PARAM_NAMES_4 = ["C", "H", "N", "O"]
PARAM_NAMES_7 = ["C", "H", "N", "O", "C_prime", "N_prime", "O_prime"]
PARAM_NAMES_HYBRID = ["C_sp3", "C_sp2", "C_sp", "H", "N_sp3", "N_sp2", "N_sp", "O_sp3", "O_sp2", "O_sp"]
PARAM_NAMES_EXTENDED = [
    "C_sp3_3H", "C_sp3_2H", "C_sp3_1H", "C_sp3_0H",
    "C_sp2_2H", "C_sp2_1H", "C_sp2_0H", "C_sp",
    "H", "N_sp3", "N_sp2", "N_sp", "O_sp3", "O_sp2", "O_sp",
]


def predict_dhf(u_kcal: float, atom_counts: dict[str, int], epsilon: dict[str, float]) -> float:
    """Predict ΔHf° = u - Σ(nl × εl)."""
    correction = sum(atom_counts.get(k, 0) * epsilon[k] for k in epsilon)
    return u_kcal - correction


def build_design_matrix(atom_counts_list: list[dict[str, int]], param_names: list[str]) -> NDArray:
    """Build the N × p design matrix of atom counts."""
    n = len(atom_counts_list)
    p = len(param_names)
    matrix = np.zeros((n, p))
    for i, counts in enumerate(atom_counts_list):
        for j, name in enumerate(param_names):
            matrix[i, j] = counts.get(name, 0)
    return matrix


def fit_atom_equivalents(
    atom_counts_list: list[dict[str, int]],
    u_values_kcal: list[float],
    exp_dhf: list[float],
    param_names: list[str],
) -> dict[str, float]:
    """Fit atom equivalent energies by linear least squares.

    Solves: N @ ε = u - ΔHf_exp
    """
    matrix = build_design_matrix(atom_counts_list, param_names)
    u = np.array(u_values_kcal)
    dhf = np.array(exp_dhf)
    target = u - dhf

    epsilon, _, _, _ = np.linalg.lstsq(matrix, target, rcond=None)
    return {name: float(val) for name, val in zip(param_names, epsilon)}


def kfold_cross_validation(
    atom_counts_list: list[dict[str, int]],
    u_values_kcal: list[float],
    exp_dhf: list[float],
    param_names: list[str],
    k: int = 10,
    seed: int = 42,
) -> dict:
    """Perform k-fold cross-validation of atom equivalent energies.

    Returns dict with cv_error (average MSD), mean_epsilon, and fold_results.
    """
    rng = np.random.default_rng(seed)
    n = len(atom_counts_list)
    indices = rng.permutation(n)
    fold_size = n // k

    fold_results = []
    all_epsilon = []

    for fold in range(k):
        start = fold * fold_size
        if fold == k - 1:
            test_idx = indices[start:]
        else:
            test_idx = indices[start : start + fold_size]
        train_idx = np.setdiff1d(indices, test_idx)

        train_counts = [atom_counts_list[i] for i in train_idx]
        train_u = [u_values_kcal[i] for i in train_idx]
        train_dhf = [exp_dhf[i] for i in train_idx]

        epsilon = fit_atom_equivalents(train_counts, train_u, train_dhf, param_names)
        all_epsilon.append(epsilon)

        test_predictions = []
        test_experimental = []
        for i in test_idx:
            pred = predict_dhf(u_values_kcal[i], atom_counts_list[i], epsilon)
            test_predictions.append(pred)
            test_experimental.append(exp_dhf[i])

        msd = float(np.mean((np.array(test_predictions) - np.array(test_experimental)) ** 2))
        fold_results.append(
            {
                "epsilon": epsilon,
                "msd": msd,
                "predictions": test_predictions,
                "experimental": test_experimental,
            }
        )

    cv_error = float(np.mean([f["msd"] for f in fold_results]))
    cv_rmsd = float(np.sqrt(cv_error))
    mean_epsilon = {name: float(np.mean([f["epsilon"][name] for f in fold_results])) for name in param_names}
    std_epsilon = {
        name: float(np.std([f["epsilon"][name] for f in fold_results], ddof=1))
        for name in param_names
    }

    return {
        "cv_error": cv_error,
        "cv_rmsd": cv_rmsd,
        "mean_epsilon": mean_epsilon,
        "std_epsilon": std_epsilon,
        "fold_results": fold_results,
    }


def rmsd(predicted: list[float], experimental: list[float]) -> float:
    """Root-mean-square deviation."""
    return float(np.sqrt(np.mean((np.array(predicted) - np.array(experimental)) ** 2)))


def mean_abs_deviation(predicted: list[float], experimental: list[float]) -> float:
    """Mean absolute deviation."""
    return float(np.mean(np.abs(np.array(predicted) - np.array(experimental))))


def max_abs_deviation(predicted: list[float], experimental: list[float]) -> float:
    """Maximum absolute deviation."""
    return float(np.max(np.abs(np.array(predicted) - np.array(experimental))))


def r_squared(predicted: list[float], experimental: list[float], p: int = 0) -> float:
    """Coefficient of determination (adjusted R² if p > 0).

    When p (number of fitted parameters) is provided, returns the adjusted R²:
        R²_adj = 1 - (1 - R²) * (n - 1) / (n - p - 1)
    """
    exp = np.array(experimental)
    n = len(exp)
    ss_res = np.sum((np.array(predicted) - exp) ** 2)
    ss_tot = np.sum((exp - np.mean(exp)) ** 2)
    r2 = 1.0 - ss_res / ss_tot
    if p > 0 and n > p + 1:
        r2 = 1.0 - (1.0 - r2) * (n - 1) / (n - p - 1)
    return float(r2)
