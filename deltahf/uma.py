"""MLIP geometry optimization via ASE (UMA, eSEN, AIMNet2)."""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from pathlib import Path

# Suppress noisy warnings from MLIP dependencies
warnings.filterwarnings("ignore", category=UserWarning, module="torchtnt")
warnings.filterwarnings("ignore", category=DeprecationWarning, message=".*warp\\.vec.*")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="warp")

from deltahf.constants import HARTREE_TO_EV  # noqa: E402

DEFAULT_MAX_STEPS = 250


@dataclass
class UmaResult:
    energy: float  # Total energy in Hartree
    optimized_xyz_path: Path | None
    converged: bool
    stdout: str  # Log/diagnostic string


def _load_dotenv() -> None:
    """Load .env from the current working directory if python-dotenv is available."""
    try:
        from dotenv import load_dotenv
        load_dotenv()
    except ImportError:
        pass


def _get_model_dir() -> str | None:
    """Return MODEL_DIR from environment (loading .env first), expanding ~ if needed."""
    _load_dotenv()
    model_dir = os.getenv("MODEL_DIR")
    if model_dir is not None:
        model_dir = os.path.expanduser(model_dir)
    return model_dir


def load_mlip_model(model: str):
    """Load and return the MLIP predictor object.

    For FAIRChem models (uma, esen): loads a predictor unit from MODEL_DIR.
    For AIMNet2: returns None (calculator is instantiated per-molecule).

    Raises ImportError or FileNotFoundError with clear install instructions if
    required packages or model files are not available.
    """
    if model == "aimnet2":
        try:
            from aimnet.calculators import AIMNet2ASE  # noqa: F401
        except ImportError:
            raise ImportError(
                "AIMNet2 not available. Install with: pip install aimnet2calc"
            )
        return None  # Calculator is created fresh per molecule

    # FAIRChem models (uma, esen)
    try:
        from fairchem.core import FAIRChemCalculator  # noqa: F401
        from fairchem.core.units.mlip_unit import load_predict_unit
    except ImportError:
        raise ImportError(
            f"FAIRChem not available (required for --optimizer {model}). "
            "Install with: pip install fairchem-core"
        )

    model_dir = _get_model_dir()
    if model_dir is None:
        raise FileNotFoundError(
            f"MODEL_DIR environment variable not set (required for --optimizer {model}). "
            "Create a .env file or set MODEL_DIR=/path/to/fairchem/models"
        )

    model_filenames = {
        "uma": "uma-s-1p1.pt",
        "esen": "esen_sm_conserving_all.pt",
    }
    filename = model_filenames[model]
    model_path = os.path.join(model_dir, filename)

    if not os.path.exists(model_path):
        raise FileNotFoundError(
            f"Model file not found: {model_path}\n"
            f"Download the {model.upper()} model and place it in MODEL_DIR={model_dir}"
        )

    device = "cuda" if _cuda_available() else "cpu"
    print(f"   Device: {device}")
    predictor = load_predict_unit(model_path, device=device)
    return predictor


def _cuda_available() -> bool:
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False


def run_mlip_optimization(
    xyz_path: Path,
    model: str,
    predictor,
    fmax: float = 0.05,
    timeout: int = 300,
    charge: int = 0,
    spin: int = 1,
) -> UmaResult:
    """Optimize geometry using an MLIP via ASE's L-BFGS optimizer.

    Reads xyz_path, attaches the appropriate ASE calculator, runs L-BFGS
    to convergence, writes the optimized geometry to uma_opt.xyz in the
    same directory, and returns energy (Hartree) and path.

    Args:
        xyz_path: Path to input XYZ file.
        model: One of "uma", "esen", "aimnet2".
        predictor: Pre-loaded FAIRChem predictor (from load_mlip_model), or
                   None for AIMNet2.
        fmax: Force convergence criterion (eV/Å).
        timeout: Not enforced by ASE; reserved for future use.
        charge: Total formal charge on the system (default 0).
        spin: Spin multiplicity (2S+1); 1 = singlet (default), 2 = doublet, etc.

    Returns:
        UmaResult with energy in Hartree and path to optimized XYZ.
    """
    try:
        from ase.io import read, write
        from ase.optimize import LBFGS
    except ImportError:
        raise ImportError("ASE not available. Install with: pip install ase")

    try:
        atoms = read(str(xyz_path))
    except Exception as e:
        return UmaResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout=f"Failed to read XYZ: {e}",
        )

    # Attach calculator
    try:
        if model == "aimnet2":
            from aimnet.calculators import AIMNet2ASE
            atoms.calc = AIMNet2ASE("aimnet2")
        else:
            from fairchem.core import FAIRChemCalculator
            atoms.info["charge"] = charge
            atoms.info["spin"] = spin
            atoms.calc = FAIRChemCalculator(predictor, task_name="omol")
    except Exception as e:
        return UmaResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout=f"Failed to set up {model} calculator: {e}",
        )

    # Run L-BFGS optimization
    log_lines = []
    try:
        opt = LBFGS(atoms, logfile=None)
        opt.run(fmax=fmax, steps=DEFAULT_MAX_STEPS)
        energy_ev = atoms.get_potential_energy()
        energy_hartree = energy_ev / HARTREE_TO_EV
    except Exception as e:
        return UmaResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout=f"{model} optimization failed: {e}",
        )

    # Write optimized geometry
    opt_xyz_path = xyz_path.parent / "uma_opt.xyz"
    try:
        write(str(opt_xyz_path), atoms, format="xyz")
    except Exception as e:
        return UmaResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout=f"Failed to write optimized XYZ: {e}",
        )

    return UmaResult(
        energy=energy_hartree,
        optimized_xyz_path=opt_xyz_path,
        converged=True,
        stdout="\n".join(log_lines),
    )
