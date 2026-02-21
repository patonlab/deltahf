"""xTB subprocess wrapper for geometry optimization and energy extraction."""

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

_ENERGY_RE = re.compile(r"total energy\s+([-\d.]+)\s+Eh", re.IGNORECASE)


@dataclass
class XtbResult:
    energy: float  # Total energy in Hartree
    optimized_xyz_path: Path | None
    converged: bool
    stdout: str
    wbo_path: Path | None = None


def find_xtb_binary() -> str:
    """Locate the xtb binary. Raises FileNotFoundError if not found."""
    path = shutil.which("xtb")
    if path is None:
        raise FileNotFoundError("xtb binary not found. Install via: conda install -c conda-forge xtb")
    return path


def find_gxtb_binary() -> str | None:
    """Locate the gxtb binary. Returns None if not found (gxtb is optional)."""
    return shutil.which("gxtb")


def build_xtb_command(
    xyz_path: str,
    opt: bool = True,
    gfn: int = 2,
    charge: int = 0,
    uhf: int = 0,
    parallel: int | None = None,
) -> list[str]:
    """Build xtb command line arguments."""
    cmd = ["xtb", str(xyz_path)]
    if opt:
        cmd.append("--opt")
    cmd.extend(["--gfn", str(gfn)])
    cmd.extend(["--chrg", str(charge)])
    if uhf > 0:
        cmd.extend(["--uhf", str(uhf)])
    if parallel is not None:
        cmd.extend(["--parallel", str(parallel)])
    return cmd


def parse_total_energy(output: str) -> float:
    """Parse the total energy (in Hartree) from xTB stdout.

    Looks for the last occurrence of 'total energy' in the SUMMARY block.
    """
    matches = _ENERGY_RE.findall(output)
    if not matches:
        excerpt = output[-500:] if len(output) > 500 else output
        raise RuntimeError(f"Could not parse xTB total energy. Output tail:\n{excerpt}")
    return float(matches[-1])


def parse_wbo_file(wbo_path: Path) -> dict[tuple[int, int], float]:
    """Parse xTB Wiberg bond order file.

    The wbo file has lines of the form: atom_i  atom_j  bond_order
    where atom indices are 1-based. Returns a dict mapping (i, j) pairs
    (converted to 0-based) to bond orders.
    """
    wbos: dict[tuple[int, int], float] = {}
    for line in wbo_path.read_text().strip().splitlines():
        parts = line.split()
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        bo = float(parts[2])
        wbos[(i, j)] = bo
        wbos[(j, i)] = bo
    return wbos


def run_xtb_optimization(
    xyz_path: Path,
    gfn: int = 2,
    charge: int = 0,
    timeout: int = 300,
    parallel: int | None = None,
) -> XtbResult:
    """Run xTB geometry optimization and return results."""
    import os
    find_xtb_binary()
    cmd = build_xtb_command(str(xyz_path), opt=True, gfn=gfn, charge=charge, parallel=parallel)

    work_dir = xyz_path.parent

    env = os.environ.copy()
    if parallel is not None:
        # Constrain both OpenMP and BLAS thread pools to avoid oversubscription
        env["OMP_NUM_THREADS"] = str(parallel)
        env["MKL_NUM_THREADS"] = str(parallel)
        env["OPENBLAS_NUM_THREADS"] = str(parallel)

    try:
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            env=env,
        )
    except subprocess.TimeoutExpired:
        return XtbResult(energy=float("nan"), optimized_xyz_path=None, converged=False, stdout="Timeout expired")

    if result.returncode != 0:
        return XtbResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout=result.stdout + result.stderr,
        )

    energy = parse_total_energy(result.stdout)
    opt_xyz = work_dir / "xtbopt.xyz"
    wbo_file = work_dir / "wbo"

    return XtbResult(
        energy=energy,
        optimized_xyz_path=opt_xyz if opt_xyz.exists() else None,
        converged=True,
        stdout=result.stdout,
        wbo_path=wbo_file if wbo_file.exists() else None,
    )


def parse_gxtb_energy_file(energy_file_path: Path) -> float:
    """Parse gxtb energy file.

    The energy file has the format:
        energy
        <value1> <energy_in_hartree>
        $end

    Returns the energy in Hartree (second column of line 1).
    """
    lines = energy_file_path.read_text().strip().splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"gxtb energy file has insufficient lines: {len(lines)}")

    data_line = lines[1].strip().split()
    if len(data_line) < 2:
        raise RuntimeError(f"gxtb energy file line 1 has insufficient columns: {len(data_line)}")

    return float(data_line[1])


def run_gxtb_single_point(
    xyz_path: Path,
    timeout: int = 300,
) -> XtbResult:
    """Run gxtb single-point energy calculation and return results.

    gxtb does not support geometry optimization (no analytical gradients).
    It reads an XYZ file and writes energy to an 'energy' file in the working directory.
    """
    gxtb_path = find_gxtb_binary()
    if gxtb_path is None:
        return XtbResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout="gxtb binary not found",
        )

    cmd = ["gxtb", "-c", xyz_path.name]
    work_dir = xyz_path.parent

    try:
        result = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return XtbResult(energy=float("nan"), optimized_xyz_path=None, converged=False, stdout="Timeout expired")

    if result.returncode != 0:
        return XtbResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout=result.stdout + result.stderr,
        )

    energy_file = work_dir / "energy"
    if not energy_file.exists():
        return XtbResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout=result.stdout + "\nError: gxtb did not produce 'energy' file",
        )

    try:
        energy = parse_gxtb_energy_file(energy_file)
    except (ValueError, RuntimeError) as e:
        return XtbResult(
            energy=float("nan"),
            optimized_xyz_path=None,
            converged=False,
            stdout=result.stdout + f"\nError parsing energy file: {e}",
        )

    return XtbResult(
        energy=energy,
        optimized_xyz_path=xyz_path,  # gxtb doesn't modify geometry
        converged=True,
        stdout=result.stdout,
    )
