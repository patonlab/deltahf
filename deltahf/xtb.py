"""xTB subprocess wrapper for geometry optimization and energy extraction."""

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

HARTREE_TO_KCAL = 627.5094740631


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


def build_xtb_command(
    xyz_path: str,
    opt: bool = True,
    gfn: int = 2,
    charge: int = 0,
    uhf: int = 0,
) -> list[str]:
    """Build xtb command line arguments."""
    cmd = ["xtb", str(xyz_path)]
    if opt:
        cmd.append("--opt")
    cmd.extend(["--gfn", str(gfn)])
    cmd.extend(["--chrg", str(charge)])
    if uhf > 0:
        cmd.extend(["--uhf", str(uhf)])
    return cmd


def parse_total_energy(output: str) -> float:
    """Parse the total energy (in Hartree) from xTB stdout.

    Looks for the last occurrence of 'total energy' in the SUMMARY block.
    """
    pattern = r"total energy\s+([-\d.]+)\s+Eh"
    matches = re.findall(pattern, output, re.IGNORECASE)
    if not matches:
        raise RuntimeError("Could not parse total energy from xTB output")
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
) -> XtbResult:
    """Run xTB geometry optimization and return results."""
    find_xtb_binary()
    cmd = build_xtb_command(str(xyz_path), opt=True, gfn=gfn, charge=charge)

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
