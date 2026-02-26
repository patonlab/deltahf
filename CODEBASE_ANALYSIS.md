# Codebase Analysis: Readability, Testing, and CI

Analysis of the `deltahf` codebase (v0.1.0) — a Python CLI tool for predicting
gas-phase heats of formation using atom equivalent energies and semi-empirical
quantum chemistry.

**Codebase:** ~2,300 lines of source across 10 modules; ~1,700 lines of tests
across 13 test files.

---

## 1. Readability

### Strengths

The codebase is well-written overall:

- **Type hints**: Comprehensive modern Python 3.10+ annotations throughout,
  including union types (`str | None`), generics (`list[str]`), and
  `TYPE_CHECKING` for forward references (`pipeline.py:34`).
- **Docstrings**: Google-style docstrings on all public functions with `Args`,
  `Returns`, and `Raises` sections. Algorithm choices are explained inline
  (e.g., `smiles.py:97-105` on Kekulization rationale).
- **Naming**: Clear, consistent `snake_case` for functions and variables,
  `UPPER_CASE` for constants, descriptive names throughout. No single-letter
  variables in business logic.
- **Module separation**: Each module has a single responsibility — `smiles.py`
  for atom classification, `conformers.py` for 3D generation, `xtb.py` for
  subprocess execution, `pipeline.py` for orchestration.
- **Import organization**: Standard library, third-party, and local imports
  are cleanly separated with blank lines.

### Suggested Improvements

#### 1.1 Extract magic numbers to named constants

Several numeric thresholds are used inline across multiple functions. Centralizing
them in `constants.py` would improve discoverability and prevent silent
inconsistencies if they need to change.

| Value | Location(s) | Suggested constant |
|-------|-------------|-------------------|
| `1.25` | `smiles.py:134`, `smiles.py:440`, `smiles.py:446` | `MULTIPLY_BONDED_THRESHOLD` |
| `0.125` | `conformers.py:50` | `RMSD_PRUNE_THRESHOLD` |
| `50` | `__main__.py:115` | `TABLE_WIDTH` |
| `82` | `__main__.py:142` | `OUTLIER_TABLE_WIDTH` |

#### 1.2 Reduce repetition in atom classification functions

The eight `classify_atoms_*` functions in `smiles.py` share a common structure:
parse SMILES → initialize zeroed counts → iterate atoms → classify by
element/H/F/Cl guards → model-specific logic → return counts. The per-function
boilerplate (element guard, H short-circuit, F/Cl handling) is duplicated ~8
times.

A possible refactor: extract the common iteration pattern into a helper, letting
each model supply only its classification function:

```python
def _classify_atoms(smiles: str, param_names: list[str], classify_fn) -> dict[str, int]:
    """Common atom classification loop."""
    mol = smiles_to_mol(smiles)
    counts = _zero_counts(param_names)
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in SUPPORTED_ELEMENTS:
            continue
        classify_fn(atom, symbol, counts)
    return counts
```

This is low priority — the current code is readable as-is and the duplication
doesn't cause bugs. But it would reduce ~150 lines and make new model additions
less error-prone.

#### 1.3 Reduce duplication in `uma.py` error handling

Four near-identical `except Exception` blocks in `run_mlip_optimization()`
(`uma.py:144-150`, `162-168`, `177-183`, `189-195`) each construct the same
`UmaResult(energy=float("nan"), ...)`. A small helper would reduce this:

```python
def _fail(msg: str) -> UmaResult:
    return UmaResult(energy=float("nan"), optimized_xyz_path=None, converged=False, stdout=msg)
```

#### 1.4 Consider splitting `__main__.py`

At 641 lines, `__main__.py` bundles argument parsing, model result printing,
and the full `cmd_fit()` / `cmd_predict()` orchestration. Splitting into
`cli.py` (parser definition) and keeping `__main__.py` as a thin entry point
would improve navigability. Low priority since the file is well-sectioned.

#### 1.5 Unused variable in `uma.py`

`log_lines` is initialized at `uma.py:171` but never appended to — the
`UmaResult.stdout` on the success path (`uma.py:201`) is always empty.
Either populate it with optimizer diagnostics or remove it.

---

## 2. Testing

### Strengths

- **13 test modules** covering all core modules with clear class-based
  organization (`TestSmilesToMol`, `TestCountAtoms`, etc.).
- **`pytest` with markers**: Integration tests requiring `xtb` CLI are properly
  marked (`@pytest.mark.integration`) and skipped when unavailable.
- **Synthetic data fixtures**: `test_atom_equivalents.py:96-116` uses
  deterministic synthetic data to verify fitting correctness.
- **Mock-based unit tests**: `test_pipeline.py` patches `run_xtb_optimization`
  to test pipeline logic without requiring xTB.
- **Cache round-trip tests**: Persistence, backward compatibility with old
  cache formats, and SMILES canonicalization are all tested.
- **Parametric molecule tests**: Multiple molecules tested per classification
  function (methane, benzene, TNT, HCN, etc.) covering sp/sp2/sp3, aromatic,
  and charge-separated cases.

### Suggested Improvements

#### 2.1 Add test coverage reporting

There is no coverage measurement configured. Adding `pytest-cov` reveals which
code paths are untested:

```toml
# pyproject.toml
[tool.pytest.ini_options]
addopts = "--cov=deltahf --cov-report=term-missing"
```

```yaml
# environment.yml — add to dependencies
- pytest-cov
```

#### 2.2 Missing test coverage for key modules

| Module | Current coverage | Gap |
|--------|-----------------|-----|
| `uma.py` | **No unit tests** | `load_mlip_model()`, `run_mlip_optimization()`, `_cuda_available()` are untested. Mock ASE/FAIRChem imports to test error paths and success flow. |
| `pipeline.py:process_csv()` | **No direct tests** | Only `process_molecule()` is tested. `process_csv()` handles tqdm iteration and DataFrame construction. |
| `__main__.py:cmd_fit()` | **No tests beyond parser** | `test_cli.py` only tests argument parsing, not actual fit/predict execution. A smoke test with a tiny 5-molecule CSV would catch regressions. |
| `conformers.py:check_connectivity()` | **Tested in `test_isomerization.py`** but only via known SMILES. Edge cases (unreadable XYZ, atom count mismatch) are not tested directly. |
| `atom_equivalents.py:r_squared()` | **No tests** for the adjusted R² branch (`p > 0`). |
| `atom_equivalents.py:mean_abs_deviation()` | **No tests**. |

#### 2.3 Add `pytest.mark.parametrize` for classification functions

The `TestClassifyAtoms*` classes have many test methods that share the same
pattern (call function, assert counts). Using `@pytest.mark.parametrize` would
be more concise and make it easier to add new test molecules:

```python
@pytest.mark.parametrize("smiles, expected_C, expected_C_prime", [
    ("CC", 2, 0),
    ("C=C", 0, 2),
    ("C#C", 0, 2),
    ("c1ccccc1", 0, 6),
])
def test_classify_7param_carbon(smiles, expected_C, expected_C_prime):
    counts = classify_atoms_7param(smiles)
    assert counts["C"] == expected_C
    assert counts["C_prime"] == expected_C_prime
```

#### 2.4 Missing tests for `bondorder`, `bondorder_ext`, `bondorder_ar`, `neighbour` classifiers

`test_smiles.py` imports and tests `count_atoms`, `classify_atoms_7param`,
`classify_atoms_hybrid`, and `classify_atoms_extended`, but does **not** import
or test:
- `classify_atoms_bondorder`
- `classify_atoms_bondorder_ext`
- `classify_atoms_bondorder_ar`
- `classify_atoms_neighbour`

These four classifiers (`smiles.py:253-427`) have zero direct unit test
coverage. `test_atom_counts.py` partially covers them but only for simple
molecules.

#### 2.5 Add edge case tests

- **Empty/whitespace SMILES** for all classification functions (not just
  `smiles_to_mol`)
- **Charged molecules** (e.g., `[NH4+]`, `[O-]c1ccccc1`) for charge-handling
  paths
- **Molecules with all supported halogens** (F, S, Cl) through all 8
  classifiers to verify halogen handling consistency
- **Single-atom molecules** (`[H][H]`, hydrogen gas) for conformer generation
  edge cases

#### 2.6 Integration tests only run in CircleCI

Integration tests (`@pytest.mark.integration`) require the `xtb` binary, which
is only available in the CircleCI environment. GitHub Actions — the primary
CI for PRs — skips them entirely. This means integration regressions can merge
into `main` undetected (see CI section below).

---

## 3. CI/CD

### Current Setup

| System | Trigger | Steps |
|--------|---------|-------|
| **GitHub Actions** (`.github/workflows/ci.yml`) | Push/PR to `main` | Lint (ruff) + unit tests only |
| **CircleCI** (`.circleci/config.yml`) | Push to `main` only | Lint (ruff) + unit tests + integration tests |

### Suggested Improvements

#### 3.1 Unify CI: run integration tests in GitHub Actions

The most critical gap: **integration tests only run _after_ merging to `main`**
via CircleCI. PRs are only validated by unit tests in GitHub Actions. This means
a broken xTB integration can merge undetected.

Fix: install `xtb` in the GitHub Actions workflow (it's already in
`environment.yml` via conda-forge) and run the full suite:

```yaml
      - name: Unit tests
        run: pytest -m "not integration" -v --tb=short

      - name: Integration tests
        run: pytest -m integration -v --tb=short
```

Consider making CircleCI redundant or using it for additional checks
(nightly, performance benchmarks).

#### 3.2 Add Python version matrix testing

The project requires `>=3.10` but only tests on one Python version (whatever
`ubuntu-latest` ships). A matrix ensures compatibility:

```yaml
    strategy:
      matrix:
        python-version: ["3.10", "3.11", "3.12"]
    steps:
      - uses: conda-incubator/setup-miniconda@v3
        with:
          python-version: ${{ matrix.python-version }}
```

#### 3.3 Add test coverage reporting

Integrate `pytest-cov` with a coverage threshold to prevent regression:

```yaml
      - name: Tests with coverage
        run: pytest --cov=deltahf --cov-report=xml --cov-fail-under=70
```

Optionally upload to Codecov or Coveralls for PR-level coverage diffs.

#### 3.4 Add pre-commit hooks

There is no `.pre-commit-config.yaml`, so linting violations can be committed
locally and only caught in CI. Adding pre-commit ensures fast local feedback:

```yaml
# .pre-commit-config.yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.8.0
    hooks:
      - id: ruff
        args: [--fix]
      - id: ruff-format
```

#### 3.5 Add `ruff format` for consistent code formatting

Currently only `ruff check` (linting) is configured. There is no formatter
enforced, so code style can drift between contributors. Adding `ruff format`
to CI (and pre-commit) ensures consistent formatting:

```yaml
      - name: Format check
        run: ruff format --check deltahf/ tests/
```

#### 3.6 Add type checking with mypy

The codebase already has comprehensive type annotations, but they are never
validated. Adding `mypy` to CI would catch type errors statically:

```toml
# pyproject.toml
[tool.mypy]
python_version = "3.10"
warn_return_any = true
warn_unused_configs = true
ignore_missing_imports = true  # for rdkit, ase, fairchem
```

```yaml
      - name: Type check
        run: mypy deltahf/
```

#### 3.7 Pin dependency versions for reproducible CI

Dependencies in `environment.yml` and `pyproject.toml` are unpinned (e.g.,
`rdkit`, `numpy`, `pandas`). A surprise upstream release could break CI
silently. Options:

- Add a `conda-lock.yml` or `requirements-lock.txt` for pinned versions
- Use Dependabot/Renovate for automated dependency update PRs
- At minimum, pin major versions in `pyproject.toml`:
  `"numpy>=1.24,<3"`, `"pandas>=2.0,<3"`

#### 3.8 Add a Makefile for common development tasks

There is no automation for common developer commands. A `Makefile` reduces
friction:

```makefile
.PHONY: lint test test-all format typecheck

lint:
	ruff check deltahf/ tests/

format:
	ruff format deltahf/ tests/

test:
	pytest -m "not integration" -v

test-all:
	pytest -v

typecheck:
	mypy deltahf/
```

#### 3.9 No automated release process

The project is at `v0.1.0` with no mechanism to bump versions, tag releases,
or publish to PyPI. If distribution is a goal, consider:

- `python-semantic-release` or manual `bump2version` for version management
- A GitHub Actions workflow triggered on tags for PyPI publishing
- GitHub Releases for changelog generation

---

## Summary: Priority Ranking

| # | Category | Improvement | Impact | Effort |
|---|----------|------------|--------|--------|
| 1 | CI | Run integration tests in GitHub Actions PRs | **High** — prevents broken merges | Low |
| 2 | Testing | Add unit tests for `uma.py` | **High** — entirely untested module | Medium |
| 3 | CI | Add `pytest-cov` with coverage threshold | **High** — visibility into gaps | Low |
| 4 | Testing | Add tests for bondorder/neighbour classifiers | **Medium** — 4 untested classifiers | Medium |
| 5 | CI | Add pre-commit hooks (ruff lint + format) | **Medium** — faster feedback loop | Low |
| 6 | CI | Python version matrix (3.10, 3.11, 3.12) | **Medium** — compatibility assurance | Low |
| 7 | CI | Add `ruff format` check | **Medium** — consistent formatting | Low |
| 8 | Readability | Extract magic numbers to constants | **Low** — minor clarity gain | Low |
| 9 | CI | Add mypy type checking | **Low** — annotations exist, just needs validation | Low |
| 10 | Readability | Reduce classify_atoms_* duplication | **Low** — code is already readable | Medium |
| 11 | Testing | CLI smoke tests (`cmd_fit`, `cmd_predict`) | **Medium** — catches CLI regressions | Medium |
| 12 | CI | Pin dependency versions | **Low** — prevents rare breakage | Low |
