"""Shared test fixtures and markers."""

import shutil

import pytest


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line("markers", "integration: marks tests that require xtb CLI")


@pytest.fixture
def xtb_available():
    """Skip test if xtb CLI is not installed."""
    if shutil.which("xtb") is None:
        pytest.skip("xtb CLI not available")
