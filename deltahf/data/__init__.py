"""Training data loading utilities."""

from importlib import resources

import pandas as pd


def load_training_data() -> pd.DataFrame:
    """Load the training dataset of 103 molecules with experimental ΔHf° values."""
    data_dir = resources.files("deltahf.data")
    csv_path = data_dir / "training_data.csv"
    return pd.read_csv(csv_path)
