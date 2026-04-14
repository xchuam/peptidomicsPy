from __future__ import annotations

from importlib.resources import files

import pandas as pd


def load_example_protein_mapping() -> pd.DataFrame:
    data_path = files("peptidomicspy.data").joinpath("peptidomics_protein_mapping_example.csv")
    return pd.read_csv(data_path)

