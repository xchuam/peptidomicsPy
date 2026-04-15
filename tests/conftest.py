from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from peptidomicspy import processPeptides

matplotlib.use("Agg")


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "v1.1.1-alpha"
REFERENCE_DIR = REPO_ROOT / "tests" / "reference" / "v1.1.1-alpha"


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    return FIXTURE_DIR


@pytest.fixture(scope="session")
def reference_dir() -> Path:
    return REFERENCE_DIR


@pytest.fixture(scope="session")
def reference_manifest(reference_dir: Path) -> dict[str, object]:
    return json.loads((reference_dir / "manifest.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="session")
def process_result(fixture_dir: Path):
    return processPeptides(
        peptides_file=fixture_dir / "Yogurtexample_QR188-205.csv",
        intensity_columns_file=fixture_dir / "Intensity_columns.csv",
        protein_mapping_file=fixture_dir / "protein_mapping.csv",
    )


@pytest.fixture(scope="session")
def reference_frame(reference_dir: Path):
    def _load(relative_path: str) -> pd.DataFrame:
        return pd.read_csv(reference_dir / relative_path)

    return _load


def normalize_frame(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy(deep=True)
    for col in out.columns:
        if isinstance(out[col].dtype, pd.CategoricalDtype):
            out[col] = out[col].astype(str)
    return out


def assert_dataframe_close(actual: pd.DataFrame, reference: pd.DataFrame) -> None:
    actual_norm = normalize_frame(actual)
    reference_norm = normalize_frame(reference)
    assert list(actual_norm.columns) == list(reference_norm.columns)

    sort_cols = list(reference_norm.columns)
    actual_norm = actual_norm.sort_values(by=sort_cols, kind="mergesort").reset_index(drop=True)
    reference_norm = reference_norm.sort_values(by=sort_cols, kind="mergesort").reset_index(drop=True)

    numeric_cols = [
        col
        for col in reference_norm.columns
        if pd.api.types.is_numeric_dtype(reference_norm[col]) or pd.api.types.is_numeric_dtype(actual_norm[col])
    ]
    other_cols = [col for col in reference_norm.columns if col not in numeric_cols]

    if other_cols:
        assert_frame_equal(
            actual_norm[other_cols],
            reference_norm[other_cols],
            check_dtype=False,
            check_like=False,
        )
    for col in numeric_cols:
        left = pd.to_numeric(actual_norm[col], errors="coerce").to_numpy()
        right = pd.to_numeric(reference_norm[col], errors="coerce").to_numpy()
        np.testing.assert_allclose(left, right, rtol=1e-9, atol=1e-9, equal_nan=True)


def get_plot_label(plot, key: str) -> str | None:
    labels = getattr(plot, "labels", None)
    if labels is None:
        return None
    for candidate in (key, "colour" if key == "color" else key):
        if hasattr(labels, candidate):
            value = getattr(labels, candidate)
            if value is not None:
                return value
        if hasattr(labels, "get"):
            value = labels.get(candidate)
            if value is not None:
                return value
    return None
