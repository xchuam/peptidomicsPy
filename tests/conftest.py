from __future__ import annotations

import json
import os
import subprocess
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
REFERENCE_DIR = REPO_ROOT / ".migration-artifacts" / "v1.1.1-alpha" / "reference_outputs"
DEFAULT_UPSTREAM_R_DIR = Path("/tmp/peptidomicsR-v1.1.1-alpha")


@pytest.fixture(scope="session")
def fixture_dir() -> Path:
    return FIXTURE_DIR


@pytest.fixture(scope="session")
def reference_dir() -> Path:
    return REFERENCE_DIR


@pytest.fixture(scope="session")
def upstream_r_dir() -> Path:
    candidate = Path(os.environ.get("PEPTIDOMICSR_UPSTREAM_DIR", DEFAULT_UPSTREAM_R_DIR))
    if not candidate.exists():
        pytest.skip(f"Upstream R package checkout not found at {candidate}")
    return candidate


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
def r_parity_dir(tmp_path_factory, upstream_r_dir: Path, fixture_dir: Path) -> Path:
    out_dir = tmp_path_factory.mktemp("r-parity")
    env = os.environ.copy()
    env.setdefault("R_LIBS_USER", str(Path.home() / "R" / "x86_64-pc-linux-gnu-library" / "4.3"))
    subprocess.run(
        [
            "Rscript",
            "tools/generate_r_parity_suite.R",
            str(upstream_r_dir),
            str(fixture_dir),
            str(out_dir),
        ],
        cwd=REPO_ROOT,
        env=env,
        check=True,
    )
    return out_dir


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


def load_rgb_image(path: Path) -> np.ndarray:
    arr = matplotlib.image.imread(path).astype(float)
    if arr.max() > 1:
        arr = arr / 255.0
    if arr.shape[-1] == 4:
        alpha = arr[..., 3:4]
        arr = arr[..., :3] * alpha + (1 - alpha)
    return arr[..., :3]


def resize_nn(image: np.ndarray, height: int, width: int) -> np.ndarray:
    ys = np.linspace(0, image.shape[0] - 1, height).round().astype(int)
    xs = np.linspace(0, image.shape[1] - 1, width).round().astype(int)
    return image[np.ix_(ys, xs)]


def image_mae(path_a: Path, path_b: Path) -> float:
    image_a = load_rgb_image(path_a)
    image_b = load_rgb_image(path_b)
    height = max(image_a.shape[0], image_b.shape[0])
    width = max(image_a.shape[1], image_b.shape[1])
    image_a = resize_nn(image_a, height, width)
    image_b = resize_nn(image_b, height, width)
    return float(np.abs(image_a - image_b).mean())
