from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def read_input_table(data: str | Path | pd.DataFrame) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data.copy(deep=True)
    path = Path(data)
    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {path}")
    return pd.read_csv(
        path,
        sep=None,
        engine="python",
        encoding="utf-8-sig",
    )


def normalize_column_names(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy(deep=True)
    out.columns = [str(col).replace(" ", ".") for col in out.columns]
    return out


def to_ordered_category(frame: pd.DataFrame, column: str, levels: Iterable[object]) -> None:
    base_levels = list(levels)
    raw_values = pd.Series(frame[column].astype("object"), index=frame.index)
    observed = [value for value in raw_values.dropna().unique().tolist() if value not in base_levels]
    frame[column] = pd.Categorical(raw_values, categories=[*base_levels, *observed], ordered=True)


def scientific_10_label(value: float | int | None) -> str:
    if value is None:
        return ""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return ""
    if np.isnan(numeric):
        return ""
    if numeric > 10:
        formatted = f"{numeric:.2e}"
        mantissa, exponent = formatted.split("e")
        mantissa = mantissa.rstrip("0").rstrip(".")
        return rf"${mantissa} \times 10^{{{int(exponent)}}}$"
    if numeric.is_integer():
        return str(int(numeric))
    return f"{numeric:g}"


def scientific_10_labels(values: Iterable[float]) -> list[str]:
    return [scientific_10_label(value) for value in values]


def split_facet_expression(expr: str | None) -> list[str]:
    if not expr or expr == ".":
        return []
    parts = [piece.strip() for piece in expr.split("+")]
    return [piece for piece in parts if piece]


def interaction_label(frame: pd.DataFrame, columns: list[str], sep: str) -> pd.Series:
    if not columns:
        return pd.Series(["Sample"] * len(frame), index=frame.index, dtype="object")
    values = frame[columns].astype(str)
    return values.agg(sep.join, axis=1)


def bh_adjust(pvalues: Iterable[float]) -> np.ndarray:
    arr = np.asarray(list(pvalues), dtype=float)
    out = np.full(arr.shape, np.nan)
    valid = np.isfinite(arr)
    if not valid.any():
        return out
    vals = arr[valid]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    cumulative = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        candidate = ranked[idx] * n / rank
        cumulative = min(cumulative, candidate)
        adjusted[idx] = cumulative
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    out[valid] = restored
    return out


def adjust_pvalues(pvalues: Iterable[float], method: str = "BH") -> np.ndarray:
    method_norm = method.lower()
    aliases = {
        "fdr": "bh",
    }
    method_norm = aliases.get(method_norm, method_norm)
    arr = np.asarray(list(pvalues), dtype=float)
    out = np.full(arr.shape, np.nan)
    valid = np.isfinite(arr)
    if not valid.any():
        return out
    vals = arr[valid]
    n = len(vals)

    if method_norm in {"bh"}:
        out[valid] = bh_adjust(vals)
        return out

    if method_norm == "bonferroni":
        out[valid] = np.clip(vals * n, 0.0, 1.0)
        return out

    if method_norm == "none":
        out[valid] = vals
        return out

    if method_norm == "holm":
        order = np.argsort(vals)
        ranked = vals[order]
        adjusted = np.empty(n, dtype=float)
        running = 0.0
        for idx, value in enumerate(ranked):
            candidate = (n - idx) * value
            running = max(running, candidate)
            adjusted[idx] = running
        adjusted = np.clip(adjusted, 0.0, 1.0)
        restored = np.empty(n, dtype=float)
        restored[order] = adjusted
        out[valid] = restored
        return out

    if method_norm == "by":
        harmonic = float(np.sum(1.0 / np.arange(1, n + 1)))
        out[valid] = np.clip(bh_adjust(vals) * harmonic, 0.0, 1.0)
        return out

    raise ValueError("Unsupported adjust method. Supported methods: BH, fdr, bonferroni, holm, BY, none.")
