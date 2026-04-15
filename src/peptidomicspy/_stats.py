from __future__ import annotations

import warnings
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats

from ._result import PeptidomicsResult
from ._utils import adjust_pvalues


def _normalize_comparisons(
    comparisons: list[tuple[str, str]] | tuple[str, str] | list[list[str]] | list[str],
) -> list[tuple[str, str]]:
    if isinstance(comparisons, tuple) and len(comparisons) == 2:
        return [(str(comparisons[0]), str(comparisons[1]))]
    if isinstance(comparisons, list) and len(comparisons) == 2 and all(
        isinstance(item, str) for item in comparisons
    ):
        return [(comparisons[0], comparisons[1])]
    normalized: list[tuple[str, str]] = []
    if isinstance(comparisons, list):
        for item in comparisons:
            if not isinstance(item, (list, tuple)) or len(item) != 2:
                raise ValueError(
                    "`comparisons` must be a list of character pairs, or a character vector of length 2."
                )
            normalized.append((str(item[0]), str(item[1])))
        return normalized
    raise ValueError(
        "`comparisons` must be a list of character pairs, or a character vector of length 2."
    )


def _parse_selector(selector: str, grp_cols: list[str], levels_by_col: dict[str, set[str]]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    used: list[str] = []
    for token in selector.split("_"):
        if ":" in token or "=" in token:
            col, value = re_split_assignment(token)
            if col not in grp_cols:
                raise ValueError(
                    f"Unknown grouping column '{col}' in selector '{selector}'. Check result.grp_cols."
                )
            mapping[col] = value
            used.append(col)
            continue
        hits = [col for col in grp_cols if token in levels_by_col[col]]
        if len(hits) == 1:
            mapping[hits[0]] = token
            used.append(hits[0])
        elif len(hits) == 0:
            raise ValueError(f"'{token}' not found in any grouping column as indicated by result.grp_cols.")
        else:
            joined = ", ".join(hits)
            raise ValueError(
                f"Token '{token}' is ambiguous (matches columns: {joined}). "
                f"Disambiguate using 'col=val', e.g. 'time={token}'."
            )
    if len(used) != len(set(used)):
        raise ValueError(f"Selector '{selector}' assigns the same column more than once.")
    return mapping


def re_split_assignment(token: str) -> tuple[str, str]:
    for sep in (":", "="):
        if sep in token:
            parts = token.split(sep, 1)
            if len(parts) != 2 or not parts[0]:
                raise ValueError(f"Malformed token: '{token}'. Use 'col=val'.")
            return parts[0], parts[1]
    raise ValueError(f"Malformed token: '{token}'. Use 'col=val'.")


def _subset_by_selector(frame: pd.DataFrame, selector_map: dict[str, str]) -> pd.DataFrame:
    out = frame
    for col, value in selector_map.items():
        out = out.loc[out[col].astype(str) == str(value)]
    return out


def _safe_values(value: object) -> np.ndarray:
    if isinstance(value, list):
        return np.asarray(value, dtype=float)
    if isinstance(value, np.ndarray):
        return value.astype(float)
    return np.asarray([], dtype=float)


def _mean_and_n(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.nan, np.nan
    return float(finite.mean()), float(finite.size)


def _plain_ttest(a: np.ndarray, b: np.ndarray, equal_var: bool) -> tuple[float, float, float]:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return np.nan, np.nan, np.nan
    mean_a = float(a.mean())
    mean_b = float(b.mean())
    var_a = float(np.var(a, ddof=1))
    var_b = float(np.var(b, ddof=1))
    delta = mean_a - mean_b

    if equal_var:
        df = a.size + b.size - 2
        pooled = (((a.size - 1) * var_a) + ((b.size - 1) * var_b)) / df
        se = math_sqrt_safe(pooled * ((1 / a.size) + (1 / b.size)))
    else:
        se_sq = (var_a / a.size) + (var_b / b.size)
        se = math_sqrt_safe(se_sq)
        numerator = se_sq**2
        denominator = 0.0
        if a.size > 1 and var_a > 0:
            denominator += ((var_a / a.size) ** 2) / (a.size - 1)
        if b.size > 1 and var_b > 0:
            denominator += ((var_b / b.size) ** 2) / (b.size - 1)
        df = numerator / denominator if denominator > 0 else np.nan

    if se == 0:
        if delta == 0:
            return 0.0, df, 1.0
        return np.sign(delta) * np.inf, df, 0.0
    statistic = delta / se
    pvalue = 2 * stats.t.sf(abs(statistic), df)
    return float(statistic), float(df), float(pvalue)


def _treat_approximation(a: np.ndarray, b: np.ndarray, lfc_thresh: float) -> tuple[float, float, float]:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return np.nan, np.nan, np.nan

    mean_a = float(a.mean())
    mean_b = float(b.mean())
    var_a = float(np.var(a, ddof=1))
    var_b = float(np.var(b, ddof=1))
    delta = mean_a - mean_b
    se_sq = (var_a / a.size) + (var_b / b.size)
    se = math_sqrt_safe(se_sq)

    numerator = se_sq**2
    denominator = 0.0
    if a.size > 1 and var_a > 0:
        denominator += ((var_a / a.size) ** 2) / (a.size - 1)
    if b.size > 1 and var_b > 0:
        denominator += ((var_b / b.size) ** 2) / (b.size - 1)
    df = numerator / denominator if denominator > 0 else np.nan

    if se == 0:
        if abs(delta) <= lfc_thresh:
            return 0.0, df, 1.0
        return np.sign(delta) * np.inf, df, 0.0

    statistic = (abs(delta) - lfc_thresh) / se
    statistic = max(float(statistic), 0.0)
    pvalue = float(stats.t.sf(statistic, df))
    return float(np.sign(delta) * statistic), float(df), pvalue


def math_sqrt_safe(value: float) -> float:
    return float(np.sqrt(value)) if value > 0 else 0.0


def ttestPeptides(
    result: PeptidomicsResult,
    comparisons: list[tuple[str, str]] | tuple[str, str] | list[list[str]] | list[str],
    pseudocount: float = 1,
    test_method: str = "treat",
    lfc_thresh: float = 1,
    alpha: float = 0.05,
    equal_var: bool = False,
    adjust: str = "BH",
    min_reps_per_side: int = 2,
) -> dict[str, pd.DataFrame]:
    if test_method not in {"treat", "plain"}:
        raise ValueError("test_method must be either 'treat' or 'plain'.")

    dt = result.dt_peptides_int_ttest.copy(deep=True)
    grp_cols = list(result.grp_cols)
    id_cols = [*result.peptides_select_col_basic, "Protein.name", "Protein.group"]
    rep_cols = [col for col in dt.columns if col.startswith("R") and col[1:].isdigit()]
    if not rep_cols:
        raise ValueError("No replicate columns 'R1..Rn' found in result.dt_peptides_int_ttest")

    long = dt.melt(
        id_vars=[*id_cols, *grp_cols],
        value_vars=rep_cols,
        var_name="Replicate",
        value_name="Intensity",
    )
    long["Intensity_adj"] = np.where(long["Intensity"] == 0, pseudocount, long["Intensity"])
    long["value"] = np.log2(long["Intensity_adj"])

    levels_by_col = {
        col: set(long[col].dropna().astype(str).unique().tolist())
        for col in grp_cols
    }

    comparison_pairs = _normalize_comparisons(comparisons)
    results: dict[str, pd.DataFrame] = {}
    for sel_a_raw, sel_b_raw in comparison_pairs:
        sel_a = _parse_selector(sel_a_raw, grp_cols, levels_by_col)
        sel_b = _parse_selector(sel_b_raw, grp_cols, levels_by_col)
        if len(sel_a) < len(grp_cols) or len(sel_b) < len(grp_cols):
            warnings.warn(
                "One or both selectors omit grouping columns; replicates will be pooled across omitted dimensions.\n"
                f"  A: {', '.join(f'{k}={v}' for k, v in sel_a.items())}\n"
                f"  B: {', '.join(f'{k}={v}' for k, v in sel_b.items())}",
                stacklevel=2,
            )

        frame_a = _subset_by_selector(long, sel_a)
        frame_b = _subset_by_selector(long, sel_b)

        sum_a = (
            frame_a.groupby(id_cols, sort=False, observed=True)["value"]
            .apply(lambda s: s[np.isfinite(s)].tolist())
            .reset_index(name="A_vals")
        )
        sum_b = (
            frame_b.groupby(id_cols, sort=False, observed=True)["value"]
            .apply(lambda s: s[np.isfinite(s)].tolist())
            .reset_index(name="B_vals")
        )

        merged = sum_a.merge(sum_b, on=id_cols, how="outer")
        merged["A_vals"] = merged["A_vals"].apply(_safe_values)
        merged["B_vals"] = merged["B_vals"].apply(_safe_values)
        merged[["mean_A", "n_A"]] = merged["A_vals"].apply(lambda x: pd.Series(_mean_and_n(x)))
        merged[["mean_B", "n_B"]] = merged["B_vals"].apply(lambda x: pd.Series(_mean_and_n(x)))

        pseudocount_log2 = float(np.log2(pseudocount))
        merged = merged.loc[
            ~(
                (merged["mean_A"] == pseudocount_log2)
                & (merged["mean_B"] == pseudocount_log2)
            )
        ].copy()

        statistics: list[tuple[float, float, float]] = []
        for _, row in merged.iterrows():
            a_vals = row["A_vals"]
            b_vals = row["B_vals"]
            if len(a_vals) < min_reps_per_side or len(b_vals) < min_reps_per_side:
                statistics.append((np.nan, np.nan, np.nan))
                continue
            if test_method == "plain":
                statistics.append(_plain_ttest(a_vals, b_vals, equal_var=equal_var))
            else:
                statistics.append(_treat_approximation(a_vals, b_vals, lfc_thresh=lfc_thresh))

        stats_frame = pd.DataFrame(statistics, columns=["t", "df", "p.value"], index=merged.index)
        merged = pd.concat([merged, stats_frame], axis=1)
        merged["log2FC"] = merged["mean_A"] - merged["mean_B"]
        merged["p.adj"] = adjust_pvalues(merged["p.value"], method=adjust)
        merged["sig"] = np.where(
            (merged["p.adj"] <= alpha) & (merged["log2FC"].abs() >= lfc_thresh),
            "yes",
            "no",
        )

        merged = merged.rename(
            columns={
                "A_vals": f"vals.{sel_a_raw}",
                "n_A": f"n.{sel_a_raw}",
                "mean_A": f"mean.{sel_a_raw}",
                "B_vals": f"vals.{sel_b_raw}",
                "n_B": f"n.{sel_b_raw}",
                "mean_B": f"mean.{sel_b_raw}",
            }
        )
        merged = merged.sort_values(by="p.adj", kind="stable").reset_index(drop=True)
        results[f"{sel_a_raw}_vs_{sel_b_raw}"] = merged

    return results
