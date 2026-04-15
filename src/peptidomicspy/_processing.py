from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from ._constants import HYDROPATHY, PEPTIDES_SELECT_COL_BASIC
from ._result import PeptidomicsResult
from ._utils import normalize_column_names, read_input_table, to_ordered_category


def calculate_gravy(peptide: str) -> float:
    if not isinstance(peptide, str) or not peptide:
        return float("nan")
    values = [HYDROPATHY.get(residue) for residue in peptide]
    if any(value is None for value in values):
        return float("nan")
    return float(sum(values) / len(values))


def _apply_group_categories(
    frame: pd.DataFrame,
    grp_cols: list[str],
    dt_int_col: pd.DataFrame,
) -> None:
    for col in grp_cols:
        levels = dt_int_col[col].drop_duplicates().tolist()
        to_ordered_category(frame, col, levels)


def _derive_type_numbers(
    dt_mean: pd.DataFrame,
    dt_reps: pd.DataFrame,
    grp_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    dt_typenum = (
        dt_mean.groupby(["Length", "Protein.name", "Protein.group", *grp_cols], sort=False, observed=True)
        .size()
        .reset_index(name="Mean.Peptides.type.number")
    )
    dt_typenum_reps = (
        dt_reps.groupby(
            ["Length", "Protein.name", "Protein.group", "Replicate", *grp_cols],
            sort=False,
            observed=True,
        )
        .size()
        .reset_index(name="Peptides.type.number")
    )
    return dt_typenum, dt_typenum_reps


def processPeptides(
    peptides_file: str | Path | pd.DataFrame,
    intensity_columns_file: str | Path | pd.DataFrame,
    protein_mapping_file: str | Path | pd.DataFrame,
) -> PeptidomicsResult:
    dt_peptides_raw = normalize_column_names(read_input_table(peptides_file))
    dt_int_col = normalize_column_names(read_input_table(intensity_columns_file))
    dt_p_map = normalize_column_names(read_input_table(protein_mapping_file))

    if "Intensity.column" not in dt_int_col.columns:
        raise KeyError("intensity_columns_file must contain 'Intensity.column'")
    dt_int_col["Intensity.column"] = dt_int_col["Intensity.column"].astype(str).str.replace(" ", ".", regex=False)

    intensity_columns = dt_int_col["Intensity.column"].tolist()
    missing_intensity_columns = [col for col in intensity_columns if col not in dt_peptides_raw.columns]
    if missing_intensity_columns:
        missing = ", ".join(missing_intensity_columns)
        raise KeyError(f"Missing intensity columns in peptides_file: {missing}")

    required_basic_columns = [col for col in PEPTIDES_SELECT_COL_BASIC if col not in dt_peptides_raw.columns]
    if required_basic_columns:
        missing = ", ".join(required_basic_columns)
        raise KeyError(f"Missing required peptide columns: {missing}")

    peptides_select_col = [*PEPTIDES_SELECT_COL_BASIC, *intensity_columns]
    dt_peptides = dt_peptides_raw.loc[:, peptides_select_col].copy()

    mask = ~dt_peptides["Leading.razor.protein"].astype(str).str.contains(
        r"^(?:CON_|REV_)",
        regex=True,
        na=False,
    )
    dt_peptides = dt_peptides.loc[mask].copy()

    mapping_cols = ["Leading.razor.protein", "Protein.name", "Protein.group"]
    dt_peptides = dt_peptides.merge(dt_p_map.loc[:, mapping_cols], on="Leading.razor.protein", how="left")
    dt_peptides["Protein.name"] = dt_peptides["Protein.name"].fillna("Others")
    dt_peptides["Protein.group"] = dt_peptides["Protein.group"].fillna("Whey")

    protein_name_levels = [*dt_p_map["Protein.name"].dropna().drop_duplicates().tolist(), "Others"]
    protein_group_levels = dt_p_map["Protein.group"].dropna().drop_duplicates().tolist()
    to_ordered_category(dt_peptides, "Protein.name", protein_name_levels)
    to_ordered_category(dt_peptides, "Protein.group", protein_group_levels)

    grp_cols = [col for col in dt_int_col.columns if col not in {"Intensity.column", "Replicate"}]
    dt_mean_parts: list[pd.DataFrame] = []
    dt_rep_parts: list[pd.DataFrame] = []
    dt_ttest_parts: list[pd.DataFrame] = []

    group_iter = dt_int_col.groupby(grp_cols, sort=False, observed=True) if grp_cols else [((), dt_int_col)]
    for group_key, sub_meta in group_iter:
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        cols = [
            *PEPTIDES_SELECT_COL_BASIC,
            "Protein.name",
            "Protein.group",
            *sub_meta["Intensity.column"].tolist(),
        ]
        dt_mean = dt_peptides.loc[:, cols].copy()
        rep_names = [f"R{i}" for i in range(1, len(sub_meta) + 1)]
        rename_map = dict(zip(sub_meta["Intensity.column"].tolist(), rep_names))
        dt_mean = dt_mean.rename(columns=rename_map)
        dt_ttest = dt_mean.copy()

        rep_values = dt_mean.loc[:, rep_names].fillna(0)
        keep_rows = (rep_values > 0).sum(axis=1) >= math.ceil(len(rep_names) / 2)
        dt_mean = dt_mean.loc[keep_rows].copy()

        rep_values = dt_mean.loc[:, rep_names].fillna(0)
        nonzero_counts = (rep_values > 0).sum(axis=1)
        dt_mean["Mean.Intensity"] = rep_values.sum(axis=1) / nonzero_counts
        dt_mean["log10.Mean.Intensity"] = np.log10(dt_mean["Mean.Intensity"])

        dt_reps = dt_mean.melt(
            id_vars=[*PEPTIDES_SELECT_COL_BASIC, "Protein.name", "Protein.group"],
            value_vars=rep_names,
            var_name="Replicate",
            value_name="Intensity",
        )
        dt_reps = dt_reps.loc[dt_reps["Intensity"].notna() & (dt_reps["Intensity"] > 0)].copy()
        dt_reps["log10.Intensity"] = np.log10(dt_reps["Intensity"])

        for col, value in zip(grp_cols, group_key):
            dt_mean[col] = value
            dt_reps[col] = value
            dt_ttest[col] = value

        dt_mean_parts.append(dt_mean)
        dt_rep_parts.append(dt_reps)
        dt_ttest_parts.append(dt_ttest)

    dt_peptides_int = pd.concat(dt_mean_parts, ignore_index=True)
    dt_peptides_int_reps = pd.concat(dt_rep_parts, ignore_index=True)
    dt_peptides_int_ttest = pd.concat(dt_ttest_parts, ignore_index=True)

    _apply_group_categories(dt_peptides_int, grp_cols, dt_int_col)
    _apply_group_categories(dt_peptides_int_reps, grp_cols, dt_int_col)
    _apply_group_categories(dt_peptides_int_ttest, grp_cols, dt_int_col)

    dt_peptides_typenum, dt_peptides_typenum_reps = _derive_type_numbers(
        dt_peptides_int,
        dt_peptides_int_reps,
        grp_cols,
    )
    _apply_group_categories(dt_peptides_typenum, grp_cols, dt_int_col)
    _apply_group_categories(dt_peptides_typenum_reps, grp_cols, dt_int_col)

    dt_peptides["GRAVY.score"] = dt_peptides["Sequence"].map(calculate_gravy)
    dt_peptides_int["GRAVY.score"] = dt_peptides_int["Sequence"].map(calculate_gravy)
    dt_peptides_int_reps["GRAVY.score"] = dt_peptides_int_reps["Sequence"].map(calculate_gravy)

    return PeptidomicsResult(
        dt_peptides=dt_peptides,
        dt_peptides_int=dt_peptides_int,
        dt_peptides_int_reps=dt_peptides_int_reps,
        dt_peptides_int_ttest=dt_peptides_int_ttest,
        dt_peptides_typenum=dt_peptides_typenum,
        dt_peptides_typenum_reps=dt_peptides_typenum_reps,
        dt_int_col=dt_int_col.copy(deep=True),
        grp_cols=grp_cols,
        peptides_select_col_basic=list(PEPTIDES_SELECT_COL_BASIC),
    )


def filterPeptides(
    result: PeptidomicsResult,
    seqs: list[str] | tuple[str, ...] | None = None,
    seq_pattern: str | None = None,
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
) -> PeptidomicsResult:
    dt_p = result.dt_peptides.copy(deep=True)
    dt_mean = result.dt_peptides_int.copy(deep=True)
    dt_mean_r = result.dt_peptides_int_reps.copy(deep=True)
    dt_ttest = result.dt_peptides_int_ttest.copy(deep=True)
    dt_int_col = result.dt_int_col.copy(deep=True)
    grp_cols = list(result.grp_cols)

    keep_seqs: list[str] = []
    all_seqs = pd.unique(dt_p["Sequence"]).tolist()

    if seqs is not None:
        keep_seqs = list(dict.fromkeys([*keep_seqs, *list(seqs)]))
    if seq_pattern is not None:
        matched = [seq for seq in all_seqs if re.search(seq_pattern, str(seq))]
        keep_seqs = list(dict.fromkeys([*keep_seqs, *matched]))
    if len(keep_seqs) == 0 and filter_params is None:
        keep_seqs = all_seqs

    dt_p = dt_p.loc[dt_p["Sequence"].isin(keep_seqs)].copy()
    dt_mean = dt_mean.loc[dt_mean["Sequence"].isin(keep_seqs)].copy()
    dt_mean_r = dt_mean_r.loc[dt_mean_r["Sequence"].isin(keep_seqs)].copy()
    dt_ttest = dt_ttest.loc[dt_ttest["Sequence"].isin(keep_seqs)].copy()

    dt_typenum, dt_typenum_r = _derive_type_numbers(dt_mean, dt_mean_r, grp_cols)

    if filter_params is not None:
        for col, values in filter_params.items():
            allowed = values if isinstance(values, (list, tuple, set, pd.Series, np.ndarray)) else [values]
            dt_mean = dt_mean.loc[dt_mean[col].isin(allowed)].copy()
            dt_mean_r = dt_mean_r.loc[dt_mean_r[col].isin(allowed)].copy()
            dt_ttest = dt_ttest.loc[dt_ttest[col].isin(allowed)].copy()
            dt_typenum = dt_typenum.loc[dt_typenum[col].isin(allowed)].copy()
            dt_typenum_r = dt_typenum_r.loc[dt_typenum_r[col].isin(allowed)].copy()
            dt_int_col = dt_int_col.loc[dt_int_col[col].isin(allowed)].copy()

        columns = [
            *result.peptides_select_col_basic,
            *dt_int_col["Intensity.column"].tolist(),
            "Protein.name",
            "Protein.group",
        ]
        dt_p = dt_p.loc[:, columns].copy()

    _apply_group_categories(dt_mean, grp_cols, dt_int_col)
    _apply_group_categories(dt_mean_r, grp_cols, dt_int_col)
    _apply_group_categories(dt_typenum, grp_cols, dt_int_col)
    _apply_group_categories(dt_typenum_r, grp_cols, dt_int_col)

    return PeptidomicsResult(
        dt_peptides=dt_p,
        dt_peptides_int=dt_mean,
        dt_peptides_int_reps=dt_mean_r,
        dt_peptides_int_ttest=dt_ttest,
        dt_peptides_typenum=dt_typenum,
        dt_peptides_typenum_reps=dt_typenum_r,
        dt_int_col=dt_int_col,
        grp_cols=grp_cols,
        peptides_select_col_basic=list(result.peptides_select_col_basic),
    )
