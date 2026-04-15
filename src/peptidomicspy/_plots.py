from __future__ import annotations

import math
import importlib
import warnings
from pathlib import Path

import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from logomaker.src.colors import get_color_dict
from matplotlib import colormaps
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
from mizani.palettes import hue_pal
from mizani.breaks import breaks_extended
from plotnine import (
    aes,
    facet_grid,
    geom_col,
    geom_density,
    geom_hline,
    geom_point,
    geom_text,
    geom_vline,
    ggplot,
    labs,
    scale_color_manual,
    scale_fill_manual,
    scale_y_continuous,
    theme,
    theme_bw,
)
from plotnine.themes.elements import element_blank, element_line, element_rect, element_text

from ._constants import AMINO_ACIDS, DEFAULT_COLOR, PROTEIN_COLOR
from ._datasets import load_example_protein_mapping
from ._result import PeptidomicsResult
from ._utils import (
    interaction_label,
    scientific_10_label,
    scientific_10_labels,
    split_facet_expression,
)


plotnine_stat_density = importlib.import_module("plotnine.stats.stat_density")
_ORIGINAL_PLOTNINE_COMPUTE_DENSITY = plotnine_stat_density.compute_density


def _compute_density_with_writable_inputs(x, weight, range_x, params):
    x_array = np.asarray(x, dtype=float).copy()
    weight_array = None if weight is None else np.asarray(weight, dtype=float).copy()
    return _ORIGINAL_PLOTNINE_COMPUTE_DENSITY(x_array, weight_array, range_x, params)


if plotnine_stat_density.compute_density is not _compute_density_with_writable_inputs:
    plotnine_stat_density.compute_density = _compute_density_with_writable_inputs


def _apply_filter_params(
    frame: pd.DataFrame,
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None,
) -> pd.DataFrame:
    out = frame
    if not filter_params:
        return out
    for col, values in filter_params.items():
        allowed = values if isinstance(values, (list, tuple, set, pd.Series, np.ndarray)) else [values]
        out = out.loc[out[col].isin(allowed)].copy()
    return out


def _palette_for_levels(levels: list[str]) -> dict[str, str]:
    palette = hue_pal()(max(1, len(levels)))
    return {level: palette[idx] for idx, level in enumerate(levels)}


def _theme_pubr() -> theme:
    return theme_bw() + theme(
        panel_border=element_blank(),
        panel_grid_major=element_blank(),
        panel_grid_minor=element_blank(),
        axis_line=element_line(color="black", size=0.5),
        axis_text=element_text(color="black"),
        legend_key=element_blank(),
        strip_background=element_rect(fill="#F2F2F2", color="black", size=0.7),
        legend_position="top",
    )


def _facet_plot(plot: ggplot, facet_rows: str | None, facet_cols: str | None) -> ggplot:
    if facet_rows is None and facet_cols is None:
        return plot
    rows = facet_rows or "."
    cols = facet_cols or "."
    return plot + facet_grid(f"{rows} ~ {cols}")


def _warn_missing_distinctions(
    frame: pd.DataFrame,
    varying_vars: list[str],
    used_vars: list[str],
    context: str = "facets",
) -> None:
    varying = [var for var in varying_vars if frame[var].nunique(dropna=False) > 1]
    missing = [var for var in varying if var not in used_vars]
    if not missing:
        return
    if context == "x_var/facets":
        message = "Grouping variable(s) {} not used in x_var/facets; data may be aggregated across them."
    else:
        message = "Grouping variable(s) {} not used in facets; data may be aggregated across them."
    warnings.warn(message.format(", ".join(missing)), stacklevel=2)


def _prepare_sample_groups(
    frame: pd.DataFrame,
    sample_vars: list[str],
) -> tuple[pd.DataFrame, str, dict[str, str]]:
    out = frame.copy(deep=True)
    out["Sample"] = interaction_label(out, sample_vars, sep="_")
    legend_title = " × ".join(sample_vars) if sample_vars else "Sample"
    levels = sorted(out["Sample"].astype(str).unique().tolist())
    out["Sample"] = pd.Categorical(out["Sample"], categories=levels, ordered=True)
    return out, legend_title, _palette_for_levels(levels)


def _scientific_formatter():
    return FuncFormatter(lambda value, _position: scientific_10_label(value))


def _ordered_interaction_levels(frame: pd.DataFrame, columns: list[str], sep: str) -> list[str]:
    if not columns:
        return ["Sample"]
    unique = frame.loc[:, columns].drop_duplicates().copy()
    unique = unique.sort_values(columns, kind="mergesort")
    return unique.astype(str).agg(sep.join, axis=1).tolist()


def plot_int(
    result: PeptidomicsResult,
    type: str = "mean",
    x_var: str | None = None,
    color_by: str = "Protein.group",
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
    facet_rows: str | None = None,
    facet_cols: str | None = None,
    scientific_10_y: bool = True,
):
    if type not in {"mean", "reps"}:
        raise ValueError("type must be either 'mean' or 'reps'.")
    if color_by not in {"Protein.group", "Protein.name", "none"}:
        raise ValueError("color_by must be 'Protein.group', 'Protein.name', or 'none'.")

    frame = result.dt_peptides_int_reps.copy(deep=True) if type == "reps" else result.dt_peptides_int.copy(deep=True)
    frame = _apply_filter_params(frame, filter_params)
    grp_cols = list(result.grp_cols)
    if x_var is None:
        x_var = grp_cols[0]
    valid_x = [*grp_cols, "Replicate"] if type == "reps" else grp_cols
    if x_var not in valid_x:
        raise ValueError(f"x_var must be one of {valid_x}.")
    if (
        type == "reps"
        and x_var != "Replicate"
        and facet_cols is None
        and (facet_rows is None or "Replicate" not in facet_rows)
    ):
        facet_cols = "Replicate"

    used_vars = [x_var, *split_facet_expression(facet_rows), *split_facet_expression(facet_cols)]
    _warn_missing_distinctions(frame, [*grp_cols, *([] if type == "mean" else ["Replicate"])], used_vars, "x_var/facets")

    y_var = "Intensity" if type == "reps" else "Mean.Intensity"
    if color_by == "none":
        plot = ggplot(frame, aes(x=x_var, y=y_var)) + geom_col(fill=DEFAULT_COLOR)
    else:
        plot = ggplot(frame, aes(x=x_var, y=y_var, fill=color_by)) + geom_col()
        plot = plot + scale_fill_manual(values=PROTEIN_COLOR)

    plot = plot + _theme_pubr()
    if scientific_10_y:
        plot = plot + scale_y_continuous(labels=scientific_10_labels)
    return _facet_plot(plot, facet_rows, facet_cols)


def _plot_type_num_impl(
    result: PeptidomicsResult,
    type: str = "mean",
    x_var: str | None = None,
    color_by: str = "Protein.group",
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
    facet_rows: str | None = None,
    facet_cols: str | None = None,
    scientific_10_y: bool = False,
):
    if type not in {"mean", "reps"}:
        raise ValueError("type must be either 'mean' or 'reps'.")
    if color_by not in {"Protein.group", "Protein.name", "none"}:
        raise ValueError("color_by must be 'Protein.group', 'Protein.name', or 'none'.")

    frame = (
        result.dt_peptides_typenum_reps.copy(deep=True)
        if type == "reps"
        else result.dt_peptides_typenum.copy(deep=True)
    )
    frame = _apply_filter_params(frame, filter_params)
    grp_cols = list(result.grp_cols)
    if x_var is None:
        x_var = grp_cols[0]
    valid_x = [*grp_cols, "Replicate"] if type == "reps" else grp_cols
    if x_var not in valid_x:
        raise ValueError(f"x_var must be one of {valid_x}.")
    if (
        type == "reps"
        and x_var != "Replicate"
        and facet_cols is None
        and (facet_rows is None or "Replicate" not in facet_rows)
    ):
        facet_cols = "Replicate"

    used_vars = [x_var, *split_facet_expression(facet_rows), *split_facet_expression(facet_cols)]
    _warn_missing_distinctions(frame, [*grp_cols, *([] if type == "mean" else ["Replicate"])], used_vars, "x_var/facets")

    y_var = "Peptides.type.number" if type == "reps" else "Mean.Peptides.type.number"
    if color_by == "none":
        plot = ggplot(frame, aes(x=x_var, y=y_var)) + geom_col(fill=DEFAULT_COLOR)
    else:
        plot = ggplot(frame, aes(x=x_var, y=y_var, fill=color_by)) + geom_col()
        plot = plot + scale_fill_manual(values=PROTEIN_COLOR)

    plot = plot + _theme_pubr()
    if scientific_10_y:
        plot = plot + scale_y_continuous(labels=scientific_10_labels)
    return _facet_plot(plot, facet_rows, facet_cols)


def plot_type_num(
    result: PeptidomicsResult,
    type: str = "mean",
    x_var: str | None = None,
    color_by: str = "Protein.group",
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
    facet_rows: str | None = None,
    facet_cols: str | None = None,
    scientific_10_y: bool = False,
):
    return _plot_type_num_impl(
        result=result,
        type=type,
        x_var=x_var,
        color_by=color_by,
        filter_params=filter_params,
        facet_rows=facet_rows,
        facet_cols=facet_cols,
        scientific_10_y=scientific_10_y,
    )


def plot_count(
    result: PeptidomicsResult,
    type: str = "mean",
    x_var: str | None = None,
    color_by: str = "Protein.group",
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
    facet_rows: str | None = None,
    facet_cols: str | None = None,
    scientific_10_y: bool = False,
):
    warnings.warn(
        "plot_count() is deprecated and has been replaced by plot_type_num(). Please use plot_type_num() instead.",
        stacklevel=2,
    )
    return _plot_type_num_impl(
        result=result,
        type=type,
        x_var=x_var,
        color_by=color_by,
        filter_params=filter_params,
        facet_rows=facet_rows,
        facet_cols=facet_cols,
        scientific_10_y=scientific_10_y,
    )


def plot_length_distribution(
    result: PeptidomicsResult,
    type: str = "mean",
    metric: str = "intensity",
    color_by: str = "Protein.group",
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
    facet_rows: str | None = None,
    facet_cols: str | None = None,
    scientific_10_y: bool = True,
    plot_mode: str = "bar",
):
    if type not in {"mean", "reps"}:
        raise ValueError("type must be either 'mean' or 'reps'.")
    if metric == "count":
        metric = "type_num"
    if metric not in {"intensity", "type_num"}:
        raise ValueError("metric must be either 'intensity' or 'type_num'.")
    if plot_mode not in {"bar", "density"}:
        raise ValueError("plot_mode must be either 'bar' or 'density'.")
    if plot_mode == "bar" and color_by not in {"Protein.group", "Protein.name", "none"}:
        raise ValueError("color_by must be 'Protein.group', 'Protein.name', or 'none'.")

    if type == "mean":
        frame = result.dt_peptides_int.copy(deep=True) if metric == "intensity" else result.dt_peptides_typenum.copy(deep=True)
        y_var = "Mean.Intensity" if metric == "intensity" else "Mean.Peptides.type.number"
    else:
        frame = result.dt_peptides_int_reps.copy(deep=True) if metric == "intensity" else result.dt_peptides_typenum_reps.copy(deep=True)
        y_var = "Intensity" if metric == "intensity" else "Peptides.type.number"
    frame = _apply_filter_params(frame, filter_params)
    grp_cols = list(result.grp_cols)

    if plot_mode == "bar":
        if type == "reps" and facet_cols is None and (facet_rows is None or "Replicate" not in facet_rows):
            facet_cols = "Replicate"
        if color_by == "none":
            plot = ggplot(frame, aes(x="Length", y=y_var)) + geom_col(fill=DEFAULT_COLOR)
        else:
            plot = ggplot(frame, aes(x="Length", y=y_var, fill=color_by)) + geom_col()
            plot = plot + scale_fill_manual(values=PROTEIN_COLOR)
        plot = plot + _theme_pubr()
        if scientific_10_y:
            plot = plot + scale_y_continuous(labels=scientific_10_labels)
        facet_vars = [*split_facet_expression(facet_rows), *split_facet_expression(facet_cols)]
        _warn_missing_distinctions(frame, [*grp_cols, *([] if type == "mean" else ["Replicate"])], facet_vars, "facets")
        return _facet_plot(plot, facet_rows, facet_cols)

    sample_vars = [*grp_cols, *([] if type == "mean" else ["Replicate"])]
    sample_vars = [var for var in sample_vars if frame[var].nunique(dropna=False) > 1]
    facet_vars = [*split_facet_expression(facet_rows), *split_facet_expression(facet_cols)]
    sample_vars = [var for var in sample_vars if var not in facet_vars]
    frame, legend_title, sample_colors = _prepare_sample_groups(frame, sample_vars)
    frame["Length"] = np.asarray(frame["Length"], dtype=float).copy()
    frame[y_var] = np.asarray(frame[y_var], dtype=float).copy()
    frame = pd.DataFrame(frame).reset_index(drop=True).copy(deep=True)

    plot = (
        ggplot(frame, aes(x="Length", color="Sample"))
        + geom_density(aes(weight=y_var), size=1)
        + scale_color_manual(values=sample_colors)
        + labs(color=legend_title, y="Density")
        + _theme_pubr()
    )
    return _facet_plot(plot, facet_rows, facet_cols)


def plot_gravy_vs_intensity(
    result: PeptidomicsResult,
    type: str = "mean",
    color_by: str = "Protein.group",
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
    facet_rows: str | None = None,
    facet_cols: str | None = None,
    alpha_value: float = 0.8,
    plot_mode: str = "scatter",
):
    if type not in {"mean", "reps"}:
        raise ValueError("type must be either 'mean' or 'reps'.")
    if plot_mode not in {"scatter", "density"}:
        raise ValueError("plot_mode must be either 'scatter' or 'density'.")
    if plot_mode == "scatter" and color_by not in {"Protein.group", "Protein.name", "none"}:
        raise ValueError("color_by must be 'Protein.group', 'Protein.name', or 'none'.")

    if type == "mean":
        frame = result.dt_peptides_int.copy(deep=True)
        y_var = "log10.Mean.Intensity"
        weight_var = "Mean.Intensity"
    else:
        frame = result.dt_peptides_int_reps.copy(deep=True)
        y_var = "log10.Intensity"
        weight_var = "Intensity"
    frame = _apply_filter_params(frame, filter_params)
    grp_cols = list(result.grp_cols)

    facet_vars = [*split_facet_expression(facet_rows), *split_facet_expression(facet_cols)]
    if plot_mode == "scatter":
        if type == "reps" and facet_cols is None and (facet_rows is None or "Replicate" not in facet_rows):
            facet_cols = "Replicate"
            facet_vars = [*split_facet_expression(facet_rows), *split_facet_expression(facet_cols)]
        if color_by == "none":
            plot = ggplot(frame, aes(x="GRAVY.score", y=y_var)) + geom_point(color=DEFAULT_COLOR, alpha=alpha_value)
        else:
            plot = ggplot(frame, aes(x="GRAVY.score", y=y_var, color=color_by)) + geom_point(alpha=alpha_value)
            plot = plot + scale_color_manual(values=PROTEIN_COLOR)
        plot = plot + _theme_pubr()
        _warn_missing_distinctions(frame, [*grp_cols, *([] if type == "mean" else ["Replicate"])], facet_vars, "facets")
        return _facet_plot(plot, facet_rows, facet_cols)

    frame = frame.loc[frame["GRAVY.score"].notna() & frame[weight_var].notna()].copy()
    sample_vars = [*grp_cols, *([] if type == "mean" else ["Replicate"])]
    sample_vars = [var for var in sample_vars if frame[var].nunique(dropna=False) > 1]
    sample_vars = [var for var in sample_vars if var not in facet_vars]
    frame, legend_title, sample_colors = _prepare_sample_groups(frame, sample_vars)
    frame["GRAVY.score"] = np.asarray(frame["GRAVY.score"], dtype=float).copy()
    frame[weight_var] = np.asarray(frame[weight_var], dtype=float).copy()
    frame = pd.DataFrame(frame).reset_index(drop=True).copy(deep=True)

    plot = (
        ggplot(frame, aes(x="GRAVY.score", color="Sample"))
        + geom_density(aes(weight=weight_var), size=1)
        + scale_color_manual(values=sample_colors)
        + labs(color=legend_title, y="Density")
        + _theme_pubr()
    )
    used_vars = [*facet_vars, *sample_vars]
    _warn_missing_distinctions(frame, [*grp_cols, *([] if type == "mean" else ["Replicate"])], used_vars, "facets")
    return _facet_plot(plot, facet_rows, facet_cols)


def plot_volcano(
    ttest_result: dict[str, pd.DataFrame],
    comparisons,
    test_method: str = "treat",
    show_threshold: bool = True,
    lfc_thresh: float = 1,
    alpha: float = 0.05,
    fill_values: dict[str, str] | None = None,
    point_size: float = 2,
    point_alpha: float = 0.85,
    label_seqs: list[str] | None = None,
    label_size: float = 3,
    label_col: str = "black",
    highlight_seqs: dict[str, list[str]] | None = None,
    highlight_size: float | None = None,
    highlight_stroke: float = 1.2,
    y_log_scale: bool = False,
):
    fill_values = fill_values or {"no": "#BFBFBF", "yes": "#FFC010"}
    selected_keys: list[str] = []
    if isinstance(comparisons, int):
        selected_keys = [list(ttest_result)[comparisons - 1]]
    elif isinstance(comparisons, (list, tuple)) and comparisons and all(isinstance(item, int) for item in comparisons):
        selected_keys = [list(ttest_result)[idx - 1] for idx in comparisons]
    elif isinstance(comparisons, (list, tuple)) and len(comparisons) == 2 and all(
        isinstance(item, str) for item in comparisons
    ):
        selected_keys = [f"{comparisons[0]}_vs_{comparisons[1]}"]
    elif isinstance(comparisons, list):
        for item in comparisons:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                selected_keys.append(f"{item[0]}_vs_{item[1]}")
            else:
                raise ValueError("comparisons must be indexes or name pairs.")
    else:
        raise ValueError("comparisons must be indexes or name pairs.")

    out: dict[str, ggplot] = {}
    for key in selected_keys:
        frame = ttest_result[key].copy(deep=True)
        padj = frame["p.adj"].astype(float).copy()
        if (padj == 0).any():
            positive = padj.loc[padj > 0]
            min_pos = float(positive.min()) if not positive.empty else np.finfo(float).tiny
            padj = np.where(padj == 0, min_pos * 0.1, padj)
        frame["yvar"] = -np.log10(np.asarray(padj, dtype=float))
        frame["sig"] = pd.Categorical(frame["sig"], categories=["yes", "no"], ordered=True)

        if y_log_scale:
            positive_y = frame.loc[frame["yvar"] > 0, "yvar"]
            y_min_pos = float(positive_y.min()) if not positive_y.empty else 1e-3
            eps = y_min_pos * 0.1
            frame["yvar_plot"] = np.where(frame["yvar"] > 0, frame["yvar"], eps)
            y_range = frame["yvar_plot"].replace([np.inf, -np.inf], np.nan).dropna()
            e_min = int(np.floor(np.log10(float(y_range.min()))))
            e_max = int(np.ceil(np.log10(float(y_range.max()))))
            y_breaks = [10**exponent for exponent in range(e_min, e_max + 1)]
            y_scale = scale_y_continuous(
                trans="log10",
                breaks=y_breaks,
                labels=lambda values: [f"{float(value):g}" for value in values],
            )
        else:
            frame["yvar_plot"] = frame["yvar"]
            y_scale = scale_y_continuous()

        plot = (
            ggplot(frame, aes(x="log2FC", y="yvar_plot"))
            + geom_point(aes(fill="sig", color="sig"), shape="o", size=point_size, alpha=point_alpha, stroke=0.6)
            + scale_fill_manual(values=fill_values, name="Significant")
            + scale_color_manual(values=fill_values, name="Significant")
            + labs(x="log2 fold change", y=r"$-\log_{10}$(adjusted p-value)", title=key)
            + y_scale
            + _theme_pubr()
            + theme(legend_position="bottom")
        )

        if show_threshold:
            plot = plot + geom_hline(yintercept=-math.log10(alpha), linetype="dashed", size=0.5)
            if test_method == "plain":
                plot = plot + geom_vline(xintercept=[-lfc_thresh, lfc_thresh], linetype="dashed", size=0.5)

        if label_seqs:
            labeled = frame.loc[frame["Sequence"].isin(label_seqs)].copy()
            if not labeled.empty:
                plot = plot + geom_text(
                    data=labeled,
                    mapping=aes(label="Sequence"),
                    color=label_col,
                    size=label_size * 2.4,
                    adjust_text={
                        "expand": (1.15, 1.35),
                        "force_text": (0.4, 0.8),
                        "force_static": (0.2, 0.5),
                        "force_pull": (0.02, 0.05),
                        "iter_lim": 300,
                        "min_arrow_len": 0,
                        "arrowprops": {"arrowstyle": "-", "color": label_col, "lw": 0.4},
                    },
                )
                plot = plot + geom_point(
                    data=labeled,
                    mapping=aes(x="log2FC", y="yvar_plot"),
                    shape="o",
                    size=point_size,
                    stroke=0.9,
                    color=label_col,
                    fill=None,
                )

        if highlight_seqs:
            ring_size = highlight_size if highlight_size is not None else point_size + 0.6
            for color, seqs in highlight_seqs.items():
                highlighted = frame.loc[frame["Sequence"].isin(seqs)].copy()
                if highlighted.empty:
                    continue
                plot = plot + geom_point(
                    data=highlighted,
                    mapping=aes(x="log2FC", y="yvar_plot"),
                    shape="o",
                    size=ring_size,
                    stroke=highlight_stroke,
                    color=color,
                    fill=None,
                )

        out[key] = plot
    return out


def _stack_alignment_rows(sample_frame: pd.DataFrame, start_col: str, end_col: str) -> pd.DataFrame:
    sample_frame = sample_frame.sort_values([start_col, end_col], ascending=[True, False]).copy()
    if sample_frame.empty:
        sample_frame["stack_row"] = []
        return sample_frame

    stack_group = [1]
    current_group = 1
    current_end = float(sample_frame.iloc[0][end_col])
    for _, peptide in sample_frame.iloc[1:].iterrows():
        start_value = float(peptide[start_col])
        end_value = float(peptide[end_col])
        if start_value <= current_end + 1:
            stack_group.append(current_group)
            current_end = max(current_end, end_value)
        else:
            current_group += 1
            stack_group.append(current_group)
            current_end = end_value
    sample_frame["stack_group"] = stack_group

    sample_frame["stack_row"] = np.nan
    for group_id in sorted(sample_frame["stack_group"].unique().tolist()):
        group = sample_frame.loc[sample_frame["stack_group"] == group_id].sort_values(start_col).copy()
        row_assignments: list[int] = []
        row_ends: list[float] = []
        for _, peptide in group.iterrows():
            start_value = float(peptide[start_col])
            end_value = float(peptide[end_col])
            eligible = [idx for idx, row_end in enumerate(row_ends) if start_value > row_end + 1]
            if eligible:
                chosen = max(eligible, key=lambda idx: row_ends[idx])
                row_assignments.append(chosen + 1)
                row_ends[chosen] = end_value
            else:
                row_assignments.append(len(row_ends) + 1)
                row_ends.append(end_value)
        sample_frame.loc[group.index, "stack_row"] = row_assignments

    sample_frame["stack_row"] = sample_frame["stack_row"].astype(int)
    return sample_frame


def plot_pep_align(
    result: PeptidomicsResult,
    protein_name: str,
    protein_seq: str | None = None,
    type: str = "mean",
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
    x_interval: float = 25,
    seq_col: str = "Sequence",
    start_col: str = "Start.position",
    end_col: str = "End.position",
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    intensity_col: str | None = None,
    label_seq: dict[str, list[str]] | None = None,
    label_col: dict[str, str] | None = None,
    auto_size: bool = True,
    save_file_location: str | None = None,
    save_file_dpi: float = 300,
    width_per_aa: float = 0.1,
    height_per_row: float = 0.18,
    min_width: float = 4,
    min_height: float = 3,
):
    if type not in {"mean", "reps"}:
        raise ValueError("type must be one of 'mean' or 'reps'.")
    frame = result.dt_peptides_int.copy(deep=True) if type == "mean" else result.dt_peptides_int_reps.copy(deep=True)
    grp_cols = list(result.grp_cols)
    if intensity_col is None:
        intensity_col = "Mean.Intensity" if type == "mean" else "Intensity"

    if filter_params:
        allowed = grp_cols if type == "mean" else [*grp_cols, "Replicate"]
        invalid = [col for col in filter_params if col not in allowed]
        if invalid:
            raise ValueError(f"Filter column(s) not found in the result: {', '.join(invalid)}")

    frame = _apply_filter_params(frame, filter_params)
    required = {"Leading.razor.protein", seq_col, start_col, end_col, intensity_col}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(sorted(missing))}")

    frame = frame.loc[frame["Leading.razor.protein"] == protein_name].copy()
    frame = frame.loc[frame[start_col].notna() & frame[end_col].notna()].copy()
    if frame.empty:
        raise ValueError("No peptides found for protein_name after filtering.")

    if protein_seq is None:
        mapping = load_example_protein_mapping()
        hit = mapping.loc[mapping["Leading.razor.protein"] == protein_name, "Protein.seq"]
        if not hit.empty and pd.notna(hit.iloc[0]):
            protein_seq = str(hit.iloc[0])
        else:
            warnings.warn(
                f"Protein {protein_name} not found in peptidomics_protein_mapping_example; using peptide endpoints for axis scaling.",
                stacklevel=2,
            )
    protein_length = len(protein_seq) if protein_seq else int(frame[end_col].max())
    if protein_length <= 0:
        raise ValueError("Unable to determine protein length; please provide protein_seq.")

    facet_vars = [*grp_cols, *([] if type == "mean" else ["Replicate"])]
    if filter_params:
        facet_vars = [var for var in facet_vars if var not in filter_params]
    facet_vars = [var for var in facet_vars if frame[var].nunique(dropna=False) > 1]
    if facet_vars:
        frame["sample_panel"] = interaction_label(frame, facet_vars, sep=" | ")
        panel_levels = _ordered_interaction_levels(frame, facet_vars, sep=" | ")
    else:
        frame["sample_panel"] = "Sample"
        panel_levels = ["Sample"]
    frame["sample_panel"] = pd.Categorical(frame["sample_panel"], categories=panel_levels, ordered=True)

    frame[start_col] = pd.to_numeric(frame[start_col])
    frame[end_col] = pd.to_numeric(frame[end_col])
    frame[intensity_col] = pd.to_numeric(frame[intensity_col])
    frame = frame.loc[frame[intensity_col].notna() & (frame[intensity_col] > 0)].copy()
    if frame.empty:
        raise ValueError("No peptides found for protein_name after dropping non-positive intensities.")

    stacked = pd.concat(
        [_stack_alignment_rows(group.copy(), start_col, end_col) for _, group in frame.groupby("sample_panel", sort=False)],
        ignore_index=True,
    )
    stacked["rect_xmax"] = stacked[end_col] + 1

    label_seq = label_seq or {}
    label_col = label_col or {}
    default_label_colors = ["black", "#4daf4a", "#984ea3", "#a65628", "#f781bf"]
    label_color_map: dict[str, str] = {}
    for idx, label_name in enumerate(label_seq):
        label_color_map[label_name] = label_col.get(label_name, default_label_colors[idx % len(default_label_colors)])

    panel_data = [stacked.loc[stacked["sample_panel"] == panel].copy() for panel in panel_levels]
    panel_heights = [
        max(int(panel["stack_row"].max()), 1) if y_range is None else int(abs(y_range[0] - y_range[1]) + 1)
        for panel in panel_data
    ]
    total_rows = float(sum(panel_heights))

    if x_range is None:
        x_range = (1.0, float(protein_length + 1))
    if auto_size:
        visible_span = abs(float(x_range[1]) - float(x_range[0]))
        width = max(min_width, visible_span * width_per_aa)
        height = max(min_height, total_rows * height_per_row)
    else:
        width = min_width
        height = min_height

    fig = plt.figure(figsize=(width, height))
    grid = GridSpec(1 + len(panel_data), 1, figure=fig, height_ratios=[1, *panel_heights], hspace=0.03)
    protein_ax = fig.add_subplot(grid[0, 0])
    peptide_axes = [fig.add_subplot(grid[idx + 1, 0], sharex=protein_ax) for idx in range(len(panel_data))]

    for position in range(1, protein_length + 1):
        protein_ax.add_patch(
            Rectangle((position, 0.5), 1, 1, facecolor="0.95", edgecolor="0.3", linewidth=0.1)
        )
    if protein_seq:
        for position, letter in enumerate(protein_seq, start=1):
            protein_ax.text(
                float(position) + 0.5,
                1.0,
                letter,
                ha="center",
                va="center",
                fontsize=6,
                color="black",
                clip_on=True,
            )
    protein_break_end = max(float(protein_length + 1), float(x_range[1]))
    protein_breaks = [1, *range(int(x_interval), int(protein_break_end + x_interval), int(x_interval))]
    protein_ax.set_xlim(*x_range)
    protein_ax.set_ylim(0.4, 1.6)
    protein_ax.xaxis.tick_top()
    protein_ax.set_xticks(protein_breaks)
    protein_ax.set_yticks([])
    protein_ax.set_ylabel("")
    protein_ax.set_xlabel("")
    protein_ax.text(-0.04, 0.2, "P", rotation=90, transform=protein_ax.transAxes, ha="center", va="center", fontsize=10)
    protein_ax.spines["right"].set_visible(False)
    protein_ax.spines["left"].set_visible(False)
    protein_ax.spines["bottom"].set_visible(False)

    intensity_min = float(stacked[intensity_col].min())
    intensity_max = float(stacked[intensity_col].max())
    norm = LogNorm(vmin=intensity_min, vmax=intensity_max)
    cmap = colormaps.get_cmap("RdYlBu_r")

    for axis, panel_name, panel_frame in zip(peptide_axes, panel_levels, panel_data):
        for _, row in panel_frame.iterrows():
            rect = Rectangle(
                (float(row[start_col]), float(row["stack_row"]) - 0.45),
                float(row["rect_xmax"] - row[start_col]),
                0.9,
                facecolor=cmap(norm(float(row[intensity_col]))),
                edgecolor="none",
                linewidth=0,
                clip_on=True,
            )
            for label_name, seqs in label_seq.items():
                if row[seq_col] in seqs:
                    rect.set_edgecolor(label_color_map[label_name])
                    rect.set_linewidth(2)
            axis.add_patch(rect)

            sequence = str(row[seq_col])
            for offset, letter in enumerate(sequence, start=1):
                axis.text(
                    float(row[start_col]) + offset - 0.5,
                    float(row["stack_row"]),
                    letter,
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="black",
                    clip_on=True,
                )

        axis.set_xlim(*x_range)
        if y_range is None:
            ymax_use = float(panel_frame["stack_row"].max()) + 0.5
            ymin_use = -0.5
        else:
            # Match the automatic panel padding so custom row windows include
            # the full peptide rectangles rather than cutting the edge rows.
            ymax_use = float(y_range[1]) + 0.5
            ymin_use = float(y_range[0]) - 0.5
        axis.set_ylim(ymax_use, ymin_use)
        axis.set_yticks([])
        axis.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_visible(False)
        axis.text(-0.04, 0.5, str(panel_name), rotation=90, transform=axis.transAxes, ha="center", va="center", fontsize=10)

    scalar = ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    colorbar = fig.colorbar(
        scalar,
        ax=[protein_ax, *peptide_axes],
        orientation="horizontal",
        fraction=0.035,
        pad=0.055,
        shrink=0.32,
    )
    colorbar.set_label(intensity_col)
    colorbar.ax.xaxis.set_major_formatter(_scientific_formatter())

    bottom_margin = 0.16
    if label_color_map:
        handles = [
            Patch(facecolor="none", edgecolor=color, linewidth=2, label=label_name)
            for label_name, color in label_color_map.items()
        ]
        fig.legend(
            handles=handles,
            title="Label",
            loc="lower center",
            bbox_to_anchor=(0.5, 0.005),
            ncol=min(len(handles), 3),
            frameon=False,
        )
        bottom_margin = 0.24

    fig.suggested_fig_width = width
    fig.suggested_fig_height = height
    fig.suggested_fig_size = (width, height)
    fig.subplots_adjust(top=0.96, bottom=bottom_margin, hspace=0.03)
    if save_file_location:
        fig.savefig(save_file_location, dpi=save_file_dpi, bbox_inches="tight")
        fig.save_file_location = save_file_location
        fig.save_file_dpi = save_file_dpi
    return fig


def plot_cleavage_site(
    result: PeptidomicsResult,
    terminal: str = "both",
    measure: str = "intensity",
    replicate_mode: str = "mean",
    filter_params: dict[str, list[str] | tuple[str, ...] | str] | None = None,
    scientific_10_y: bool = True,
    drop_constant_groups: bool = True,
):
    if terminal not in {"both", "N", "C"}:
        raise ValueError("terminal must be one of 'both', 'N', or 'C'.")
    if measure not in {"intensity", "type_num"}:
        raise ValueError("measure must be one of 'intensity' or 'type_num'.")
    if replicate_mode not in {"mean", "reps"}:
        raise ValueError("replicate_mode must be one of 'mean' or 'reps'.")

    frame = result.dt_peptides_int_reps.copy(deep=True)
    frame = _apply_filter_params(frame, filter_params)
    grp_cols = list(result.grp_cols)
    group_cols = [*grp_cols, *([] if replicate_mode == "mean" else ["Replicate"])]
    if drop_constant_groups:
        active_group_cols = [col for col in group_cols if frame[col].nunique(dropna=False) > 1]
    else:
        active_group_cols = group_cols

    def build_logo_matrix(side: str) -> tuple[pd.DataFrame, list[str]]:
        aa_col = "First.amino.acid" if side == "N" else "Last.amino.acid"
        rep_group_cols = [aa_col, *active_group_cols]

        if measure == "intensity":
            grouped_rep = (
                frame.groupby(rep_group_cols, sort=False, observed=True)["Intensity"]
                .sum()
                .reset_index(name="value")
            )
        else:
            grouped_rep = (
                frame.loc[frame["Intensity"].notna() & (frame["Intensity"] > 0)]
                .groupby(rep_group_cols, sort=False, observed=True)["Sequence"]
                .nunique()
                .reset_index(name="value")
            )

        if replicate_mode == "mean":
            grouped = (
                grouped_rep.groupby([aa_col, *active_group_cols], sort=False, observed=True)["value"]
                .mean()
                .reset_index()
            )
        else:
            grouped = grouped_rep

        if active_group_cols:
            grouped["sample_label"] = interaction_label(grouped, active_group_cols, sep="_")
            sample_labels = _ordered_interaction_levels(grouped, active_group_cols, sep="_")
        else:
            grouped["sample_label"] = "Sample"
            sample_labels = ["Sample"]

        pivot = (
            grouped.pivot_table(
                index="sample_label",
                columns=aa_col,
                values="value",
                fill_value=0,
                sort=False,
            )
            .reindex(index=sample_labels, fill_value=0)
            .reindex(columns=AMINO_ACIDS, fill_value=0)
        )
        pivot.index = range(1, len(pivot) + 1)
        return pivot, sample_labels

    terminals = ["N", "C"] if terminal == "both" else [terminal]
    color_dict = get_color_dict("chemistry", AMINO_ACIDS)
    fig, axes = plt.subplots(1, len(terminals), figsize=(6 * len(terminals), 4), squeeze=False)

    for axis, side in zip(axes.flatten(), terminals):
        logo_frame, sample_labels = build_logo_matrix(side)
        logo = logomaker.Logo(logo_frame, ax=axis, color_scheme="chemistry")
        logo.style_spines(visible=False)
        axis.spines["left"].set_visible(True)
        axis.spines["bottom"].set_visible(True)
        axis.set_xlim(0.5, len(sample_labels) + 0.5)
        axis.set_xticks(range(1, len(sample_labels) + 1))
        axis.set_xticklabels(sample_labels, rotation=45, ha="right")
        axis.set_xlabel("Sample group")
        ylabel = f"{side} terminal {'Intensity' if measure == 'intensity' else 'Peptide type number'}"
        axis.set_ylabel(ylabel)
        ymax = float(axis.get_ylim()[1])
        breaks = breaks_extended(n=4 if measure == "intensity" else 5)((0, ymax))
        breaks = [float(value) for value in breaks if np.isfinite(value)]
        if breaks:
            axis.set_yticks(breaks)
            if breaks[-1] > ymax:
                axis.set_ylim(0, breaks[-1])
        if scientific_10_y:
            axis.yaxis.set_major_formatter(_scientific_formatter())

    chemistry_groups = [
        ("Acidic", color_dict["D"]),
        ("Basic", color_dict["K"]),
        ("Hydrophobic", color_dict["A"]),
        ("Neutral", color_dict["Q"]),
        ("Polar", color_dict["G"]),
    ]
    handles = [
        Patch(facecolor=color, edgecolor=color, label=label)
        for label, color in chemistry_groups
    ]
    if len(terminals) == 1:
        axes[0, 0].legend(
            handles=handles,
            title="chemistry",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
        )
        fig.tight_layout(rect=(0, 0, 0.86, 1))
    else:
        fig.legend(handles=handles, title="chemistry", loc="lower center", ncol=len(handles), frameon=False)
        fig.tight_layout(rect=(0, 0.08, 1, 1))
    return fig
