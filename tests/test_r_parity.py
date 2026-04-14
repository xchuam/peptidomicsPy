from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure
from plotnine import ggplot

from peptidomicspy import (
    filterPeptides,
    plot_cleavage_site,
    plot_count,
    plot_gravy_vs_intensity,
    plot_int,
    plot_length_distribution,
    plot_pep_align,
    plot_type_num,
    plot_volcano,
    ttestPeptides,
)

from conftest import assert_dataframe_close, image_mae


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def _compare_ttest_core(actual: pd.DataFrame, reference: pd.DataFrame) -> None:
    assert len(actual) == len(reference)
    merged = actual.merge(
        reference[["Sequence", "t", "df", "p.value", "log2FC", "p.adj", "sig"]],
        on="Sequence",
        suffixes=("_py", "_r"),
    )
    np.testing.assert_allclose(merged["t_py"], merged["t_r"], rtol=1e-7, atol=1e-9, equal_nan=True)
    np.testing.assert_allclose(merged["df_py"], merged["df_r"], rtol=1e-7, atol=1e-9, equal_nan=True)
    np.testing.assert_allclose(merged["p.value_py"], merged["p.value_r"], rtol=1e-7, atol=1e-9, equal_nan=True)
    np.testing.assert_allclose(merged["log2FC_py"], merged["log2FC_r"], rtol=1e-8, atol=1e-8, equal_nan=True)
    np.testing.assert_allclose(merged["p.adj_py"], merged["p.adj_r"], rtol=1e-7, atol=1e-9, equal_nan=True)
    assert (merged["sig_py"] == merged["sig_r"]).all()


def _save_plot_object(obj, path: Path, width: float, height: float) -> None:
    if isinstance(obj, ggplot):
        obj.save(path, width=width, height=height, dpi=150, verbose=False)
        return
    if isinstance(obj, Figure):
        obj.set_size_inches(width, height)
        obj.savefig(path, dpi=150)
        plt.close(obj)
        return
    raise TypeError(f"Unsupported plot object type: {type(obj)!r}")


def test_filter_r_parity_cases(process_result, r_parity_dir: Path) -> None:
    cases = {
        "exact": {"seqs": ["AAGGPGAPADPGRPT", "AAIDEASKKLNAQ"]},
        "regex": {"seq_pattern": "^AA"},
        "grouped": {"seq_pattern": "DQAM", "filter_params": {"Yogurt": "Y1", "Digest.stage": "G120"}},
        "union": {"seqs": ["AAGGPGAPADPGRPT"], "seq_pattern": "^AAI"},
    }
    table_names = [
        "dt_peptides",
        "dt_peptides_int",
        "dt_peptides_int_reps",
        "dt_peptides_typenum",
        "dt_peptides_typenum_reps",
        "dt_int_col",
    ]

    for case_name, kwargs in cases.items():
        actual = filterPeptides(process_result, **kwargs)
        case_dir = r_parity_dir / "tables" / "filter" / case_name
        for table_name in table_names:
            assert_dataframe_close(getattr(actual, table_name), _read_csv(case_dir / f"{table_name}.csv.gz"))


def test_filter_updates_ttest_table(process_result) -> None:
    filtered = filterPeptides(
        process_result,
        seq_pattern="DQAM",
        filter_params={"Yogurt": "Y1", "Digest.stage": "G120"},
    )
    assert set(filtered.dt_peptides_int_ttest["Yogurt"].astype(str).unique()) == {"Y1"}
    assert set(filtered.dt_peptides_int_ttest["Digest.stage"].astype(str).unique()) == {"G120"}
    assert set(filtered.dt_peptides_int_ttest["Sequence"].astype(str).unique()).issubset(
        set(filtered.dt_peptides["Sequence"].astype(str).unique())
    )


def test_ttest_plain_r_parity_cases(process_result, r_parity_dir: Path) -> None:
    cases = {
        "plain_main": {
            "comparisons": [("G120_Y1", "I120_Y1")],
            "result_key": "G120_Y1_vs_I120_Y1",
        },
        "plain_reordered": {
            "comparisons": [("Y1_G120", "Y1_I120")],
            "result_key": "Y1_G120_vs_Y1_I120",
        },
        "plain_explicit": {
            "comparisons": [("Digest.stage=G120_Yogurt=Y1", "Digest.stage=I120_Yogurt=Y1")],
            "result_key": "Digest.stage=G120_Yogurt=Y1_vs_Digest.stage=I120_Yogurt=Y1",
        },
        "plain_pooled": {
            "comparisons": [("Y1", "Y2")],
            "result_key": "Y1_vs_Y2",
        },
        "plain_equalvar": {
            "comparisons": [("G120_Y1", "I120_Y1")],
            "result_key": "G120_Y1_vs_I120_Y1",
            "equal_var": True,
        },
    }

    for case_name, kwargs in cases.items():
        result_key = kwargs["result_key"]
        test_kwargs = {key: value for key, value in kwargs.items() if key != "result_key"}
        actual = ttestPeptides(process_result, test_method="plain", **test_kwargs)[result_key]
        reference = _read_csv(r_parity_dir / "tables" / "ttest" / f"{case_name}.csv.gz")
        _compare_ttest_core(actual, reference)


def test_ttest_treat_r_parity_main(process_result, r_parity_dir: Path) -> None:
    actual = ttestPeptides(
        process_result,
        comparisons=[("G120_Y1", "I120_Y1")],
        test_method="treat",
    )["G120_Y1_vs_I120_Y1"]
    reference = _read_csv(r_parity_dir / "tables" / "ttest" / "treat_main.csv.gz")
    assert len(actual) == len(reference)
    merged = actual.merge(reference[["Sequence", "log2FC", "sig"]], on="Sequence", suffixes=("_py", "_r"))
    np.testing.assert_allclose(merged["log2FC_py"], merged["log2FC_r"], rtol=1e-8, atol=1e-8, equal_nan=True)
    assert (merged["sig_py"] == merged["sig_r"]).mean() >= 0.95


@pytest.mark.parametrize(
    ("case_name", "threshold"),
    [
        ("plot_int_mean_named", 0.12),
        ("plot_int_reps_default", 0.07),
        ("plot_type_num_reps", 0.08),
        ("plot_count_mean_none", 0.03),
        ("plot_length_bar", 0.08),
        ("plot_length_density", 0.05),
        ("plot_gravy_scatter", 0.08),
        ("plot_gravy_density", 0.05),
        ("plot_volcano_plain", 0.05),
        ("plot_volcano_treat_annotated", 0.06),
        ("plot_pep_align_mean", 0.10),
        ("plot_cleavage_N", 0.20),
        ("plot_cleavage_both_filtered", 0.23),
    ],
)
def test_plot_r_parity(case_name: str, threshold: float, process_result, r_parity_dir: Path, tmp_path: Path) -> None:
    plain = ttestPeptides(
        process_result,
        comparisons=[("G120_Y1", "I120_Y1")],
        test_method="plain",
    )
    treat_reference = _read_csv(r_parity_dir / "tables" / "ttest" / "treat_main.csv.gz")
    treat_labels = list(treat_reference["Sequence"].head(2))
    treat_key = "G120_Y1_vs_I120_Y1"
    treat_reference_dict = {treat_key: treat_reference}

    builders = {
        "plot_int_mean_named": lambda: (plot_int(process_result, type="mean", x_var="Yogurt", filter_params={"Digest.stage": "G120"}, color_by="Protein.name"), 8, 5),
        "plot_int_reps_default": lambda: (plot_int(process_result, type="reps"), 8, 5),
        "plot_type_num_reps": lambda: (plot_type_num(process_result, type="reps", x_var="Yogurt", facet_rows="Digest.stage", color_by="Protein.group"), 8, 5),
        "plot_count_mean_none": lambda: (plot_count(process_result, type="mean", x_var="Yogurt", color_by="none"), 8, 5),
        "plot_length_bar": lambda: (plot_length_distribution(process_result, metric="intensity", filter_params={"Digest.stage": "G120"}, facet_rows="Yogurt"), 8, 5),
        "plot_length_density": lambda: (plot_length_distribution(process_result, type="reps", metric="type_num", plot_mode="density", facet_rows="Digest.stage"), 8, 5),
        "plot_gravy_scatter": lambda: (plot_gravy_vs_intensity(process_result), 8, 5),
        "plot_gravy_density": lambda: (plot_gravy_vs_intensity(process_result, type="reps", plot_mode="density", facet_rows="Digest.stage", filter_params={"Yogurt": "Y1"}), 8, 5),
        "plot_volcano_plain": lambda: (plot_volcano(plain, comparisons=[("G120_Y1", "I120_Y1")], test_method="plain")["G120_Y1_vs_I120_Y1"], 8, 5),
        "plot_volcano_treat_annotated": lambda: (plot_volcano(treat_reference_dict, comparisons=[("G120_Y1", "I120_Y1")], test_method="treat", label_seqs=treat_labels, highlight_seqs={"black": [treat_labels[0]]})[treat_key], 8, 5),
        "plot_pep_align_mean": lambda: (plot_pep_align(process_result, protein_name="P81265", protein_seq="A" * 650, auto_size=False), 12, 8),
        "plot_cleavage_N": lambda: (plot_cleavage_site(process_result, terminal="N"), 8, 5),
        "plot_cleavage_both_filtered": lambda: (plot_cleavage_site(process_result, terminal="both", measure="type_num", replicate_mode="reps", filter_params={"Yogurt": "Y1"}), 12, 5),
    }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj, width, height = builders[case_name]()
    py_path = tmp_path / f"{case_name}.png"
    _save_plot_object(obj, py_path, width, height)
    r_path = r_parity_dir / "plots" / f"{case_name}.png"
    assert image_mae(py_path, r_path) <= threshold
