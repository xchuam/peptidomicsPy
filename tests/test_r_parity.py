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
    calculate_gravy,
    filterPeptides,
    load_example_protein_mapping,
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


def _compare_treat_core(actual: pd.DataFrame, reference: pd.DataFrame, min_sig_agreement: float = 0.95) -> None:
    assert len(actual) == len(reference)
    merged = actual.merge(reference[["Sequence", "log2FC", "sig"]], on="Sequence", suffixes=("_py", "_r"))
    np.testing.assert_allclose(merged["log2FC_py"], merged["log2FC_r"], rtol=1e-8, atol=1e-8, equal_nan=True)
    assert (merged["sig_py"] == merged["sig_r"]).mean() >= min_sig_agreement


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


def _pep_align_custom_result(process_result):
    custom = process_result.copy()
    custom.dt_peptides_int_reps = custom.dt_peptides_int_reps.rename(
        columns={
            "Sequence": "Seq.custom",
            "Start.position": "Start.custom",
            "End.position": "End.custom",
            "Intensity": "Intensity.custom",
        }
    )
    return custom


def _pep_align_metadata_row(r_parity_dir: Path, case_name: str) -> pd.Series:
    metadata = pd.read_csv(r_parity_dir / "tables" / "helpers" / "plot_pep_align_metadata.csv")
    row = metadata.loc[metadata["case"] == case_name]
    assert not row.empty
    return row.iloc[0]


def test_mapping_and_gravy_helpers_r_parity(r_parity_dir: Path) -> None:
    mapping_reference = _read_csv(r_parity_dir / "tables" / "helpers" / "mapping_example_dataset.csv.gz")
    assert_dataframe_close(load_example_protein_mapping(), mapping_reference)

    gravy_reference = _read_csv(r_parity_dir / "tables" / "helpers" / "calculate_gravy_examples.csv.gz")
    actual = pd.DataFrame(
        {
            "peptide": gravy_reference["peptide"],
            "gravy": gravy_reference["peptide"].map(calculate_gravy),
        }
    )
    assert_dataframe_close(actual, gravy_reference)


def test_process_default_r_parity(process_result, r_parity_dir: Path) -> None:
    process_dir = r_parity_dir / "tables" / "process" / "default"
    table_names = [
        "dt_peptides",
        "dt_peptides_int",
        "dt_peptides_int_reps",
        "dt_peptides_typenum",
        "dt_peptides_typenum_reps",
        "dt_int_col",
    ]
    for table_name in table_names:
        assert_dataframe_close(getattr(process_result, table_name), _read_csv(process_dir / f"{table_name}.csv.gz"))

    metadata = pd.read_csv(process_dir / "metadata.csv").iloc[0]
    assert process_result.grp_cols == metadata["grp_cols"].split("|")
    assert process_result.peptides_select_col_basic == metadata["peptides_select_col_basic"].split("|")


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
        "plain_custom_adjust": {
            "comparisons": [("G120_Y1", "I120_Y1")],
            "result_key": "G120_Y1_vs_I120_Y1",
            "pseudocount": 5,
            "adjust": "bonferroni",
        },
        "plain_min_reps": {
            "comparisons": [("G120_Y1", "I120_Y1")],
            "result_key": "G120_Y1_vs_I120_Y1",
            "min_reps_per_side": 3,
        },
    }

    for case_name, kwargs in cases.items():
        result_key = kwargs["result_key"]
        test_kwargs = {key: value for key, value in kwargs.items() if key != "result_key"}
        actual = ttestPeptides(process_result, test_method="plain", **test_kwargs)[result_key]
        reference = _read_csv(r_parity_dir / "tables" / "ttest" / f"{case_name}.csv.gz")
        _compare_ttest_core(actual, reference)


@pytest.mark.parametrize(
    ("case_name", "kwargs"),
    [
        ("treat_main", {"comparisons": [("G120_Y1", "I120_Y1")], "result_key": "G120_Y1_vs_I120_Y1"}),
        (
            "treat_custom_threshold",
            {
                "comparisons": [("G120_Y1", "I120_Y1")],
                "result_key": "G120_Y1_vs_I120_Y1",
                "lfc_thresh": 0.5,
                "alpha": 0.1,
            },
        ),
    ],
)
def test_ttest_treat_r_parity_cases(case_name: str, kwargs: dict[str, object], process_result, r_parity_dir: Path) -> None:
    result_key = kwargs["result_key"]
    test_kwargs = {key: value for key, value in kwargs.items() if key != "result_key"}
    actual = ttestPeptides(process_result, test_method="treat", **test_kwargs)[result_key]
    reference = _read_csv(r_parity_dir / "tables" / "ttest" / f"{case_name}.csv.gz")
    _compare_treat_core(actual, reference)


def test_plot_pep_align_metadata_r_parity(process_result, r_parity_dir: Path) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lookup_plot = plot_pep_align(process_result, protein_name="P02662")
        custom_plot = plot_pep_align(
            _pep_align_custom_result(process_result),
            protein_name="P02662",
            type="reps",
            filter_params={"Yogurt": "Y1"},
            seq_col="Seq.custom",
            start_col="Start.custom",
            end_col="End.custom",
            intensity_col="Intensity.custom",
        )

    for case_name, figure in {
        "plot_pep_align_lookup_default": lookup_plot,
        "plot_pep_align_reps_custom_cols": custom_plot,
    }.items():
        expected = _pep_align_metadata_row(r_parity_dir, case_name)
        assert figure.suggested_fig_width == pytest.approx(expected["suggested_fig_width"])
        assert figure.suggested_fig_height == pytest.approx(expected["suggested_fig_height"])
        plt.close(figure)


@pytest.mark.parametrize(
    ("case_name", "threshold"),
    [
        ("plot_int_mean_named", 0.12),
        ("plot_int_reps_default", 0.07),
        ("plot_int_faceted_linear", 0.08),
        ("plot_type_num_reps", 0.08),
        ("plot_type_num_filtered_scientific", 0.08),
        ("plot_count_mean_none", 0.03),
        ("plot_count_reps_faceted", 0.08),
        ("plot_length_bar", 0.08),
        ("plot_length_density", 0.05),
        ("plot_length_bar_named_linear", 0.11),
        ("plot_gravy_scatter", 0.08),
        ("plot_gravy_density", 0.05),
        ("plot_gravy_scatter_named", 0.08),
        ("plot_volcano_plain", 0.05),
        ("plot_volcano_treat_annotated", 0.06),
        ("plot_volcano_custom_ylog", 0.07),
        ("plot_pep_align_lookup_default", 0.24),
        ("plot_pep_align_mean", 0.10),
        ("plot_pep_align_ranged_labels", 0.18),
        ("plot_pep_align_reps_custom_cols", 0.22),
        ("plot_cleavage_N", 0.20),
        ("plot_cleavage_both_filtered", 0.23),
        ("plot_cleavage_C_keep_constant_linear", 0.24),
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
    plain_labels = list(plain[treat_key]["Sequence"].head(2))

    builders = {
        "plot_int_mean_named": lambda: (
            plot_int(process_result, type="mean", x_var="Yogurt", filter_params={"Digest.stage": "G120"}, color_by="Protein.name"),
            8,
            5,
        ),
        "plot_int_reps_default": lambda: (plot_int(process_result, type="reps"), 8, 5),
        "plot_int_faceted_linear": lambda: (
            plot_int(
                process_result,
                type="reps",
                x_var="Digest.stage",
                color_by="Protein.group",
                facet_rows="Yogurt",
                facet_cols="Replicate",
                scientific_10_y=False,
            ),
            10,
            6,
        ),
        "plot_type_num_reps": lambda: (
            plot_type_num(process_result, type="reps", x_var="Yogurt", facet_rows="Digest.stage", color_by="Protein.group"),
            8,
            5,
        ),
        "plot_type_num_filtered_scientific": lambda: (
            plot_type_num(
                process_result,
                type="reps",
                x_var="Digest.stage",
                filter_params={"Yogurt": "Y1"},
                facet_cols="Replicate",
                scientific_10_y=True,
            ),
            9,
            5,
        ),
        "plot_count_mean_none": lambda: (plot_count(process_result, type="mean", x_var="Yogurt", color_by="none"), 8, 5),
        "plot_count_reps_faceted": lambda: (
            plot_count(
                process_result,
                type="reps",
                x_var="Yogurt",
                color_by="Protein.group",
                filter_params={"Digest.stage": "G120"},
                facet_rows="Digest.stage",
                facet_cols="Replicate",
                scientific_10_y=True,
            ),
            10,
            6,
        ),
        "plot_length_bar": lambda: (
            plot_length_distribution(process_result, metric="intensity", filter_params={"Digest.stage": "G120"}, facet_rows="Yogurt"),
            8,
            5,
        ),
        "plot_length_density": lambda: (
            plot_length_distribution(process_result, type="reps", metric="type_num", plot_mode="density", facet_rows="Digest.stage"),
            8,
            5,
        ),
        "plot_length_bar_named_linear": lambda: (
            plot_length_distribution(
                process_result,
                type="reps",
                metric="intensity",
                color_by="Protein.name",
                facet_rows="Digest.stage",
                facet_cols="Yogurt+Replicate",
                scientific_10_y=False,
            ),
            10,
            6,
        ),
        "plot_gravy_scatter": lambda: (plot_gravy_vs_intensity(process_result), 8, 5),
        "plot_gravy_density": lambda: (
            plot_gravy_vs_intensity(process_result, type="reps", plot_mode="density", facet_rows="Digest.stage", filter_params={"Yogurt": "Y1"}),
            8,
            5,
        ),
        "plot_gravy_scatter_named": lambda: (
            plot_gravy_vs_intensity(process_result, color_by="Protein.name", facet_cols="Digest.stage", alpha_value=0.4),
            10,
            5,
        ),
        "plot_volcano_plain": lambda: (
            plot_volcano(plain, comparisons=[("G120_Y1", "I120_Y1")], test_method="plain")["G120_Y1_vs_I120_Y1"],
            8,
            5,
        ),
        "plot_volcano_treat_annotated": lambda: (
            plot_volcano(
                treat_reference_dict,
                comparisons=[("G120_Y1", "I120_Y1")],
                test_method="treat",
                label_seqs=treat_labels,
                highlight_seqs={"black": [treat_labels[0]]},
            )[treat_key],
            8,
            5,
        ),
        "plot_volcano_custom_ylog": lambda: (
            plot_volcano(
                plain,
                comparisons=1,
                test_method="plain",
                show_threshold=True,
                lfc_thresh=0.5,
                alpha=0.1,
                fill_values={"no": "#AAAAAA", "yes": "#E41A1C"},
                point_size=2.5,
                point_alpha=0.6,
                label_seqs=plain_labels,
                label_size=4,
                label_col="navy",
                highlight_seqs={"black": [plain_labels[0]], "#377EB8": [plain_labels[1]]},
                highlight_size=3,
                highlight_stroke=2,
                y_log_scale=True,
            )[treat_key],
            8,
            5,
        ),
        "plot_pep_align_lookup_default": lambda: (plot_pep_align(process_result, protein_name="P02662"), 8, 5),
        "plot_pep_align_mean": lambda: (
            plot_pep_align(process_result, protein_name="P81265", protein_seq="A" * 650, auto_size=False),
            12,
            8,
        ),
        "plot_pep_align_ranged_labels": lambda: (
            plot_pep_align(
                process_result,
                protein_name="P02662",
                filter_params={"Digest.stage": "G120"},
                x_interval=5,
                x_range=(80, 165),
                y_range=(0, 12),
                label_seq={
                    "highlight_1": ["DQAMEDIKQ", "DQAMEDIKQM"],
                    "highlight_2": ["AMEDIKQM"],
                },
                label_col={"highlight_1": "red", "highlight_2": "blue"},
            ),
            9,
            6,
        ),
        "plot_pep_align_reps_custom_cols": lambda: (
            plot_pep_align(
                _pep_align_custom_result(process_result),
                protein_name="P02662",
                type="reps",
                filter_params={"Yogurt": "Y1"},
                seq_col="Seq.custom",
                start_col="Start.custom",
                end_col="End.custom",
                intensity_col="Intensity.custom",
            ),
            9,
            6,
        ),
        "plot_cleavage_N": lambda: (plot_cleavage_site(process_result, terminal="N"), 8, 5),
        "plot_cleavage_both_filtered": lambda: (
            plot_cleavage_site(process_result, terminal="both", measure="type_num", replicate_mode="reps", filter_params={"Yogurt": "Y1"}),
            12,
            5,
        ),
        "plot_cleavage_C_keep_constant_linear": lambda: (
            plot_cleavage_site(
                process_result,
                terminal="C",
                filter_params={"Yogurt": "Y1", "Digest.stage": "G120"},
                scientific_10_y=False,
                drop_constant_groups=False,
            ),
            8,
            5,
        ),
    }

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj, width, height = builders[case_name]()
    py_path = tmp_path / f"{case_name}.png"
    _save_plot_object(obj, py_path, width, height)
    r_path = r_parity_dir / "plots" / f"{case_name}.png"
    assert image_mae(py_path, r_path) <= threshold


def test_plot_pep_align_saved_r_parity(process_result, r_parity_dir: Path, tmp_path: Path) -> None:
    save_path = tmp_path / "plot_pep_align_saved.png"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        figure = plot_pep_align(
            process_result,
            protein_name="P02662",
            filter_params={"Digest.stage": "G120"},
            auto_size=True,
            save_file_location=str(save_path),
            save_file_dpi=200,
            width_per_aa=0.12,
            height_per_row=0.2,
            min_width=5,
            min_height=4,
        )
    expected = _pep_align_metadata_row(r_parity_dir, "plot_pep_align_saved")
    assert figure.suggested_fig_width == pytest.approx(expected["suggested_fig_width"])
    assert figure.suggested_fig_height == pytest.approx(expected["suggested_fig_height"])
    assert getattr(figure, "save_file_location") == str(save_path)
    assert getattr(figure, "save_file_dpi") == 200
    assert save_path.exists()
    assert image_mae(save_path, r_parity_dir / "plots" / "plot_pep_align_saved.png") <= 0.12
    plt.close(figure)
