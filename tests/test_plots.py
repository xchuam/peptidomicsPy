from __future__ import annotations

from pathlib import Path

from matplotlib.figure import Figure
from plotnine import ggplot

from peptidomicspy import (
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

from conftest import get_plot_label


def test_plot_int_and_type_num_return_ggplot(process_result) -> None:
    assert isinstance(plot_int(process_result), ggplot)
    assert isinstance(plot_type_num(process_result, type="reps"), ggplot)


def test_plot_count_warns_and_returns_ggplot(process_result) -> None:
    import warnings

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        plot = plot_count(process_result, type="reps")
    assert isinstance(plot, ggplot)
    assert any("deprecated" in str(item.message) for item in captured)


def test_density_plot_labels_match_reference(process_result, reference_manifest) -> None:
    expected = {item["case"]: item["colour_label"] for item in reference_manifest["plot_reference"]}

    length_plot = plot_length_distribution(
        process_result,
        type="reps",
        metric="count",
        plot_mode="density",
        facet_rows="Digest.stage",
    )
    gravy_plot = plot_gravy_vs_intensity(
        process_result,
        type="reps",
        plot_mode="density",
        facet_rows="Digest.stage",
        filter_params={"Yogurt": "Y1"},
    )
    assert get_plot_label(length_plot, "color") == expected["plot_length_distribution_density"]
    assert get_plot_label(gravy_plot, "color") == expected["plot_gravy_vs_intensity_density"]


def test_plot_volcano_returns_named_ggplots(process_result) -> None:
    stats_result = ttestPeptides(
        process_result,
        comparisons=[("G120_Y1", "I120_Y1")],
        test_method="plain",
    )
    plots = plot_volcano(stats_result, comparisons=[("G120_Y1", "I120_Y1")], test_method="plain")
    assert list(plots) == ["G120_Y1_vs_I120_Y1"]
    assert isinstance(plots["G120_Y1_vs_I120_Y1"], ggplot)


def test_plot_pep_align_saves_and_sets_autosize(process_result, tmp_path: Path) -> None:
    output_path = tmp_path / "pep_align.png"
    figure = plot_pep_align(
        process_result,
        protein_name="P81265",
        protein_seq="A" * 650,
        save_file_location=str(output_path),
    )
    assert isinstance(figure, Figure)
    assert output_path.exists()
    assert hasattr(figure, "suggested_fig_width")
    assert hasattr(figure, "suggested_fig_height")
    assert hasattr(figure, "suggested_fig_size")


def test_plot_cleavage_site_supports_modes(process_result) -> None:
    assert isinstance(plot_cleavage_site(process_result, terminal="N"), Figure)
    assert isinstance(plot_cleavage_site(process_result, terminal="C"), Figure)
    assert isinstance(plot_cleavage_site(process_result, terminal="both"), Figure)
