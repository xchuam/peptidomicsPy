from __future__ import annotations

from peptidomicspy import filterPeptides

from conftest import assert_dataframe_close


def test_filter_exact_matches_reference(process_result, reference_frame) -> None:
    actual = filterPeptides(
        process_result,
        seqs=["AAGGPGAPADPGRPT", "AAIDEASKKLNAQ"],
    )
    assert_dataframe_close(actual.dt_peptides, reference_frame("filterPeptides/exact/dt_peptides.csv.gz"))
    assert_dataframe_close(actual.dt_peptides_int, reference_frame("filterPeptides/exact/dt_peptides_int.csv.gz"))
    assert_dataframe_close(
        actual.dt_peptides_int_reps,
        reference_frame("filterPeptides/exact/dt_peptides_int_reps.csv.gz"),
    )


def test_filter_regex_matches_reference(process_result, reference_frame) -> None:
    actual = filterPeptides(process_result, seq_pattern="^AA")
    reference = reference_frame("filterPeptides/regex/dt_peptides.csv.gz")
    assert sorted(actual.dt_peptides["Sequence"].unique().tolist()) == sorted(
        reference["Sequence"].unique().tolist()
    )


def test_filter_grouping_filters_and_columns(process_result, reference_frame) -> None:
    actual = filterPeptides(
        process_result,
        seq_pattern="DQAM",
        filter_params={"Yogurt": "Y1", "Digest.stage": "G120"},
    )
    assert actual.dt_peptides_int["Yogurt"].astype(str).unique().tolist() == ["Y1"]
    assert actual.dt_peptides_int["Digest.stage"].astype(str).unique().tolist() == ["G120"]
    assert_dataframe_close(
        actual.dt_int_col,
        reference_frame("filterPeptides/grouped/dt_int_col.csv.gz"),
    )
    assert_dataframe_close(
        actual.dt_peptides_int,
        reference_frame("filterPeptides/grouped/dt_peptides_int.csv.gz"),
    )
