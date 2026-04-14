from __future__ import annotations

import pandas as pd

from conftest import assert_dataframe_close


def test_process_row_counts(process_result) -> None:
    assert process_result.grp_cols == ["Digest.stage", "Yogurt"]
    assert len(process_result.dt_peptides) == 3155
    assert len(process_result.dt_peptides_int) == 7671
    assert len(process_result.dt_peptides_int_reps) == 21825
    assert len(process_result.dt_peptides_int_ttest) == 18930
    assert len(process_result.dt_peptides_typenum) == 512
    assert len(process_result.dt_peptides_typenum_reps) == 1512


def test_process_group_and_category_order(process_result, fixture_dir) -> None:
    assert list(process_result.dt_peptides_int["Digest.stage"].cat.categories) == ["G120", "I120"]
    assert list(process_result.dt_peptides_int["Yogurt"].cat.categories) == ["Y1", "Y2", "Y3"]

    mapping = pd.read_csv(fixture_dir / "protein_mapping.csv", sep=";", encoding="utf-8-sig")
    expected_names = [*mapping["Protein.name"].drop_duplicates().tolist(), "Others"]
    expected_groups = mapping["Protein.group"].drop_duplicates().tolist()
    assert list(process_result.dt_peptides["Protein.name"].cat.categories) == expected_names
    assert list(process_result.dt_peptides["Protein.group"].cat.categories) == expected_groups


def test_process_reference_tables_match(process_result, reference_frame) -> None:
    table_map = {
        "dt_int_col": "processPeptides/dt_int_col.csv.gz",
        "dt_peptides": "processPeptides/dt_peptides.csv.gz",
        "dt_peptides_int": "processPeptides/dt_peptides_int.csv.gz",
        "dt_peptides_int_reps": "processPeptides/dt_peptides_int_reps.csv.gz",
        "dt_peptides_int_ttest": "processPeptides/dt_peptides_int_ttest.csv.gz",
        "dt_peptides_typenum": "processPeptides/dt_peptides_typenum.csv.gz",
        "dt_peptides_typenum_reps": "processPeptides/dt_peptides_typenum_reps.csv.gz",
    }
    for attr, relpath in table_map.items():
        actual = getattr(process_result, attr)
        reference = reference_frame(relpath)
        assert_dataframe_close(actual, reference)
