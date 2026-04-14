from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from peptidomicspy import ttestPeptides


def test_ttest_plain_matches_reference(process_result, reference_frame) -> None:
    actual = ttestPeptides(
        process_result,
        comparisons=[("G120_Y1", "I120_Y1")],
        test_method="plain",
    )["G120_Y1_vs_I120_Y1"]
    reference = reference_frame("ttestPeptides/plain_G120_Y1_vs_I120_Y1.csv.gz")
    assert len(actual) == 2666

    merged = actual.merge(reference[["Sequence", "p.value", "p.adj", "log2FC"]], on="Sequence", suffixes=("_py", "_r"))
    np.testing.assert_allclose(merged["log2FC_py"], merged["log2FC_r"], rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(merged["p.value_py"], merged["p.value_r"], rtol=1e-7, atol=1e-9, equal_nan=True)
    np.testing.assert_allclose(merged["p.adj_py"], merged["p.adj_r"], rtol=1e-7, atol=1e-9, equal_nan=True)


def test_ttest_selector_variants(process_result) -> None:
    direct = ttestPeptides(
        process_result,
        comparisons=[("G120_Y1", "I120_Y1")],
        test_method="plain",
    )["G120_Y1_vs_I120_Y1"]
    reordered = ttestPeptides(
        process_result,
        comparisons=[("Y1_G120", "Y1_I120")],
        test_method="plain",
    )["Y1_G120_vs_Y1_I120"]
    explicit = ttestPeptides(
        process_result,
        comparisons=[("Digest.stage=G120_Yogurt=Y1", "Digest.stage=I120_Yogurt=Y1")],
        test_method="plain",
    )["Digest.stage=G120_Yogurt=Y1_vs_Digest.stage=I120_Yogurt=Y1"]

    direct_map = direct.set_index("Sequence")["log2FC"]
    reordered_map = reordered.set_index("Sequence")["log2FC"]
    explicit_map = explicit.set_index("Sequence")["log2FC"]
    common = direct_map.index.intersection(reordered_map.index).intersection(explicit_map.index)
    np.testing.assert_allclose(direct_map.loc[common], reordered_map.loc[common], rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(direct_map.loc[common], explicit_map.loc[common], rtol=1e-8, atol=1e-8)


def test_ttest_partial_selectors_warn(process_result) -> None:
    with pytest.warns(UserWarning):
        ttestPeptides(
            process_result,
            comparisons=[("Y1", "Y2")],
            test_method="plain",
        )


def test_ttest_treat_approximation_agreement(process_result, reference_frame) -> None:
    actual = ttestPeptides(
        process_result,
        comparisons=[("G120_Y1", "I120_Y1")],
        test_method="treat",
    )["G120_Y1_vs_I120_Y1"]
    reference = reference_frame("ttestPeptides/treat_G120_Y1_vs_I120_Y1.csv.gz")
    assert len(actual) == len(reference) == 2666

    merged = actual.merge(reference[["Sequence", "log2FC", "sig"]], on="Sequence", suffixes=("_py", "_r"))
    np.testing.assert_allclose(merged["log2FC_py"], merged["log2FC_r"], rtol=1e-8, atol=1e-8)
    agreement = (merged["sig_py"] == merged["sig_r"]).mean()
    assert agreement >= 0.95

