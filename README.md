# peptidomicspy

Python migration of `peptidomicsR` for peptidomics and protein digestion analysis.

This branch targets the upstream R tag `v1.1.1-alpha` and keeps the first Python release aligned to that migration target as package version `1.1.1a0`.

## What is included

- `processPeptides()` for MaxQuant peptide processing
- `filterPeptides()` for sequence and grouping filters
- `ttestPeptides()` for plain and thresholded comparisons
- `plot_int()`, `plot_type_num()`, `plot_length_distribution()`, `plot_gravy_vs_intensity()`, `plot_volcano()`
- `plot_pep_align()` and `plot_cleavage_site()` with Matplotlib-based implementations
- `load_example_protein_mapping()` for the shipped example mapping dataset

The public function names intentionally follow the R package for this migration target.

## Install

```bash
python -m venv .venv
.venv/bin/pip install -e '.[dev]'
```

## Quick start

```python
from peptidomicspy import processPeptides, plot_int, ttestPeptides

result = processPeptides(
    peptides_file="tests/fixtures/v1.1.1-alpha/Yogurtexample_QR188-205.csv",
    intensity_columns_file="tests/fixtures/v1.1.1-alpha/Intensity_columns.csv",
    protein_mapping_file="tests/fixtures/v1.1.1-alpha/protein_mapping.csv",
)

plot = plot_int(result, type="mean", x_var="Yogurt")

stats = ttestPeptides(
    result,
    comparisons=[("G120_Y1", "I120_Y1")],
    test_method="plain",
)
```

## Data layout

- Test fixtures live under `tests/fixtures/v1.1.1-alpha/`
- Frozen R parity outputs live under `.migration-artifacts/v1.1.1-alpha/`
- Package data lives under `src/peptidomicspy/data/`

## Status

This is a migration branch implementation, not the final user-facing `main` branch release. Package-facing code, tests, and documentation are present here so the `v1.1.1-alpha` port can be validated before promotion.
