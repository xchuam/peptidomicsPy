from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path


VERSION_LABEL = "v1.1.1-alpha"


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def gzip_csv_summary(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        rows = list(reader)
    header = rows[0] if rows else []
    data_rows = rows[1:] if len(rows) > 1 else []
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "path": str(path.relative_to(path.parents[3])),
        "rows": len(data_rows),
        "columns": header,
        "sha256": digest,
    }


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    upstream_dir = Path("/tmp/peptidomicsR-v1.1.1-alpha")
    fixtures_dir = repo_root / "tests" / "fixtures" / VERSION_LABEL
    artifacts_dir = repo_root / ".migration-artifacts" / VERSION_LABEL
    reference_dir = artifacts_dir / "reference_outputs"
    reference_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "Rscript",
            str(repo_root / "tools" / "generate_r_reference_outputs.R"),
            str(upstream_dir),
            str(fixtures_dir),
            str(reference_dir),
        ],
        check=True,
    )

    description = read_text(upstream_dir / "DESCRIPTION")
    namespace = read_text(upstream_dir / "NAMESPACE")
    readme = read_text(upstream_dir / "README.md")
    commit_sha = (
        subprocess.check_output(["git", "-C", str(upstream_dir), "rev-parse", "HEAD"], text=True).strip()
    )

    exports = []
    for line in namespace.splitlines():
        if line.startswith("export("):
            exports.append(line.removeprefix("export(").removesuffix(")"))

    (artifacts_dir / "EXPORTED_API.csv").write_text(
        "function_name\n" + "\n".join(exports) + "\n",
        encoding="utf-8",
    )

    dependency_lines = []
    capture = None
    for line in description.splitlines():
        if line.startswith("License:"):
            license_line = line.split(":", 1)[1].strip()
        if line.startswith("Imports:"):
            capture = "Imports"
            dependency_lines.append("## Imports")
            continue
        if line.startswith("Suggests:"):
            capture = "Suggests"
            dependency_lines.append("\n## Suggests")
            continue
        if line and not line.startswith(" ") and capture:
            capture = None
        if capture and line.strip():
            dependency_lines.append(f"- {line.strip().rstrip(',')}")

    dependency_map = "\n".join(
        [
            "# Dependency Map",
            "",
            f"Upstream license: `{license_line}`",
            "",
            *dependency_lines,
            "",
            "Key Python runtime counterparts:",
            "- `data.table` -> `pandas`",
            "- `ggplot2` and wrappers -> `plotnine`",
            "- `limma::treat` -> documented pure-Python approximation",
            "- `ggseqlogo` -> `logomaker`",
        ]
    )
    (artifacts_dir / "DEPENDENCY_MAP.md").write_text(dependency_map, encoding="utf-8")

    dossier = "\n".join(
        [
            "# R Package Dossier",
            "",
            f"- Upstream repository: `xchuam/peptidomicsR`",
            f"- Frozen tag: `{VERSION_LABEL}`",
            f"- Frozen commit: `{commit_sha}`",
            "- Version in DESCRIPTION: `1.1.1.9000`",
            "- Repo canonical migration label: `v1.1.1-alpha`",
            "- Description-to-branch mapping: `1.1.1.9000` -> `v1.1.1-alpha`",
            f"- Exported functions: {', '.join(exports)}",
            "- Missing paths verified upstream: `vignettes/`, `src/`, `tests/`",
            "- Present paths verified upstream: `R/`, `man/`, `test/`, `data/`, `README.md`, `README.Rmd`",
            "",
            "## README Summary",
            "",
            readme.split("## 📦 Installation", 1)[0].strip(),
        ]
    )
    (artifacts_dir / "R_PACKAGE_DOSSIER.md").write_text(dossier + "\n", encoding="utf-8")

    data_inventory = "\n".join(
        [
            "# Data Asset Inventory",
            "",
            "- Upstream packaged R data:",
            "  - `data/peptidomics_protein_mapping_example.rda`",
            "- Local parity fixtures:",
            "  - `tests/fixtures/v1.1.1-alpha/Yogurtexample_QR188-205.csv`",
            "  - `tests/fixtures/v1.1.1-alpha/Intensity_columns.csv`",
            "  - `tests/fixtures/v1.1.1-alpha/protein_mapping.csv`",
            "- Runtime package data:",
            "  - `src/peptidomicspy/data/peptidomics_protein_mapping_example.csv`",
        ]
    )
    (artifacts_dir / "DATA_ASSET_INVENTORY.md").write_text(data_inventory + "\n", encoding="utf-8")

    (artifacts_dir / "COMPILED_CODE_INVENTORY.md").write_text(
        "# Compiled Code Inventory\n\nNo upstream `src/` directory or compiled code is present in `v1.1.1-alpha`.\n",
        encoding="utf-8",
    )

    (artifacts_dir / "MIGRATION_MODE_MAP.md").write_text(
        "\n".join(
            [
                "# Migration Mode Map",
                "",
                "| Subsystem | Mode | Notes |",
                "| --- | --- | --- |",
                "| `processPeptides` | native rewrite | Preserve result schema and grouping semantics with pandas. |",
                "| `filterPeptides` | faithful port | Preserve selector and regrouping behavior, including current sequence-first filtering semantics. |",
                "| `ttestPeptides(plain)` | faithful port | Match selector parsing and Welch/Student behavior on log2 intensities. |",
                "| `ttestPeptides(treat)` | native rewrite with documented deviation | Pure-Python thresholded t approximation; not limma parity. |",
                "| Bar and density plots | native rewrite | Preserve parameter surface while switching to plotnine. |",
                "| `plot_pep_align` | native rewrite | Matplotlib implementation for stacked peptide alignment layout. |",
                "| `plot_cleavage_site` | native rewrite | Matplotlib plus logomaker implementation. |",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (artifacts_dir / "PARITY_PRIORITY_MATRIX.md").write_text(
        "\n".join(
            [
                "# Parity Priority Matrix",
                "",
                "| Priority | Surface | Reason |",
                "| --- | --- | --- |",
                "| P0 | `processPeptides` result tables and group ordering | All downstream API depends on this shape. |",
                "| P0 | `filterPeptides` sequence and grouping semantics | Drives downstream subset workflows. |",
                "| P0 | `ttestPeptides(plain)` | Numeric parity check for core statistics. |",
                "| P1 | `ttestPeptides(treat)` significance agreement | Documented approximation target. |",
                "| P1 | Plot return types and default labels | Required for user workflow continuity. |",
                "| P2 | Plot visual styling | Backend changes are acceptable if semantics remain aligned. |",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (artifacts_dir / "SEMANTIC_RULES.md").write_text(
        "\n".join(
            [
                "# Semantic Rules",
                "",
                "- Replace spaces with dots in all peptide column names before downstream processing.",
                "- Parse the local fixture CSVs as semicolon-delimited UTF-8 with BOM handling.",
                "- Remove peptides whose `Leading.razor.protein` starts with `CON_` or `REV_`.",
                "- Use `Intensity.column` order from the metadata file to define `R1..Rn` replicate columns.",
                "- Keep peptides in mean tables only if at least half of replicates are strictly greater than zero.",
                "- Compute `Mean.Intensity` as the sum of replicate intensities divided by the count of non-zero replicates.",
                "- Preserve grouping-column order using categorical dtypes derived from the metadata file.",
                "- Fill unmapped proteins with `Protein.name = Others` and `Protein.group = Whey`.",
                "- Compute GRAVY scores from the Kyte-Doolittle scale.",
                "- `filterPeptides()` applies sequence filtering before regrouping, matching the upstream R implementation.",
                "- `ttestPeptides()` replaces zero intensities with the pseudocount before log2 transformation.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (artifacts_dir / "FIXTURE_CATALOG.md").write_text(
        "\n".join(
            [
                "# Fixture Catalog",
                "",
                "| File | Role |",
                "| --- | --- |",
                "| `tests/fixtures/v1.1.1-alpha/Yogurtexample_QR188-205.csv` | Canonical peptide input for parity and examples. |",
                "| `tests/fixtures/v1.1.1-alpha/Intensity_columns.csv` | Canonical grouping and replicate metadata. |",
                "| `tests/fixtures/v1.1.1-alpha/protein_mapping.csv` | Canonical protein-name/group lookup for parity tests. |",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (artifacts_dir / "API_MAPPING.md").write_text(
        "\n".join(
            [
                "# API Mapping",
                "",
                "| R export | Python surface | Notes |",
                "| --- | --- | --- |",
                "| `processPeptides` | `processPeptides` | Returns `PeptidomicsResult` instead of an R list. |",
                "| `filterPeptides` | `filterPeptides` | Accepts and returns `PeptidomicsResult`. |",
                "| `ttestPeptides` | `ttestPeptides` | Returns `dict[str, pandas.DataFrame]`. |",
                "| `plot_int` | `plot_int` | Returns plotnine `ggplot`. |",
                "| `plot_type_num` | `plot_type_num` | Returns plotnine `ggplot`. |",
                "| `plot_count` | `plot_count` | Deprecated wrapper around `plot_type_num`. |",
                "| `plot_length_distribution` | `plot_length_distribution` | Returns plotnine `ggplot`. |",
                "| `plot_gravy_vs_intensity` | `plot_gravy_vs_intensity` | Returns plotnine `ggplot`. |",
                "| `plot_pep_align` | `plot_pep_align` | Returns matplotlib `Figure`. |",
                "| `plot_cleavage_site` | `plot_cleavage_site` | Returns matplotlib `Figure`. |",
                "| `plot_volcano` | `plot_volcano` | Returns `dict[str, ggplot]`. |",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (artifacts_dir / "DEVIATION_LOG.md").write_text(
        "\n".join(
            [
                "# Deviation Log",
                "",
                "- Python package metadata uses `1.1.1a0` while migration naming uses `v1.1.1-alpha`.",
                "- `ttestPeptides(test_method=\"treat\")` is implemented as a pure-Python thresholded t approximation, not a limma empirical-Bayes port.",
                "- Plot backends change from ggplot2/grid objects to plotnine `ggplot` and matplotlib `Figure` objects.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (artifacts_dir / "PARITY_REPORT.md").write_text(
        "\n".join(
            [
                "# Parity Report",
                "",
                "- Status: upstream reference frozen and baseline outputs generated.",
                "- Core expected row counts frozen from R reference output.",
                "- `plain` statistics reference output generated for `G120_Y1` vs `I120_Y1`.",
                "- `treat` reference output generated for significance-agreement benchmarking.",
                "- Plot smoke labels frozen for density legends.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (artifacts_dir / "RELEASE_CHECKLIST.md").write_text(
        "\n".join(
            [
                "# Release Checklist",
                "",
                "- [ ] Core process parity checks pass.",
                "- [ ] Filter parity checks pass.",
                "- [ ] `plain` statistics parity checks pass.",
                "- [ ] `treat` significance agreement >= 95%.",
                "- [ ] Plot smoke tests pass.",
                "- [ ] README rewritten for Python usage.",
                "- [ ] Package metadata aligned to `1.1.1a0` / `v1.1.1-alpha` mapping.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summaries = []
    for path in sorted(reference_dir.rglob("*.csv.gz")):
        summaries.append(gzip_csv_summary(path))
    manifest = {
        "version_label": VERSION_LABEL,
        "description_version": "1.1.1.9000",
        "upstream_commit": commit_sha,
        "reference_tables": summaries,
        "plot_reference": [
            {
                "case": "plot_length_distribution_density",
                "colour_label": "Yogurt × Replicate",
            },
            {
                "case": "plot_gravy_vs_intensity_density",
                "colour_label": "Replicate",
            },
        ],
    }
    (reference_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
