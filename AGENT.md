# Repository Agent Guide

This repository is for migrating the R package `xchuam/peptidomicsR` into a native Python package.

Primary upstream reference:

- https://github.com/xchuam/peptidomicsR

The upstream R package describes itself as tools for peptidomics analysis of digesta from protein digestion and exposes functions for:

- data processing via `processPeptides()`
- filtering and statistics via `pcaPeptides()`, `filterPeptides()`, and `ttestPeptides()`
- visualization via `plot_int()`, `plot_count()`, `plot_length_distribution()`, `plot_cleavage_site()`, `plot_gravy_vs_intensity()`, `plot_pcaPeptides()`, and `plot_volcano()`

This repo should use both of these installed skills together:

- `cran-to-native-python-migration`
- `python-pypi-package-builder`

## Mission

Build a real Python package that preserves the useful external behavior of `peptidomicsR` while following modern Python packaging standards.

Default target:

- Python package name: `peptidomicspy`
- import package name: `peptidomicspy`

If later repo decisions require a different distribution name, update this file and all migration artifacts explicitly rather than drifting implicitly.

Version policy:

- The Python package version must align with the upstream R package version being migrated.
- Read the target version from the upstream R package metadata, normally `DESCRIPTION`.
- For this upstream package, treat R development versions ending in `.9000` as equivalent to the GitHub prerelease label `-alpha`.
- Example equivalence: `1.1.1.9000` in `DESCRIPTION` maps to `v1.1.1-alpha` in GitHub tagging and migration naming.
- Use the `alpha` form as the canonical human-facing migration label in this repository.
- Do not invent an independent Python package version line during parity migration.
- When branch names, migration artifacts, or status notes need a version label, prefer the `alpha` form rather than the `.9000` form.
- If Python package metadata needs a prerelease version string, keep it semantically aligned with the same upstream target and use a Python-packaging-valid prerelease representation when needed.
- If the Python package must temporarily diverge from the upstream R version, document the reason explicitly in the migration notes and release checklist.

Branch policy:

- Do not perform migration work directly on `main`.
- Maintain one long-lived core migration branch that tracks ongoing migration state across versions.
- Preferred core branch name: `migration/core`.
- For each upstream R package version, create a version-specific migration branch from the core migration branch.
- Preferred version branch format: `migration/v<r-package-version>`.
- For upstream `.9000` development versions, convert the branch name to the equivalent `alpha` label.
- Example: upstream `1.1.1.9000` should use branch `migration/v1.1.1-alpha`.
- Example stable branch: `migration/v0.3.0`.
- Use the core migration branch to accumulate migration artifacts, audit history, and shared migration progress.
- Use the version branch to capture the work and decisions specific to a single upstream version target.
- Treat `main` as the user-facing stable line for the Python package, and only merge package-ready results into `main`.

## Required Working Style

- Treat the R package as the behavioral specification first.
- Do not do blind line-by-line syntax translation.
- Default migration mode is `native rewrite`, but classify each subsystem as exactly one of:
  - `native rewrite`
  - `faithful port`
  - `temporary parity bridge`
- Preserve API parity, test parity, and data-model parity unless a deviation is explicitly documented.
- Do not claim parity without pytest coverage, fixtures, or both.
- Before each substantial edit, state which migration stage the work belongs to.
- When making a tradeoff, state whether it is driven by API parity, test parity, or data-model parity.

## Skills To Use

### 1. `cran-to-native-python-migration`

This skill is the source of truth for migration order, parity rules, reference freeze, classification, fixture capture, and staged execution.

Always apply it for:

- auditing the upstream R package
- freezing exported API and package behavior
- deciding rewrite versus faithful port versus parity bridge
- defining semantic rules around missing values, factors, indexing, attributes, grouped summaries, and plotting behavior
- creating or updating migration artifacts

### 2. `python-pypi-package-builder`

This companion skill is required once work reaches Python package construction, implementation packaging, and release-readiness.

Always apply it for:

- `pyproject.toml`
- `src/` layout decisions
- package metadata
- test and dev extras
- linting and typing setup
- CI and release workflow
- build, install, wheel, and sdist verification

For this repository, prefer:

- `src/` layout
- modern PyPA packaging
- pytest
- Ruff
- mypy when type coverage is being introduced

Use `setuptools` with `setuptools_scm` only if release/version flow is intentionally tied to git tags. Otherwise use the packaging decision rules from the skill and document the choice.

## Migration Stages

Follow these stages unless the user explicitly redirects the work:

1. Stage 1: reference freeze and audit
2. Stage 2: migration classification
3. Stage 3: Python package scaffold
4. Stage 4: semantic normalization and fixture capture
5. Stage 5: interface and test parity
6. Stage 6: implementation
7. Stage 7: docs sync
8. Stage 8: release readiness

Do not skip Stage 1 just because implementation seems obvious.

At the start of a migration cycle, confirm the target upstream R package version, normalize `.9000` to the equivalent `alpha` label when applicable, and verify that the current branch name uses that canonical label before making substantial edits.

## First Move For New Work

Before proposing implementation, inspect the upstream R package surface at minimum:

- `DESCRIPTION`
- `NAMESPACE`
- `R/`
- `man/`
- `test/` or `tests/`
- `data/`
- `README.md`
- `README.Rmd`

Also record when expected directories are absent. At the time this guide was written, the upstream repository visibly contains:

- `R/`
- `data/`
- `man/`
- `test/`
- `DESCRIPTION`
- `NAMESPACE`
- `README.md`
- `README.Rmd`

No assumption should be made that `vignettes/` or `src/` exist until verified from the upstream source.

## Required Migration Artifacts

This repository is the Python package repository itself, but migration artifacts are still valuable project history and may be tracked in git.

Create and maintain migration-only artifacts under:

- `.migration-artifacts/`

This directory is for agent workflow and audit records. These files may be committed on migration branches so the migration process stays reproducible and traceable.

Branch expectations for these files:

- `.migration-artifacts/` should live on `migration/core` and version-specific migration branches.
- These files should generally not be merged into `main`, because `main` is the user-facing package branch.
- When promoting code from a migration branch to `main`, keep package-facing code, tests, docs, and release files, but leave migration-only records behind unless the user explicitly wants some of them preserved on `main`.

Store migration artifacts there as the audit advances:

- `.migration-artifacts/R_PACKAGE_DOSSIER.md`
- `.migration-artifacts/EXPORTED_API.csv`
- `.migration-artifacts/DEPENDENCY_MAP.md`
- `.migration-artifacts/DATA_ASSET_INVENTORY.md`
- `.migration-artifacts/COMPILED_CODE_INVENTORY.md`
- `.migration-artifacts/MIGRATION_MODE_MAP.md`
- `.migration-artifacts/PARITY_PRIORITY_MATRIX.md`
- `.migration-artifacts/SEMANTIC_RULES.md`
- `.migration-artifacts/FIXTURE_CATALOG.md`
- `.migration-artifacts/API_MAPPING.md`
- `.migration-artifacts/DEVIATION_LOG.md`
- `.migration-artifacts/PARITY_REPORT.md`
- `.migration-artifacts/RELEASE_CHECKLIST.md`

If these artifacts do not exist yet, scaffold and then fill them deliberately rather than writing ad hoc notes.

Only package-facing files that belong to the distributed Python library should be carried forward into `main`.

## Package Design Guidance

The Python package should feel native in implementation but still recognizable to users of the R package.

Prefer:

- explicit tabular data handling
- deterministic data transformations
- separation between data processing, statistical analysis, and plotting layers
- small composable functions for filtering and summarization
- clear model objects or typed records only where they genuinely improve maintainability

Be especially careful with:

- missing-value semantics
- grouped summaries and replicate aggregation
- peptide and protein identifier mapping
- regex and subset filtering behavior
- statistical comparison defaults
- PCA preprocessing assumptions
- plot parameter defaults and naming

If Python behavior intentionally differs from R behavior, record it in `.migration-artifacts/DEVIATION_LOG.md` and reflect it in `.migration-artifacts/PARITY_REPORT.md`.

## Testing Expectations

- Convert test intent, not just test syntax.
- Prefer pytest for all new Python tests.
- Add parity-focused fixtures when behavior is ambiguous or data-dependent.
- For every function modified during the current migration process, add or update explicit R-versus-Python comparison coverage.
- For functions that return tables, vectors, dicts, dataclass-like objects, or other structured results, compare the Python outputs directly against the corresponding R outputs and confirm the same calculations, column content, row counts, ordering rules, and key metadata assumptions.
- For functions that return figures, compare Python renders against R renders using the same plotting inputs and review figure behavior as well as style details that matter for user interpretation, including labels, legends, panel order, axis scale, axis units, tick formatting, color use, and annotation placement.
- For statistical and plotting code, verify both numeric outputs and metadata assumptions where practical.
- Do not merge major migration steps without corresponding tests or documented fixture-based verification.
- When a previously migrated function is changed again, update its parity tests in the same change rather than leaving the comparison suite stale.

## Manual Validation Notebook

- Maintain a version-scoped Jupyter notebook under `tests/` for manual inspection of the current migration target.
- The notebook must include manual checking cells for every function modified in the current migration branch.
- The notebook should cover both structured-output functions and plotting functions.
- For structured-output functions, include cells that load the R reference outputs, run the Python implementation, and show side-by-side or assertion-backed comparisons that make mismatches easy to inspect.
- For plotting functions, include cells that render the R reference figure, the Python figure, and a visual comparison view so manual review is straightforward.
- When a function is added, ported, or behaviorally changed during migration, update the notebook in the same branch so manual checking stays complete.
- Prefer a version-specific notebook name tied to the active migration target, for example `tests/manual_validation_v1_1_1_alpha.ipynb` or an equivalent repo-approved naming convention.

## Packaging Expectations

By the time Stage 3 begins, the repo should become a real Python package and not just a notebook-style translation workspace.

Prefer:

- `pyproject.toml`
- `src/peptidomicspy/`
- `tests/`
- optional example data only when licensing and package size are acceptable

Packaging, metadata, CI, and release decisions must follow `python-pypi-package-builder`.

## Output Style For Future Agents

When reporting migration work, structure updates in this order when useful:

1. what was inspected
2. risks or semantic traps
3. recommended edits
4. exact files changed
5. next highest-value step

Keep decisions reusable for future agents. Favor explicit tables, checklists, and artifact updates over scattered prose.
