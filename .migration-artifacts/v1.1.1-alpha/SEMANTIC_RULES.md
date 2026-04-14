# Semantic Rules

- Replace spaces with dots in all peptide column names before downstream processing.
- Parse the local fixture CSVs as semicolon-delimited UTF-8 with BOM handling.
- Remove peptides whose `Leading.razor.protein` starts with `CON_` or `REV_`.
- Use `Intensity.column` order from the metadata file to define `R1..Rn` replicate columns.
- Keep peptides in mean tables only if at least half of replicates are strictly greater than zero.
- Compute `Mean.Intensity` as the sum of replicate intensities divided by the count of non-zero replicates.
- Preserve grouping-column order using categorical dtypes derived from the metadata file.
- Fill unmapped proteins with `Protein.name = Others` and `Protein.group = Whey`.
- Compute GRAVY scores from the Kyte-Doolittle scale.
- `filterPeptides()` applies sequence filtering before regrouping, matching the upstream R implementation.
- `ttestPeptides()` replaces zero intensities with the pseudocount before log2 transformation.
