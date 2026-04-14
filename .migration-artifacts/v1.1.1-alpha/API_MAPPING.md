# API Mapping

| R export | Python surface | Notes |
| --- | --- | --- |
| `processPeptides` | `processPeptides` | Returns `PeptidomicsResult` instead of an R list. |
| `filterPeptides` | `filterPeptides` | Accepts and returns `PeptidomicsResult`. |
| `ttestPeptides` | `ttestPeptides` | Returns `dict[str, pandas.DataFrame]`. |
| `plot_int` | `plot_int` | Returns plotnine `ggplot`. |
| `plot_type_num` | `plot_type_num` | Returns plotnine `ggplot`. |
| `plot_count` | `plot_count` | Deprecated wrapper around `plot_type_num`. |
| `plot_length_distribution` | `plot_length_distribution` | Returns plotnine `ggplot`. |
| `plot_gravy_vs_intensity` | `plot_gravy_vs_intensity` | Returns plotnine `ggplot`. |
| `plot_pep_align` | `plot_pep_align` | Returns matplotlib `Figure`. |
| `plot_cleavage_site` | `plot_cleavage_site` | Returns matplotlib `Figure`. |
| `plot_volcano` | `plot_volcano` | Returns `dict[str, ggplot]`. |
