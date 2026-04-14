# Migration Mode Map

| Subsystem | Mode | Notes |
| --- | --- | --- |
| `processPeptides` | native rewrite | Preserve result schema and grouping semantics with pandas. |
| `filterPeptides` | faithful port | Preserve selector and regrouping behavior, including current sequence-first filtering semantics. |
| `ttestPeptides(plain)` | faithful port | Match selector parsing and Welch/Student behavior on log2 intensities. |
| `ttestPeptides(treat)` | native rewrite with documented deviation | Pure-Python thresholded t approximation; not limma parity. |
| Bar and density plots | native rewrite | Preserve parameter surface while switching to plotnine. |
| `plot_pep_align` | native rewrite | Matplotlib implementation for stacked peptide alignment layout. |
| `plot_cleavage_site` | native rewrite | Matplotlib plus logomaker implementation. |
