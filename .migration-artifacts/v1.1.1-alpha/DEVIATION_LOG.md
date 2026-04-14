# Deviation Log

- Python package metadata uses `1.1.1a0` while migration naming uses `v1.1.1-alpha`.
- `ttestPeptides(test_method="treat")` is implemented as a pure-Python thresholded t approximation, not a limma empirical-Bayes port.
- Plot backends change from ggplot2/grid objects to plotnine `ggplot` and matplotlib `Figure` objects.
