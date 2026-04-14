# Parity Priority Matrix

| Priority | Surface | Reason |
| --- | --- | --- |
| P0 | `processPeptides` result tables and group ordering | All downstream API depends on this shape. |
| P0 | `filterPeptides` sequence and grouping semantics | Drives downstream subset workflows. |
| P0 | `ttestPeptides(plain)` | Numeric parity check for core statistics. |
| P1 | `ttestPeptides(treat)` significance agreement | Documented approximation target. |
| P1 | Plot return types and default labels | Required for user workflow continuity. |
| P2 | Plot visual styling | Backend changes are acceptable if semantics remain aligned. |
