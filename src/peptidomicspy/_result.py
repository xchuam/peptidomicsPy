from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(slots=True)
class PeptidomicsResult:
    dt_peptides: pd.DataFrame
    dt_peptides_int: pd.DataFrame
    dt_peptides_int_reps: pd.DataFrame
    dt_peptides_int_ttest: pd.DataFrame
    dt_peptides_typenum: pd.DataFrame
    dt_peptides_typenum_reps: pd.DataFrame
    dt_int_col: pd.DataFrame
    grp_cols: list[str]
    peptides_select_col_basic: list[str]

    _r_aliases = {
        "dt.peptides": "dt_peptides",
        "dt.peptides.int": "dt_peptides_int",
        "dt.peptides.int.reps": "dt_peptides_int_reps",
        "dt.peptides.int.ttest": "dt_peptides_int_ttest",
        "dt.peptides.typenum": "dt_peptides_typenum",
        "dt.peptides.typenum.reps": "dt_peptides_typenum_reps",
        "dt.int_col": "dt_int_col",
        "grp_cols": "grp_cols",
        "peptides_select_col.basic": "peptides_select_col_basic",
    }

    def __getitem__(self, key: str):
        attr = self._r_aliases.get(key, key)
        return getattr(self, attr)

    def copy(self) -> "PeptidomicsResult":
        return PeptidomicsResult(
            dt_peptides=self.dt_peptides.copy(deep=True),
            dt_peptides_int=self.dt_peptides_int.copy(deep=True),
            dt_peptides_int_reps=self.dt_peptides_int_reps.copy(deep=True),
            dt_peptides_int_ttest=self.dt_peptides_int_ttest.copy(deep=True),
            dt_peptides_typenum=self.dt_peptides_typenum.copy(deep=True),
            dt_peptides_typenum_reps=self.dt_peptides_typenum_reps.copy(deep=True),
            dt_int_col=self.dt_int_col.copy(deep=True),
            grp_cols=list(self.grp_cols),
            peptides_select_col_basic=list(self.peptides_select_col_basic),
        )

