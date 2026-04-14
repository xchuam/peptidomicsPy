from __future__ import annotations

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")

HYDROPATHY = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

PROTEIN_NAME_COLOR = {
    "Beta-casein": "#E41A1C",
    "Kappa-casein": "#377EB8",
    "Alpha-S1-casein": "#4DAF4A",
    "Alpha-S2-casein": "#984EA3",
    "Alpha-lactalbumin": "#FF7F00",
    "Beta-lactoglobulin": "#FCC88F",
    "Others": "#666666",
}

PROTEIN_GROUP_COLOR = {
    "Casein": "#E8A03F",
    "Whey": "#0073B4",
}

PROTEIN_COLOR = {**PROTEIN_NAME_COLOR, **PROTEIN_GROUP_COLOR}
DEFAULT_COLOR = "#9AC9DB"

PEPTIDES_SELECT_COL_BASIC = [
    "Sequence",
    "Leading.razor.protein",
    "Length",
    "Start.position",
    "End.position",
    "Amino.acid.before",
    "First.amino.acid",
    "Last.amino.acid",
    "Amino.acid.after",
]
