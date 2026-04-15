from __future__ import annotations

import inspect
import re
import subprocess
from pathlib import Path

import pytest
from matplotlib.colors import to_hex

import peptidomicspy

from parity_contract import (
    ALL_PARITY_CASES,
    FUNCTION_PARAMETER_COVERAGE,
    PYTHON_ONLY_PARITY_HELPERS,
    R_FUNCTIONS_WITH_DIRECT_COUNTERPARTS,
)


_REQUIRED = object()
_R_NAMED_STRING_VECTOR = re.compile(r'([A-Za-z0-9_.]+)\s*=\s*"([^"]*)"')


def _normalize_color_string(value: str) -> str:
    lowered = value.lower()
    match = re.fullmatch(r"gr[ae]y(\d{1,3})", lowered)
    if match:
        amount = max(0, min(100, int(match.group(1)))) / 100
        grey = int(round(amount * 255))
        return f"#{grey:02x}{grey:02x}{grey:02x}"
    return to_hex(value, keep_alpha=False)


def _load_r_formals(upstream_r_dir: Path) -> dict[str, list[tuple[str, str]]]:
    script = """
args <- commandArgs(trailingOnly = TRUE)
pkg <- normalizePath(args[[1]], mustWork = TRUE)
env <- new.env(parent = baseenv())
for (f in list.files(file.path(pkg, "R"), pattern = "\\\\.R$", full.names = TRUE)) {
  sys.source(f, envir = env)
}
funs <- c(%s)
for (nm in funs) {
  fm <- formals(get(nm, envir = env))
  for (p in names(fm)) {
    raw <- paste(deparse(fm[[p]]), collapse = " ")
    if (identical(raw, "name")) raw <- ""
    cat(nm, "\\t", p, "\\t", raw, "\\n", sep = "")
  }
}
""" % ", ".join(f'"{name}"' for name in R_FUNCTIONS_WITH_DIRECT_COUNTERPARTS)
    completed = subprocess.run(
        ["Rscript", "-e", script, str(upstream_r_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    out: dict[str, list[tuple[str, str]]] = {name: [] for name in R_FUNCTIONS_WITH_DIRECT_COUNTERPARTS}
    for raw_line in completed.stdout.splitlines():
        function_name, param_name, raw_default = raw_line.split("\t", 2)
        out[function_name].append((param_name, raw_default.strip()))
    return out


def _normalize_r_default(raw: str):
    if raw == "":
        return _REQUIRED
    if raw == "NULL":
        return None
    if raw == "TRUE":
        return True
    if raw == "FALSE":
        return False
    if raw.startswith('c("') and raw.endswith(")"):
        choices = re.findall(r'"([^"]+)"', raw)
        if choices:
            return choices[0]
    if raw.startswith('"') and raw.endswith('"'):
        return raw[1:-1]
    named_matches = _R_NAMED_STRING_VECTOR.findall(raw)
    if named_matches:
        return {
            key: _normalize_color_string(value)
            if value.startswith("#") or value.isalpha() or value.lower().startswith("grey")
            else value
            for key, value in named_matches
        }
    try:
        numeric = float(raw)
    except ValueError:
        return raw
    if numeric.is_integer():
        return int(numeric)
    return numeric


def _normalize_python_default(function_name: str, param_name: str, default):
    if default is inspect._empty:
        return _REQUIRED
    if function_name == "plot_volcano" and param_name == "fill_values":
        return {"no": _normalize_color_string("grey75"), "yes": _normalize_color_string("#FFC010")}
    return default


@pytest.fixture(scope="session")
def r_formals(upstream_r_dir: Path) -> dict[str, list[tuple[str, str]]]:
    return _load_r_formals(upstream_r_dir)


@pytest.mark.parametrize("function_name", R_FUNCTIONS_WITH_DIRECT_COUNTERPARTS)
def test_python_parameter_order_matches_r_formals(function_name: str, r_formals: dict[str, list[tuple[str, str]]]) -> None:
    python_function = getattr(peptidomicspy, function_name)
    python_params = list(inspect.signature(python_function).parameters)
    r_params = [name for name, _raw_default in r_formals[function_name]]
    assert python_params == r_params


@pytest.mark.parametrize("function_name", R_FUNCTIONS_WITH_DIRECT_COUNTERPARTS)
def test_python_effective_defaults_match_r(function_name: str, r_formals: dict[str, list[tuple[str, str]]]) -> None:
    python_signature = inspect.signature(getattr(peptidomicspy, function_name))
    for param_name, raw_default in r_formals[function_name]:
        expected = _normalize_r_default(raw_default)
        actual = _normalize_python_default(function_name, param_name, python_signature.parameters[param_name].default)
        assert actual == expected


def test_python_only_helper_has_parity_target_documented() -> None:
    for function_name, note in PYTHON_ONLY_PARITY_HELPERS.items():
        assert hasattr(peptidomicspy, function_name)
        assert note


def test_parameter_coverage_matrix_references_known_cases() -> None:
    for function_name, parameter_map in FUNCTION_PARAMETER_COVERAGE.items():
        assert hasattr(peptidomicspy, function_name)
        signature = inspect.signature(getattr(peptidomicspy, function_name))
        if "__call__" in parameter_map:
            assert len(signature.parameters) == 0
        for param_name, case_names in parameter_map.items():
            if param_name != "__call__":
                assert param_name in signature.parameters
            assert case_names
            for case_name in case_names:
                assert case_name in ALL_PARITY_CASES


def test_every_public_parameter_has_named_parity_coverage() -> None:
    functions_to_check = set(R_FUNCTIONS_WITH_DIRECT_COUNTERPARTS) | set(PYTHON_ONLY_PARITY_HELPERS)
    for function_name in functions_to_check:
        signature = inspect.signature(getattr(peptidomicspy, function_name))
        parameter_map = FUNCTION_PARAMETER_COVERAGE[function_name]
        if len(signature.parameters) == 0:
            assert "__call__" in parameter_map
            continue
        covered = {name for name in parameter_map if name != "__call__"}
        assert covered == set(signature.parameters)
