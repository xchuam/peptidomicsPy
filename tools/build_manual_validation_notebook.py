from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = ROOT / "tests" / "manual_validation_v1_1_1_alpha.ipynb"


def md(text: str):
    return nbf.v4.new_markdown_cell(dedent(text).strip() + "\n")


def code(text: str):
    return nbf.v4.new_code_cell(dedent(text).strip() + "\n")


cells = [
    md(
        """
        # Manual R-vs-Python validation for `v1.1.1-alpha`

        This notebook is for interactive checking of the migrated Python package against the upstream R package.

        What it does:
        - builds a live R parity suite from the upstream `peptidomicsR` checkout
        - focuses on the plotting functions currently under active parity review
        - shows three-panel figure comparisons for every selected plotting case

        Suggested use:
        1. Run the setup cells first.
        2. Run sections `7C`, `7D`, and `7E`.
        3. Compare the Python render against the R reference visually.
        """
    ),
    code(
        """
        %matplotlib inline

        import gc
        from pathlib import Path
        import os
        import shutil
        import subprocess
        import sys
        import tempfile
        import warnings

        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        import numpy as np
        import pandas as pd
        from IPython.display import display, Markdown
        """
    ),
    code(
        """
        def find_repo_root(start: Path) -> Path:
            start = start.resolve()
            for candidate in [start, *start.parents]:
                if (candidate / "pyproject.toml").exists():
                    return candidate
            raise FileNotFoundError("Could not find repo root from current working directory.")


        REPO_ROOT = find_repo_root(Path.cwd())
        SRC_DIR = REPO_ROOT / "src"
        TESTS_DIR = REPO_ROOT / "tests"
        FIXTURE_DIR = TESTS_DIR / "fixtures" / "v1.1.1-alpha"
        DEFAULT_UPSTREAM_R_DIR = Path(os.environ.get("PEPTIDOMICSR_UPSTREAM_DIR", "/tmp/peptidomicsR-v1.1.1-alpha"))
        R_PARITY_DIR = REPO_ROOT / ".tmp-manual-r-parity"
        NOTEBOOK_SMOKE_MODE = os.environ.get("PEPTIDOMICSPY_NOTEBOOK_SMOKE") == "1"

        for path in [SRC_DIR, TESTS_DIR]:
            if str(path) not in sys.path:
                sys.path.insert(0, str(path))

        from conftest import assert_dataframe_close, image_mae, load_rgb_image, resize_nn
        from parity_contract import ALL_PARITY_CASES, FUNCTION_PARAMETER_COVERAGE, PLOT_PARITY_CASES, STRUCTURED_PARITY_CASES
        from peptidomicspy import (
            calculate_gravy,
            filterPeptides,
            load_example_protein_mapping,
            plot_cleavage_site,
            plot_count,
            plot_gravy_vs_intensity,
            plot_int,
            plot_length_distribution,
            plot_pep_align,
            plot_type_num,
            plot_volcano,
            processPeptides,
            ttestPeptides,
        )
        """
    ),
    code(
        """
        from matplotlib.figure import Figure
        from plotnine import ggplot


        def run_r_parity_suite(force: bool = False) -> Path:
            if force and R_PARITY_DIR.exists():
                shutil.rmtree(R_PARITY_DIR)
            if R_PARITY_DIR.exists():
                return R_PARITY_DIR
            if not DEFAULT_UPSTREAM_R_DIR.exists():
                raise FileNotFoundError(
                    f"Upstream R checkout not found at {DEFAULT_UPSTREAM_R_DIR}. "
                    "Set PEPTIDOMICSR_UPSTREAM_DIR if needed."
                )
            env = os.environ.copy()
            env.setdefault("R_LIBS_USER", str(Path.home() / "R" / "x86_64-pc-linux-gnu-library" / "4.3"))
            env["PEPTIDOMICSPY_PARITY_SCOPE"] = "plot_focus"
            subprocess.run(
                [
                    "Rscript",
                    "tools/generate_r_parity_suite.R",
                    str(DEFAULT_UPSTREAM_R_DIR),
                    str(FIXTURE_DIR),
                    str(R_PARITY_DIR),
                ],
                cwd=REPO_ROOT,
                env=env,
                check=True,
            )
            return R_PARITY_DIR


        def read_csv(path: Path) -> pd.DataFrame:
            return pd.read_csv(path)


        CASE_OBSERVATIONS: dict[str, list[dict[str, object]]] = {}


        def record_case_observation(
            case_name: str,
            label: str,
            status: str,
            details: str = "",
            metric_name: str | None = None,
            metric_value: float | int | None = None,
        ) -> None:
            CASE_OBSERVATIONS.setdefault(case_name, []).append(
                {
                    "label": label,
                    "status": status,
                    "details": details,
                    "metric_name": metric_name,
                    "metric_value": metric_value,
                }
            )


        def summarize_case(case_name: str) -> dict[str, object]:
            observations = CASE_OBSERVATIONS.get(case_name, [])
            if not observations:
                return {"status": "NOT RUN", "details": "No comparison has been run yet.", "metrics": ""}
            statuses = {str(item["status"]) for item in observations}
            if "FAIL" in statuses:
                status = "FAIL"
            elif statuses == {"PASS"}:
                status = "PASS"
            elif statuses <= {"REVIEW"}:
                status = "REVIEW"
            else:
                status = "MIXED"
            details = "; ".join(
                f'{item["label"]}: {item["status"]}' + (f' ({item["details"]})' if item["details"] else "")
                for item in observations
            )
            metrics = "; ".join(
                f'{item["label"]} {item["metric_name"]}={item["metric_value"]:.4f}'
                for item in observations
                if item["metric_name"] is not None and item["metric_value"] is not None
            )
            return {"status": status, "details": details, "metrics": metrics}


        def compare_frame_case(
            case_name: str,
            actual: pd.DataFrame,
            reference: pd.DataFrame,
            head: int = 5,
            label: str = "frame",
        ) -> None:
            try:
                assert_dataframe_close(actual, reference)
                record_case_observation(
                    case_name,
                    label,
                    "PASS",
                    details=f"rows={len(actual)}, cols={len(actual.columns)}",
                )
                print(f"PASS: {case_name}::{label}")
            except AssertionError as exc:
                message = str(exc).splitlines()
                detail = message[0] if message else "DataFrame mismatch detected."
                record_case_observation(case_name, label, "FAIL", details=detail)
                print(f"FAIL: {case_name}::{label}")
                print(exc)
                print("Actual head:")
                display(actual.head(head))
                print("Reference head:")
                display(reference.head(head))


        def save_plot_object(obj, path: Path, width: float, height: float) -> None:
            if isinstance(obj, ggplot):
                obj.save(path, width=width, height=height, dpi=150, verbose=False)
                plt.close("all")
                return
            if isinstance(obj, Figure):
                obj.set_size_inches(width, height)
                obj.savefig(path, dpi=150)
                plt.close(obj)
                return
            if isinstance(obj, np.ndarray):
                plt.imsave(path, obj)
                return
            raise TypeError(f"Unsupported plot object type: {type(obj)!r}")


        def show_plot_comparison(case_name: str, obj, width: float, height: float) -> None:
            py_path = Path(tempfile.mkdtemp(prefix="peptidomicspy-manual-")) / f"{case_name}.png"
            save_plot_object(obj, py_path, width, height)
            r_path = R_PARITY_DIR / "plots" / f"{case_name}.png"
            if not r_path.exists():
                run_r_parity_suite(force=True)
            if not r_path.exists():
                raise FileNotFoundError(f"Missing R parity image for {case_name}: {r_path}")
            mae = image_mae(py_path, r_path)
            record_case_observation(case_name, "plot", "REVIEW", details=f"Visual comparison rendered. MAE={mae:.4f}", metric_name="mae", metric_value=mae)
            if NOTEBOOK_SMOKE_MODE:
                print(f"{case_name}: MAE={mae:.4f}")
                return
            py_img = load_rgb_image(py_path)
            r_img = load_rgb_image(r_path)
            h = max(py_img.shape[0], r_img.shape[0])
            w = max(py_img.shape[1], r_img.shape[1])
            py_img = resize_nn(py_img, h, w)
            r_img = resize_nn(r_img, h, w)
            diff = np.abs(py_img - r_img)

            fig, axes = plt.subplots(1, 3, figsize=(16, 5))
            fig.suptitle(f"{case_name}: {ALL_PARITY_CASES[case_name]}", fontsize=12)
            axes[0].imshow(r_img)
            axes[0].set_title(f"R reference\\n{case_name}")
            axes[1].imshow(py_img)
            axes[1].set_title(f"Python render\\nMAE={mae:.4f}")
            axes[2].imshow(diff)
            axes[2].set_title("Absolute difference")
            for axis in axes:
                axis.axis("off")
            plt.tight_layout(rect=(0, 0, 1, 0.94))
            display(fig)
            plt.close(fig)
            plt.close("all")
            gc.collect()


        def make_custom_pep_align_result(process_result):
            custom = process_result.copy()
            custom.dt_peptides_int_reps = custom.dt_peptides_int_reps.rename(
                columns={
                    "Sequence": "Seq.custom",
                    "Start.position": "Start.custom",
                    "End.position": "End.custom",
                    "Intensity": "Intensity.custom",
                }
            )
            return custom


        def show_function_coverage(function_names: list[str]) -> None:
            rows = []
            for function_name in function_names:
                for parameter_name, case_names in FUNCTION_PARAMETER_COVERAGE[function_name].items():
                    rows.append(
                        {
                            "function": function_name,
                            "parameter": "(call)" if parameter_name == "__call__" else parameter_name,
                            "cases": ", ".join(case_names),
                        }
                    )
            display(pd.DataFrame(rows))


        def show_case_catalog(function_names: list[str]) -> None:
            case_names: list[str] = []
            for function_name in function_names:
                for names in FUNCTION_PARAMETER_COVERAGE[function_name].values():
                    for case_name in names:
                        if case_name not in case_names:
                            case_names.append(case_name)
            rows = [
                {"case": case_name, "description": ALL_PARITY_CASES[case_name]}
                for case_name in case_names
            ]
            display(pd.DataFrame(rows))


        def show_case_results(function_names: list[str]) -> None:
            case_names: list[str] = []
            for function_name in function_names:
                for names in FUNCTION_PARAMETER_COVERAGE[function_name].values():
                    for case_name in names:
                        if case_name not in case_names:
                            case_names.append(case_name)
            rows = []
            for case_name in case_names:
                summary = summarize_case(case_name)
                rows.append(
                    {
                        "case": case_name,
                        "status": summary["status"],
                        "description": ALL_PARITY_CASES[case_name],
                        "details": summary["details"],
                        "metrics": summary["metrics"],
                    }
                )
            display(pd.DataFrame(rows))


        def show_parameter_results(function_names: list[str]) -> None:
            rows = []
            for function_name in function_names:
                for parameter_name, case_names in FUNCTION_PARAMETER_COVERAGE[function_name].items():
                    for case_name in case_names:
                        summary = summarize_case(case_name)
                        rows.append(
                            {
                                "function": function_name,
                                "parameter": "(call)" if parameter_name == "__call__" else parameter_name,
                                "case": case_name,
                                "status": summary["status"],
                                "metrics": summary["metrics"],
                                "details": summary["details"],
                            }
                        )
            display(pd.DataFrame(rows))


        def iter_case_items(case_map: dict[str, tuple], smoke_cases: list[str] | None = None):
            if NOTEBOOK_SMOKE_MODE and smoke_cases is not None:
                keep = set(smoke_cases)
                return [(name, spec) for name, spec in case_map.items() if name in keep]
            return list(case_map.items())
        """
    ),
    md(
        """
        ## 1. Build the live R parity suite

        Run the next cell first. It regenerates the comparison suite from the upstream R package into `.tmp-manual-r-parity/`.
        """
    ),
    code(
        """
        run_r_parity_suite(force=True)
        print(f"Live R parity suite: {R_PARITY_DIR}")
        """
    ),
    code(
        """
        PROCESS_RESULT = processPeptides(
            peptides_file=FIXTURE_DIR / "Yogurtexample_QR188-205.csv",
            intensity_columns_file=FIXTURE_DIR / "Intensity_columns.csv",
            protein_mapping_file=FIXTURE_DIR / "protein_mapping.csv",
        )
        CUSTOM_PEP_ALIGN_RESULT = make_custom_pep_align_result(PROCESS_RESULT)
        TTEST_PLAIN = ttestPeptides(PROCESS_RESULT, comparisons=[("G120_Y1", "I120_Y1")], test_method="plain")
        TREAT_KEY = "G120_Y1_vs_I120_Y1"
        TREAT_REFERENCE = read_csv(R_PARITY_DIR / "tables" / "ttest" / "treat_main.csv.gz")
        TREAT_REFERENCE_DICT = {TREAT_KEY: TREAT_REFERENCE}
        TREAT_LABELS = list(TREAT_REFERENCE["Sequence"].head(2))
        PLAIN_LABELS = list(TTEST_PLAIN[TREAT_KEY]["Sequence"].head(2))
        print("Prepared Python fixtures and comparison tables.")
        """
    ),
    md(
        """
        ## 2. Build the Python comparison inputs

        These objects are reused by the plotting sections below.
        """
    ),
    md(
        """
        ## 7. Plot parity

        Each plot cell below:
        - runs every named plot case against the R reference
        - shows one three-panel figure per case
        - uses figure titles to identify the compare setting being shown
        """
    ),
    md(
        """
        ### 7C. `plot_volcano()`
        """
    ),
    md(
        """
        `plot_volcano_plain`
        """
    ),
    code(
        """
        show_plot_comparison(
            "plot_volcano_plain",
            plot_volcano(TTEST_PLAIN, comparisons=[("G120_Y1", "I120_Y1")], test_method="plain")[TREAT_KEY],
            8,
            5,
        )
        """
    ),
    md(
        """
        `plot_volcano_treat_annotated`
        """
    ),
    code(
        """
        show_plot_comparison(
            "plot_volcano_treat_annotated",
            plot_volcano(
                TREAT_REFERENCE_DICT,
                comparisons=[("G120_Y1", "I120_Y1")],
                test_method="treat",
                label_seqs=TREAT_LABELS,
                highlight_seqs={"black": [TREAT_LABELS[0]]},
            )[TREAT_KEY],
            8,
            5,
        )
        """
    ),
    md(
        """
        `plot_volcano_custom_ylog`
        """
    ),
    code(
        """
        show_plot_comparison(
            "plot_volcano_custom_ylog",
            plot_volcano(
                TTEST_PLAIN,
                comparisons=1,
                test_method="plain",
                show_threshold=True,
                lfc_thresh=0.5,
                alpha=0.1,
                fill_values={"no": "#AAAAAA", "yes": "#E41A1C"},
                point_size=2.5,
                point_alpha=0.6,
                label_seqs=PLAIN_LABELS,
                label_size=4,
                label_col="navy",
                highlight_seqs={"black": [PLAIN_LABELS[0]], "#377EB8": [PLAIN_LABELS[1]]},
                highlight_size=3,
                highlight_stroke=2,
                y_log_scale=True,
            )[TREAT_KEY],
            8,
            5,
        )
        """
    ),
    md(
        """
        ### 7D. `plot_pep_align()`

        This section is the main manual check for the range, custom-column, and save-path parameters.
        Pay special attention to panel order, plotted peptide set, x/y range clipping, label placement, and saved-output behavior.
        """
    ),
    md(
        """
        `plot_pep_align_lookup_default`
        """
    ),
    code(
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            show_plot_comparison(
                "plot_pep_align_lookup_default",
                plot_pep_align(PROCESS_RESULT, protein_name="P02662"),
                8,
                5,
            )
        """
    ),
    md(
        """
        `plot_pep_align_mean`
        """
    ),
    code(
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            show_plot_comparison(
                "plot_pep_align_mean",
                plot_pep_align(PROCESS_RESULT, protein_name="P81265", protein_seq="A" * 650, auto_size=False),
                12,
                8,
            )
        """
    ),
    md(
        """
        `plot_pep_align_ranged_labels`
        """
    ),
    code(
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            show_plot_comparison(
                "plot_pep_align_ranged_labels",
                plot_pep_align(
                    PROCESS_RESULT,
                    protein_name="P02662",
                    filter_params={"Digest.stage": "G120"},
                    x_interval=5,
                    x_range=(80, 165),
                    y_range=(0, 12),
                    label_seq={
                        "highlight_1": ["SVEQKHIQKE", "VEQKHIQKEDVPSERY"],
                        "highlight_2": ["QKHIQKEDVPSERY"],
                    },
                    label_col={"highlight_1": "red", "highlight_2": "blue"},
                ),
                9,
                6,
            )
        """
    ),
    md(
        """
        `plot_pep_align_reps_custom_cols`
        """
    ),
    code(
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            show_plot_comparison(
                "plot_pep_align_reps_custom_cols",
                plot_pep_align(
                    CUSTOM_PEP_ALIGN_RESULT,
                    protein_name="P02662",
                    type="reps",
                    filter_params={"Yogurt": "Y1"},
                    seq_col="Seq.custom",
                    start_col="Start.custom",
                    end_col="End.custom",
                    intensity_col="Intensity.custom",
                ),
                9,
                6,
            )
        """
    ),
    md(
        """
        `plot_pep_align_saved`
        """
    ),
    code(
        """
        if NOTEBOOK_SMOKE_MODE:
            print("Smoke mode: skipping the high-cost saved-file parity render for plot_pep_align().")
            record_case_observation(
                "plot_pep_align_saved",
                "saved_output",
                "NOT RUN",
                details="Skipped in smoke mode. Run this cell normally for the saved-file comparison.",
            )
        else:
            saved_path = Path(tempfile.mkdtemp(prefix="peptidomicspy-saved-")) / "plot_pep_align_saved.png"
            saved_plot = plot_pep_align(
                PROCESS_RESULT,
                protein_name="P02662",
                filter_params={"Digest.stage": "G120"},
                auto_size=True,
                save_file_location=str(saved_path),
                save_file_dpi=200,
                width_per_aa=0.12,
                height_per_row=0.2,
                min_width=5,
                min_height=4,
            )
            print("saved file exists:", saved_path.exists())
            print("suggested size:", getattr(saved_plot, "suggested_fig_size", None))
            print("save attrs:", getattr(saved_plot, "save_file_location", None), getattr(saved_plot, "save_file_dpi", None))
            record_case_observation(
                "plot_pep_align_saved",
                "saved_output",
                "PASS" if saved_path.exists() else "FAIL",
                details=f"saved={saved_path.exists()}, suggested_fig_size={getattr(saved_plot, 'suggested_fig_size', None)}",
            )
            plt.close(saved_plot)
        """
    ),
    code(
        """
        if not NOTEBOOK_SMOKE_MODE:
            show_plot_comparison("plot_pep_align_saved", mpimg.imread(saved_path), 8, 5)
        """
    ),
    md(
        """
        ### 7E. `plot_cleavage_site()`
        """
    ),
    md(
        """
        `plot_cleavage_N`
        """
    ),
    code(
        """
        show_plot_comparison(
            "plot_cleavage_N",
            plot_cleavage_site(PROCESS_RESULT, terminal="N"),
            8,
            5,
        )
        """
    ),
    md(
        """
        `plot_cleavage_both_filtered`
        """
    ),
    code(
        """
        show_plot_comparison(
            "plot_cleavage_both_filtered",
            plot_cleavage_site(
                PROCESS_RESULT,
                terminal="both",
                measure="type_num",
                replicate_mode="reps",
                filter_params={"Yogurt": "Y1"},
            ),
            12,
            5,
        )
        """
    ),
    md(
        """
        `plot_cleavage_C_keep_constant_linear`
        """
    ),
    code(
        """
        show_plot_comparison(
            "plot_cleavage_C_keep_constant_linear",
            plot_cleavage_site(
                PROCESS_RESULT,
                terminal="C",
                filter_params={"Yogurt": "Y1", "Digest.stage": "G120"},
                scientific_10_y=False,
                drop_constant_groups=False,
            ),
            8,
            5,
        )
        """
    ),
    md(
        """
        ## 8. Final status

        This notebook is now focused on the three plotting sections under active parity review:
        - `7C` `plot_volcano()`
        - `7D` `plot_pep_align()`
        - `7E` `plot_cleavage_site()`
        """
    ),
]


notebook = nbf.v4.new_notebook()
notebook.cells = cells
notebook.metadata = {
    "kernelspec": {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    },
    "language_info": {
        "name": "python",
        "version": "3.11",
    },
}

NOTEBOOK_PATH.write_text(nbf.writes(notebook), encoding="utf-8")
print(f"Wrote {NOTEBOOK_PATH}")
