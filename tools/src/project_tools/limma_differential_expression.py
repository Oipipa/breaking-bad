import subprocess

from project_tools._paths import INTERNAL_CONFIGURATION
from project_tools.expression_preprocessing import _r_environment, ensure_bioconductor_packages


def run_limma_differential_expression(path_configuration):
    ensure_bioconductor_packages(path_configuration)

    if not path_configuration.master_sample_sheet_path.exists():
        raise FileNotFoundError(
            f"Missing master sample sheet: {path_configuration.master_sample_sheet_path}. "
            "Run sample sheet generation first."
        )

    if not path_configuration.expression_bundle_rds.exists():
        raise FileNotFoundError(
            f"Missing expression bundle: {path_configuration.expression_bundle_rds}. "
            "Run RMA preprocessing first."
        )

    path_configuration.limma_results_dir.mkdir(parents=True, exist_ok=True)

    script_path = path_configuration.r_dir / "limma_secondary_outcomes.R"
    subprocess.run(
        [
            "Rscript",
            str(script_path),
            str(path_configuration.master_sample_sheet_path),
            str(path_configuration.expression_bundle_rds),
            str(path_configuration.limma_results_dir),
            str(path_configuration.limma_comparison_summary_csv),
        ],
        env=_r_environment(path_configuration),
        check=True,
    )

    return {
        "status": "ok",
        "mode": "public",
        "limma_results_dir": path_configuration.limma_results_dir,
        "limma_comparison_summary_csv": path_configuration.limma_comparison_summary_csv,
    }


if __name__ == "__main__":
    run_limma_differential_expression(INTERNAL_CONFIGURATION)
