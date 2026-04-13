import subprocess
import os
from project_tools._paths import INTERNAL_CONFIGURATION

REQUIRED_R_PACKAGES = (
    "affy",
    "AnnotationDbi",
    "hgu133plus2.db",
    "hgu133plus2cdf"
)


def _r_environment(path_configuration):
    r_library_dir = path_configuration.r_library_dir
    r_library_dir.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["R_LIBS_USER"] = str(r_library_dir)
    return environment


def ensure_bioconductor_packages(path_configuration):
    deps_script_path = path_configuration.r_dir / "_deps.R"
    subprocess.run(
        [
            "Rscript",
            str(deps_script_path),
            "true",
            *REQUIRED_R_PACKAGES
        ],
        env=_r_environment(path_configuration),
        check=True
    )


def run_rma_preprocessing(path_configuration):
    raw_dir = path_configuration.raw_dir
    artifacts_dir = path_configuration.artifacts_dir
    r_dir = path_configuration.r_dir
    master_sheet_path = artifacts_dir / "master_sample_sheet.csv"
    quality_control_script_path = r_dir / "quality_control.R"
    expression_matrix_csv = artifacts_dir / "expression_matrix_rma.csv"
    preprocessing_qc_csv = artifacts_dir / "preprocessing_qc.csv"
    expression_bundle_rds = artifacts_dir / "expression_bundle_rma.rds"

    ensure_bioconductor_packages(path_configuration)

    artifacts_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "Rscript",
            str(quality_control_script_path),
            str(master_sheet_path),
            str(raw_dir),
            str(expression_matrix_csv),
            str(preprocessing_qc_csv),
            str(expression_bundle_rds)
        ],
        env=_r_environment(path_configuration),
        check=True
    )

    return {
        "expression_matrix_csv": expression_matrix_csv,
        "preprocessing_qc_csv": preprocessing_qc_csv,
        "expression_bundle_rds": expression_bundle_rds
    }


if __name__ == "__main__":
    run_rma_preprocessing(INTERNAL_CONFIGURATION)
