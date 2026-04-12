import subprocess
from pathlib import Path

REQUIRED_R_PACKAGES = (
    "affy",
    "AnnotationDbi",
    "hgu133plus2.db",
    "hgu133plus2cdf"
)


def ensure_bioconductor_packages(r_dir, packages = REQUIRED_R_PACKAGES, *, install_missing = True):
    deps_script_path = Path(r_dir) / "_deps.R"
    subprocess.run(
        [
            "Rscript",
            str(deps_script_path),
            str(install_missing).lower(),
            *packages
        ],
        check=True
    )


def run_rma_quality_control(raw_dir, artifacts_dir, r_dir):
    master_sheet_path = artifacts_dir / "master_sample_sheet.csv"
    quality_control_script_path = r_dir / "quality_control.R"
    expression_matrix_csv = artifacts_dir / "expression_matrix_rma.csv"
    preprocessing_qc_csv = artifacts_dir / "preprocessing_qc.csv"
    expression_bundle_rds = artifacts_dir / "expression_bundle_rma.rds"

    ensure_bioconductor_packages(r_dir=r_dir)

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
        check=True
    )

    return {
        "expression_matrix_csv": expression_matrix_csv,
        "preprocessing_qc_csv": preprocessing_qc_csv,
        "expression_bundle_rds": expression_bundle_rds
    }
