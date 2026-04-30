import subprocess

from project_tools._paths import INTERNAL_CONFIGURATION
from project_tools.expression_preprocessing import _r_environment, ensure_bioconductor_packages


def run_pathway_analysis(path_configuration):
    ensure_bioconductor_packages(path_configuration)

    if not path_configuration.limma_results_dir.exists():
        raise FileNotFoundError(
            f"Missing limma results directory: {path_configuration.limma_results_dir}. "
            "Run limma differential expression first."
        )

    if not any(path_configuration.limma_results_dir.glob("*_ranked.rnk")):
        raise FileNotFoundError(
            f"No ranked limma files found in: {path_configuration.limma_results_dir}. "
            "Run limma differential expression first."
        )

    path_configuration.pathway_results_dir.mkdir(parents=True, exist_ok=True)

    script_path = path_configuration.r_dir / "pathway_analysis.R"
    subprocess.run(
        [
            "Rscript",
            str(script_path),
            str(path_configuration.limma_results_dir),
            str(path_configuration.pathway_results_dir),
            str(path_configuration.pathway_summary_csv),
        ],
        env=_r_environment(path_configuration),
        check=True,
    )

    return {
        "status": "ok",
        "mode": "public",
        "pathway_results_dir": path_configuration.pathway_results_dir,
        "pathway_summary_csv": path_configuration.pathway_summary_csv,
    }


if __name__ == "__main__":
    run_pathway_analysis(INTERNAL_CONFIGURATION)
