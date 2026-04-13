import csv
import random
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from project_tools._paths import INTERNAL_CONFIGURATION

matplotlib.use("Agg")

MAX_PLOT_PROBES = 5000
TIMEPOINT_MARKERS = {
    "baseline": "o",
    "week8": "^",
    "week12": "s"
}
TIMEPOINT_LABELS = {
    "baseline": "Baseline",
    "week8": "Week 8",
    "week12": "Week 12"
}

def _read_master_rows(path_configuration):
    master_sheet_path = path_configuration.artifacts_dir / "master_sample_sheet.csv"
    with master_sheet_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _read_preprocessing_qc_rows(path_configuration):
    preprocessing_qc_csv = path_configuration.artifacts_dir / "preprocessing_qc.csv"
    with preprocessing_qc_csv.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _read_expression_matrix(path_configuration, max_probes=None, seed=0):
    expression_matrix_csv = path_configuration.artifacts_dir / "expression_matrix_rma.csv"
    rng = random.Random(seed)
    with expression_matrix_csv.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        sample_ids = header[4:]
        values = []
        total_probes = 0
        for row in reader:
            total_probes += 1
            vector = list(map(float, row[4:]))
            if max_probes is None:
                values.append(vector)
                continue
            if len(values) < max_probes:
                values.append(vector)
            else:
                index = rng.randrange(total_probes)
                if index < max_probes:
                    values[index] = vector
    return sample_ids, np.asarray(values, dtype=float), total_probes


def _build_qc_summary(master_rows, preprocessing_qc_rows, sample_ids, expression_matrix):
    master_by_sample = {row["geo_sample_id"]: row for row in master_rows}
    preprocessing_qc_by_sample = {
        row["geo_sample_id"]: float(row["present_pct"]) for row in preprocessing_qc_rows
    }

    centered_matrix = expression_matrix.T - expression_matrix.T.mean(axis=0, keepdims=True)
    left_singular_vectors, singular_values, _ = np.linalg.svd(centered_matrix, full_matrices=False)
    scores = left_singular_vectors * singular_values
    variance_pct = (singular_values ** 2) / np.sum(singular_values ** 2) * 100

    rma_median = np.median(expression_matrix, axis=0)
    rma_iqr = np.percentile(expression_matrix, 75, axis=0) - np.percentile(expression_matrix, 25, axis=0)

    qc_summary = []
    for index, sample_id in enumerate(sample_ids):
        row = dict(master_by_sample[sample_id])
        row["present_pct"] = preprocessing_qc_by_sample[sample_id]
        row["rma_median"] = float(rma_median[index])
        row["rma_iqr"] = float(rma_iqr[index])
        row["pc1"] = float(scores[index, 0])
        row["pc2"] = float(scores[index, 1])
        row["pc1_variance_pct"] = float(variance_pct[0])
        row["pc2_variance_pct"] = float(variance_pct[1])
        qc_summary.append(row)

    return qc_summary


def _write_qc_summary(path_configuration, qc_summary):
    qc_summary_csv = path_configuration.artifacts_dir / "qc_summary.csv"
    fieldnames = list(qc_summary[0].keys())
    with qc_summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(qc_summary)


def _write_qc_metrics(path_configuration, qc_summary, total_probes):
    qc_metrics_csv = path_configuration.artifacts_dir / "qc_metrics.csv"
    present_pct = np.asarray([row["present_pct"] for row in qc_summary], dtype=float)
    rma_median = np.asarray([row["rma_median"] for row in qc_summary], dtype=float)
    rma_iqr = np.asarray([row["rma_iqr"] for row in qc_summary], dtype=float)
    metrics = [
        ("n_samples", len(qc_summary)),
        ("n_subjects", len({row["subject_id"] for row in qc_summary})),
        ("n_probe_sets", total_probes),
        ("n_topiramate", sum(row["treatment_arm"] == "Topiramate" for row in qc_summary)),
        ("n_placebo", sum(row["treatment_arm"] == "Placebo" for row in qc_summary)),
        ("n_baseline", sum(row["timepoint"] == "baseline" for row in qc_summary)),
        ("n_week8", sum(row["timepoint"] == "week8" for row in qc_summary)),
        ("n_week12", sum(row["timepoint"] == "week12" for row in qc_summary)),
        ("present_pct_min", float(np.min(present_pct))),
        ("present_pct_median", float(np.median(present_pct))),
        ("present_pct_max", float(np.max(present_pct))),
        ("rma_median_min", float(np.min(rma_median))),
        ("rma_median_median", float(np.median(rma_median))),
        ("rma_median_max", float(np.max(rma_median))),
        ("rma_iqr_min", float(np.min(rma_iqr))),
        ("rma_iqr_median", float(np.median(rma_iqr))),
        ("rma_iqr_max", float(np.max(rma_iqr))),
        ("pc1_variance_pct", float(qc_summary[0]["pc1_variance_pct"])),
        ("pc2_variance_pct", float(qc_summary[0]["pc2_variance_pct"])),
    ]

    with qc_metrics_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        for metric, value in metrics:
            writer.writerow({"metric": metric, "value": value})


def _save_pca_plot(path_configuration, qc_summary):
    pca_plot_png = path_configuration.artifacts_dir / "pca_rma.png"
    fig, ax = plt.subplots(figsize=(10, 8))
    for treatment_arm in ("Topiramate", "Placebo"):
        for timepoint in TIMEPOINT_MARKERS:
            rows = [
                row for row in qc_summary
                if row["treatment_arm"] == treatment_arm and row["timepoint"] == timepoint
            ]
            if not rows:
                continue
            ax.scatter(
                [row["pc1"] for row in rows],
                [row["pc2"] for row in rows],
                marker=TIMEPOINT_MARKERS[timepoint],
                label=f"{treatment_arm} {TIMEPOINT_LABELS[timepoint]}"
            )

    ax.set_xlabel(f"PC1 ({qc_summary[0]['pc1_variance_pct']:.2f}%)")
    ax.set_ylabel(f"PC2 ({qc_summary[0]['pc2_variance_pct']:.2f}%)")
    ax.set_title("PCA")
    ax.legend(frameon=False, ncol=2)

    fig.tight_layout()
    fig.savefig(pca_plot_png, dpi=150)
    plt.close(fig)


def _save_boxplot(path_configuration, qc_summary, sample_ids, expression_matrix):
    rma_boxplot_png = path_configuration.artifacts_dir / "rma_boxplot.png"
    sample_index = {sample_id: index for index, sample_id in enumerate(sample_ids)}
    panel_order = [
        ("Topiramate", "baseline"),
        ("Topiramate", "week8"),
        ("Topiramate", "week12"),
        ("Placebo", "baseline"),
        ("Placebo", "week8"),
        ("Placebo", "week12")
    ]
    labels = [f"{arm}\n{TIMEPOINT_LABELS[tp]}" for arm, tp in panel_order]
    group_medians = []
    group_iqrs = []

    for panel in panel_order:
        treatment_arm, timepoint = panel
        panel_rows = [
            row for row in qc_summary
            if row["treatment_arm"] == treatment_arm and row["timepoint"] == timepoint
        ]
        panel_indices = [sample_index[row["geo_sample_id"]] for row in panel_rows]
        panel_matrix = expression_matrix[:, panel_indices]
        group_medians.append(np.median(panel_matrix, axis=0))
        group_iqrs.append(
            np.percentile(panel_matrix, 75, axis=0) - np.percentile(panel_matrix, 25, axis=0)
        )

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    axes[0].boxplot(group_medians, showfliers=False)
    axes[0].set_title("RMA Median (per sample)")
    axes[0].set_ylabel("log2 expression")

    axes[1].boxplot(group_iqrs, showfliers=False)
    axes[1].set_title("RMA IQR (per sample)")
    axes[1].set_ylabel("log2 expression")
    axes[1].set_xticks(range(1, len(labels) + 1))
    axes[1].set_xticklabels(labels, rotation=0)

    fig.suptitle("RMA Boxplots by Group", y=0.98)
    fig.tight_layout()
    fig.savefig(rma_boxplot_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_density_plot(path_configuration, qc_summary, sample_ids, expression_matrix):
    rma_density_png = path_configuration.artifacts_dir / "rma_density.png"
    sample_index = {sample_id: index for index, sample_id in enumerate(sample_ids)}
    panel_order = [
        ("Topiramate", "baseline"),
        ("Topiramate", "week8"),
        ("Topiramate", "week12"),
        ("Placebo", "baseline"),
        ("Placebo", "week8"),
        ("Placebo", "week12")
    ]
    bins = np.linspace(np.min(expression_matrix), np.max(expression_matrix), 200)
    centers = (bins[:-1] + bins[1:]) / 2

    fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharex=True, sharey=True)
    axes = axes.ravel()

    for ax, panel in zip(axes, panel_order):
        treatment_arm, timepoint = panel
        panel_rows = [
            row for row in qc_summary
            if row["treatment_arm"] == treatment_arm and row["timepoint"] == timepoint
        ]
        panel_indices = [sample_index[row["geo_sample_id"]] for row in panel_rows]
        panel_matrix = expression_matrix[:, panel_indices]
        density_curves = []
        for index in range(panel_matrix.shape[1]):
            density, edges = np.histogram(panel_matrix[:, index], bins=bins, density=True)
            density_curves.append(density)

        density_curves = np.asarray(density_curves)
        ax.plot(centers, density_curves.mean(axis=0))
        ax.set_title(f"{treatment_arm} {TIMEPOINT_LABELS[timepoint]}")

    axes[0].set_ylabel("Density")
    axes[3].set_ylabel("Density")
    for ax in axes[3:]:
        ax.set_xlabel("log2 expression")

    fig.suptitle("RMA Density", y=0.98)

    fig.tight_layout()
    fig.savefig(rma_density_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def quality_control_summary(path_configuration):
    master_rows = _read_master_rows(path_configuration)
    preprocessing_qc_rows = _read_preprocessing_qc_rows(path_configuration)
    sample_ids, expression_matrix, total_probes = _read_expression_matrix(
        path_configuration,
        max_probes=MAX_PLOT_PROBES
    )
    qc_summary = _build_qc_summary(master_rows, preprocessing_qc_rows, sample_ids, expression_matrix)

    _write_qc_summary(path_configuration, qc_summary)
    _write_qc_metrics(path_configuration, qc_summary, total_probes)
    _save_pca_plot(path_configuration, qc_summary)
    _save_boxplot(path_configuration, qc_summary, sample_ids, expression_matrix)
    _save_density_plot(path_configuration, qc_summary, sample_ids, expression_matrix)


if __name__ == "__main__":
    quality_control_summary(INTERNAL_CONFIGURATION)
