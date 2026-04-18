suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Missing package 'ggplot2'. Install it first, e.g. install.packages('ggplot2').")
  }
})

find_latest_run_dir <- function(base_runs_dir) {
  candidates <- list.dirs(base_runs_dir, recursive = FALSE, full.names = TRUE)
  if (!length(candidates)) {
    stop(sprintf("No run directories found under: %s", base_runs_dir))
  }
  info <- file.info(candidates)
  candidates[which.max(info$mtime)]
}

args <- commandArgs(trailingOnly = TRUE)
base_runs_dir <- file.path("artifacts", "deg_runs")
run_dir <- if (length(args) >= 1) args[[1]] else find_latest_run_dir(base_runs_dir)

if (!dir.exists(run_dir)) {
  stop(sprintf("Run directory not found: %s", run_dir))
}

all_files <- list.files(run_dir, pattern = "_all\\.csv$", full.names = TRUE)
if (!length(all_files)) {
  stop(sprintf("No *_all.csv files found in: %s", run_dir))
}

out_dir <- file.path(run_dir, "volcano")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fc_cutoff <- log2(1.2)
p_cutoff <- 0.05

summary_rows <- list()

for (csv_path in all_files) {
  contrast_name <- sub("_all\\.csv$", "", basename(csv_path))
  message(sprintf("Building volcano plot: %s", contrast_name))

  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  required_cols <- c("logFC", "adj.P.Val")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required columns in %s: %s",
      csv_path,
      paste(missing_cols, collapse = ", ")
    ))
  }

  df$neg_log10_adj_p <- -log10(pmax(df$adj.P.Val, .Machine$double.xmin))
  df$status <- "NS"
  df$status[df$adj.P.Val < p_cutoff & df$logFC >= fc_cutoff] <- "Up"
  df$status[df$adj.P.Val < p_cutoff & df$logFC <= -fc_cutoff] <- "Down"

  n_total <- nrow(df)
  n_up <- sum(df$status == "Up", na.rm = TRUE)
  n_down <- sum(df$status == "Down", na.rm = TRUE)
  n_sig <- n_up + n_down

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = logFC, y = neg_log10_adj_p, color = status)
  ) +
    ggplot2::geom_point(alpha = 0.7, size = 1.0, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    ggplot2::scale_color_manual(values = c(NS = "grey70", Up = "firebrick3", Down = "steelblue3")) +
    ggplot2::labs(
      title = sprintf("Volcano: %s", contrast_name),
      subtitle = sprintf("Significant probes: %d (Up: %d, Down: %d)", n_sig, n_up, n_down),
      x = "log2 fold-change",
      y = "-log10 adjusted p-value",
      color = "Class"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  png_out <- file.path(out_dir, paste0(contrast_name, "_volcano.png"))
  pdf_out <- file.path(out_dir, paste0(contrast_name, "_volcano.pdf"))

  ggplot2::ggsave(filename = png_out, plot = p, width = 8, height = 6, dpi = 150)
  ggplot2::ggsave(filename = pdf_out, plot = p, width = 8, height = 6)

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    contrast = contrast_name,
    n_total = n_total,
    n_sig = n_sig,
    n_up = n_up,
    n_down = n_down,
    fc_cutoff = fc_cutoff,
    adj_p_cutoff = p_cutoff,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file = file.path(out_dir, "volcano_summary.csv"), row.names = FALSE)

message(sprintf("Done. Volcano plots written to: %s", out_dir))
