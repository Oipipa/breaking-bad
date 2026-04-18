suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    paste(
      "Usage: limma_secondary_outcomes.R",
      "<master_sheet.csv>",
      "<expression_bundle.rds>",
      "<output_dir>",
      "<comparison_summary.csv>"
    )
  )
}

master_sheet_path <- args[[1]]
expression_bundle_path <- args[[2]]
output_dir <- args[[3]]
comparison_summary_path <- args[[4]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

annotation_table <- function(eset, probe_ids) {
  feature_data <- as.data.frame(fData(eset), stringsAsFactors = FALSE)
  feature_data$probe_set_id <- rownames(feature_data)
  if (!("gene_symbol" %in% names(feature_data))) {
    feature_data$gene_symbol <- ""
  }
  if (!("entrez_id" %in% names(feature_data))) {
    feature_data$entrez_id <- ""
  }
  if (!("gene_name" %in% names(feature_data))) {
    feature_data$gene_name <- ""
  }

  merged <- merge(
    data.frame(probe_set_id = probe_ids, stringsAsFactors = FALSE),
    feature_data[, c("probe_set_id", "gene_symbol", "entrez_id", "gene_name")],
    by = "probe_set_id",
    all.x = TRUE,
    sort = FALSE
  )
  merged[is.na(merged)] <- ""
  merged
}

first_sample_by_subject <- function(df, label) {
  if (!nrow(df)) {
    return(data.frame())
  }

  split_rows <- split(df, df$subject_id)
  selected <- lapply(split_rows, function(group_df) {
    if (nrow(group_df) > 1) {
      message(sprintf(
        "Subject %s has %d %s samples; using first sample %s",
        group_df$subject_id[[1]],
        nrow(group_df),
        label,
        group_df$geo_sample_id[[1]]
      ))
    }
    group_df[1, , drop = FALSE]
  })
  do.call(rbind, selected)
}

build_delta_dataset <- function(expression_matrix, metadata, target_timepoint) {
  baseline_samples <- metadata[metadata$timepoint == "baseline", , drop = FALSE]
  target_samples <- metadata[metadata$timepoint == target_timepoint, , drop = FALSE]

  baseline_samples <- first_sample_by_subject(baseline_samples, "baseline")
  target_samples <- first_sample_by_subject(target_samples, target_timepoint)

  baseline_by_subject <- split(baseline_samples, baseline_samples$subject_id)
  target_by_subject <- split(target_samples, target_samples$subject_id)
  subjects <- intersect(names(baseline_by_subject), names(target_by_subject))

  delta_columns <- list()
  delta_metadata <- list()
  skipped_subjects <- 0

  for (subject_id in subjects) {
    baseline_row <- baseline_by_subject[[subject_id]][1, , drop = FALSE]
    target_row <- target_by_subject[[subject_id]][1, , drop = FALSE]

    if (baseline_row$treatment_arm != target_row$treatment_arm) {
      skipped_subjects <- skipped_subjects + 1
      next
    }

    baseline_id <- baseline_row$geo_sample_id
    target_id <- target_row$geo_sample_id

    if (!(baseline_id %in% colnames(expression_matrix)) || !(target_id %in% colnames(expression_matrix))) {
      skipped_subjects <- skipped_subjects + 1
      next
    }

    delta_vector <- expression_matrix[, target_id] - expression_matrix[, baseline_id]
    delta_name <- paste(subject_id, target_timepoint, sep = "_")
    delta_columns[[delta_name]] <- delta_vector

    delta_metadata[[delta_name]] <- data.frame(
      delta_sample_id = delta_name,
      subject_id = subject_id,
      treatment_arm = target_row$treatment_arm,
      timepoint = target_timepoint,
      baseline_sample_id = baseline_id,
      target_sample_id = target_id,
      stringsAsFactors = FALSE
    )
  }

  if (!length(delta_columns)) {
    return(list(matrix = NULL, metadata = data.frame(), skipped_subjects = skipped_subjects))
  }

  delta_matrix <- do.call(cbind, delta_columns)
  rownames(delta_matrix) <- rownames(expression_matrix)

  delta_metadata_df <- do.call(rbind, delta_metadata)
  rownames(delta_metadata_df) <- delta_metadata_df$delta_sample_id

  list(matrix = delta_matrix, metadata = delta_metadata_df, skipped_subjects = skipped_subjects)
}

prepare_table <- function(eset, fit, coef_name) {
  table <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  table$probe_set_id <- rownames(table)
  annotations <- annotation_table(eset, table$probe_set_id)
  table <- merge(table, annotations, by = "probe_set_id", all.x = TRUE, sort = FALSE)

  ordered_columns <- c(
    "probe_set_id",
    "gene_symbol",
    "entrez_id",
    "gene_name",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
  table[, ordered_columns]
}

write_pathway_inputs <- function(table, output_dir, prefix) {
  ranked <- table[!is.na(table$gene_symbol) & nzchar(table$gene_symbol), c("gene_symbol", "logFC", "P.Value", "adj.P.Val")]
  if (!nrow(ranked)) {
    write.csv(ranked, file.path(output_dir, paste0(prefix, "_pathway_input_all.csv")), row.names = FALSE)
    write.table(ranked[, c("gene_symbol"), drop = FALSE], file.path(output_dir, paste0(prefix, "_ranked.rnk")), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    return(invisible(NULL))
  }

  ranked <- ranked[order(ranked$P.Value, -abs(ranked$logFC)), ]
  ranked <- ranked[!duplicated(ranked$gene_symbol), ]
  write.csv(ranked, file.path(output_dir, paste0(prefix, "_pathway_input_all.csv")), row.names = FALSE)

  sig_nominal <- ranked[ranked$P.Value < 0.05, , drop = FALSE]
  sig_fdr <- ranked[ranked$adj.P.Val < 0.15, , drop = FALSE]
  write.csv(sig_nominal, file.path(output_dir, paste0(prefix, "_pathway_input_p05.csv")), row.names = FALSE)
  write.csv(sig_fdr, file.path(output_dir, paste0(prefix, "_pathway_input_fdr015.csv")), row.names = FALSE)

  rank_score <- sign(ranked$logFC) * (-log10(pmax(ranked$P.Value, 1e-300)))
  ranked_rnk <- data.frame(gene_symbol = ranked$gene_symbol, rank_score = rank_score, stringsAsFactors = FALSE)
  write.table(ranked_rnk, file.path(output_dir, paste0(prefix, "_ranked.rnk")), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

run_public_comparison <- function(eset, delta_matrix, delta_metadata, timepoint, output_dir) {
  comparison_name <- paste("public", tolower(timepoint), "topiramate_vs_placebo", sep = "_")
  comparison_rows <- delta_metadata[delta_metadata$timepoint == timepoint, , drop = FALSE]

  n_topiramate <- sum(comparison_rows$treatment_arm == "Topiramate")
  n_placebo <- sum(comparison_rows$treatment_arm == "Placebo")

  if (nrow(comparison_rows) == 0 || n_topiramate < 2 || n_placebo < 2) {
    return(list(
      success = FALSE,
      reason = "insufficient_group_size",
      summary = data.frame(
        mode = "public",
        comparison = comparison_name,
        n_subjects = nrow(comparison_rows),
        n_topiramate = n_topiramate,
        n_placebo = n_placebo,
        n_nominal_p_lt_0_05 = 0,
        n_fdr_lt_0_15 = 0,
        n_fdr_lt_0_05 = 0,
        stringsAsFactors = FALSE
      )
    ))
  }

  sample_ids <- comparison_rows$delta_sample_id
  comparison_matrix <- delta_matrix[, sample_ids, drop = FALSE]
  arm <- factor(comparison_rows$treatment_arm, levels = c("Placebo", "Topiramate"))
  design <- model.matrix(~ arm)

  fit <- lmFit(comparison_matrix, design)
  fit <- eBayes(fit)
  table <- prepare_table(eset, fit, "armTopiramate")

  output_path <- file.path(output_dir, paste0("de_", comparison_name, "_limma.csv"))
  write.csv(table, output_path, row.names = FALSE)
  write_pathway_inputs(table, output_dir, paste0("de_", comparison_name))

  summary <- data.frame(
    mode = "public",
    comparison = comparison_name,
    n_subjects = nrow(comparison_rows),
    n_topiramate = n_topiramate,
    n_placebo = n_placebo,
    n_nominal_p_lt_0_05 = sum(table$P.Value < 0.05, na.rm = TRUE),
    n_fdr_lt_0_15 = sum(table$adj.P.Val < 0.15, na.rm = TRUE),
    n_fdr_lt_0_05 = sum(table$adj.P.Val < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(success = TRUE, reason = "ok", table = table, summary = summary)
}

message("Reading master sample sheet")
master_sheet <- read.csv(master_sheet_path, stringsAsFactors = FALSE)
required_master_columns <- c("geo_sample_id", "subject_id", "treatment_arm", "timepoint")
missing_master_columns <- setdiff(required_master_columns, names(master_sheet))
if (length(missing_master_columns)) {
  stop(sprintf("Missing required columns in master sample sheet: %s", paste(missing_master_columns, collapse = ", ")))
}

message("Reading expression bundle")
eset <- readRDS(expression_bundle_path)
expression_matrix <- exprs(eset)

sample_order <- colnames(expression_matrix)
metadata <- master_sheet[match(sample_order, master_sheet$geo_sample_id), , drop = FALSE]
if (any(is.na(metadata$geo_sample_id))) {
  stop("Expression bundle sample IDs could not be fully mapped to master sample sheet")
}

metadata$subject_id <- as.character(metadata$subject_id)
if (!all(metadata$timepoint %in% c("baseline", "week8", "week12"))) {
  unexpected <- unique(metadata$timepoint[!(metadata$timepoint %in% c("baseline", "week8", "week12"))])
  stop(sprintf("Unexpected timepoint values found: %s", paste(unexpected, collapse = ", ")))
}

message("Building baseline-normalized delta matrices")
week8_delta <- build_delta_dataset(expression_matrix, metadata, "week8")
week12_delta <- build_delta_dataset(expression_matrix, metadata, "week12")

if (is.null(week8_delta$matrix) && is.null(week12_delta$matrix)) {
  stop("No analyzable delta samples were created. Check baseline and week labels.")
}

combined_matrix <- do.call(cbind, Filter(Negate(is.null), list(week8_delta$matrix, week12_delta$matrix)))
combined_metadata <- rbind(week8_delta$metadata, week12_delta$metadata)

results <- list()
summary_rows <- list()

message("Running public no-label limma comparisons")
public_comparisons <- list("week8", "week12")
for (tp in public_comparisons) {
  result <- run_public_comparison(eset, combined_matrix, combined_metadata, tp, output_dir)
  results[[tp]] <- result
  summary_rows[[tp]] <- cbind(result$summary, status = ifelse(result$success, "ok", result$reason), stringsAsFactors = FALSE)
}

week8_table <- if (!is.null(results[["week8"]]$table)) results[["week8"]]$table else data.frame()
week12_table <- if (!is.null(results[["week12"]]$table)) results[["week12"]]$table else data.frame()

if (nrow(week8_table) && nrow(week12_table)) {
  week8_sel <- week8_table[, c("probe_set_id", "gene_symbol", "entrez_id", "gene_name", "logFC", "P.Value", "adj.P.Val")]
  names(week8_sel) <- c("probe_set_id", "gene_symbol", "entrez_id", "gene_name", "logFC_week8", "P_week8", "FDR_week8")
  week12_sel <- week12_table[, c("probe_set_id", "logFC", "P.Value", "adj.P.Val")]
  names(week12_sel) <- c("probe_set_id", "logFC_week12", "P_week12", "FDR_week12")

  overlap <- merge(week8_sel, week12_sel, by = "probe_set_id", all = FALSE)
  overlap$direction_consistent <- sign(overlap$logFC_week8) == sign(overlap$logFC_week12)
  overlap$significant_nominal_both <- overlap$P_week8 < 0.05 & overlap$P_week12 < 0.05
  overlap$significant_fdr_both <- overlap$FDR_week8 < 0.15 & overlap$FDR_week12 < 0.15
  write.csv(overlap, file.path(output_dir, "overlap_public_week8_week12.csv"), row.names = FALSE)

  overlap_pathway <- overlap[overlap$significant_nominal_both & nzchar(overlap$gene_symbol), c("gene_symbol", "logFC_week8", "P_week8", "FDR_week8", "logFC_week12", "P_week12", "FDR_week12", "direction_consistent")]
  overlap_pathway <- overlap_pathway[order(overlap_pathway$P_week8 + overlap_pathway$P_week12), ]
  overlap_pathway <- overlap_pathway[!duplicated(overlap_pathway$gene_symbol), ]
  write.csv(overlap_pathway, file.path(output_dir, "overlap_public_week8_week12_pathway_input.csv"), row.names = FALSE)
}

metadata_summary <- data.frame(
  mode = "public",
  metric = c("n_delta_week8", "n_delta_week12", "n_unique_subjects_week8", "n_unique_subjects_week12"),
  value = c(
    nrow(week8_delta$metadata),
    nrow(week12_delta$metadata),
    length(unique(week8_delta$metadata$subject_id)),
    length(unique(week12_delta$metadata$subject_id))
  ),
  stringsAsFactors = FALSE
)

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, comparison_summary_path, row.names = FALSE)
write.csv(metadata_summary, file.path(output_dir, "limma_metadata_summary.csv"), row.names = FALSE)

message("limma workflow complete")
