suppressPackageStartupMessages({
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Missing package 'fgsea'. Install it first (BiocManager::install('fgsea')).")
  }
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Missing package 'msigdbr'. Install it first (install.packages('msigdbr')).")
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

collapse_to_gene_ranks <- function(df) {
  required_cols <- c("gene_symbol", "t")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  df <- df[!is.na(df$gene_symbol) & nzchar(df$gene_symbol), , drop = FALSE]
  if (!nrow(df)) {
    return(numeric(0))
  }

  df <- df[order(abs(df$t), decreasing = TRUE), , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol), , drop = FALSE]

  ranks <- df$t
  names(ranks) <- df$gene_symbol
  sort(ranks, decreasing = TRUE)
}

build_metabolic_pathways <- function(species_name = "Homo sapiens") {
  msig <- msigdbr::msigdbr(species = species_name, category = "C2")
  if (!nrow(msig)) {
    stop("No pathways returned from msigdbr.")
  }

  metabolic_mask <- grepl(
    "METABOL|GLYCOL|LIPID|FATTY|OXIDATIVE|TCA|CITRATE|AMINO",
    msig$gs_name,
    ignore.case = TRUE
  )

  msig <- msig[metabolic_mask, c("gs_name", "gene_symbol")]
  msig <- msig[!is.na(msig$gene_symbol) & nzchar(msig$gene_symbol), , drop = FALSE]

  pathways <- split(msig$gene_symbol, msig$gs_name)
  pathways <- lapply(pathways, unique)

  pathways[sapply(pathways, length) >= 10]
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

pathways <- build_metabolic_pathways("Homo sapiens")
if (!length(pathways)) {
  stop("No metabolic pathways found after filtering.")
}

out_dir <- file.path(run_dir, "pathways")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1)
summary_rows <- list()

for (csv_path in all_files) {
  contrast_name <- sub("_all\\.csv$", "", basename(csv_path))
  message(sprintf("Running ranked enrichment: %s", contrast_name))

  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  ranks <- collapse_to_gene_ranks(df)

  if (length(ranks) < 50) {
    message(sprintf("Skipping %s: too few ranked genes (%d).", contrast_name, length(ranks)))
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      contrast = contrast_name,
      n_ranked_genes = length(ranks),
      n_sig_pathways = 0,
      note = "skipped_too_few_genes",
      stringsAsFactors = FALSE
    )
    next
  }

  fg <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 1e-10
  )

  fg <- as.data.frame(fg)
  if (nrow(fg)) {
    fg <- fg[order(fg$padj, -abs(fg$NES)), , drop = FALSE]
    fg$leadingEdge <- vapply(
      fg$leadingEdge,
      function(x) paste(x, collapse = ";"),
      FUN.VALUE = character(1)
    )
  }

  sig <- fg[!is.na(fg$padj) & fg$padj < 0.05, , drop = FALSE]

  write.csv(
    fg,
    file = file.path(out_dir, paste0(contrast_name, "_metabolic_gsea_all.csv")),
    row.names = FALSE
  )
  write.csv(
    sig,
    file = file.path(out_dir, paste0(contrast_name, "_metabolic_gsea_sig.csv")),
    row.names = FALSE
  )

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    contrast = contrast_name,
    n_ranked_genes = length(ranks),
    n_sig_pathways = nrow(sig),
    note = "ok",
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file = file.path(out_dir, "metabolic_gsea_summary.csv"), row.names = FALSE)

message(sprintf("Done. Ranked metabolic enrichment written to: %s", out_dir))
