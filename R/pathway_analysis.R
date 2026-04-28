suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    paste(
      "Usage: pathway_analysis.R",
      "<limma_results_dir>",
      "<output_dir>",
      "<summary_csv>"
    )
  )
}

limma_results_dir <- args[[1]]
output_dir <- args[[2]]
summary_csv_path <- args[[3]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

MIN_PATHWAY_SIZE <- 15
MAX_PATHWAY_SIZE <- 500
FDR_CUTOFF <- 0.25

pathway_collection_specs <- list(
  hallmark = list(collection = "H", subcollections = list(NULL)),
  reactome = list(collection = "C2", subcollections = list("CP:REACTOME", "REACTOME")),
  kegg_medicus = list(collection = "C2", subcollections = list("CP:KEGG_MEDICUS", "KEGG_MEDICUS")),
  kegg_legacy = list(collection = "C2", subcollections = list("CP:KEGG_LEGACY", "CP:KEGG", "KEGG_LEGACY", "KEGG"))
)

empty_gsea_table <- function() {
  data.frame(
    contrast = character(),
    collection = character(),
    pathway = character(),
    pathway_description = character(),
    pval = numeric(),
    padj = numeric(),
    log2err = numeric(),
    ES = numeric(),
    NES = numeric(),
    size = integer(),
    leading_edge = character(),
    db_version = character(),
    stringsAsFactors = FALSE
  )
}

empty_ora_table <- function() {
  data.frame(
    contrast = character(),
    collection = character(),
    threshold = character(),
    pathway = character(),
    pathway_description = character(),
    overlap_size = integer(),
    significant_size = integer(),
    pathway_size = integer(),
    universe_size = integer(),
    gene_ratio = numeric(),
    background_ratio = numeric(),
    pval = numeric(),
    padj = numeric(),
    overlap_genes = character(),
    db_version = character(),
    stringsAsFactors = FALSE
  )
}

safe_column <- function(df, column_name, default = "") {
  if (column_name %in% names(df)) {
    return(df[[column_name]])
  }
  rep(default, nrow(df))
}

load_msigdb_table <- function(collection, subcollection = NULL) {
  formal_names <- names(formals(msigdbr::msigdbr))

  if ("collection" %in% formal_names) {
    msigdbr_args <- list(
      db_species = "HS",
      species = "Homo sapiens",
      collection = collection
    )
    if (!is.null(subcollection)) {
      msigdbr_args$subcollection <- subcollection
    }
  } else {
    msigdbr_args <- list(
      species = "Homo sapiens",
      category = collection
    )
    if (!is.null(subcollection)) {
      msigdbr_args$subcategory <- subcollection
    }
  }

  as.data.frame(do.call(msigdbr::msigdbr, msigdbr_args), stringsAsFactors = FALSE)
}

load_collection <- function(collection_name, spec) {
  for (subcollection in spec$subcollections) {
    subcollection_suffix <- if (is.null(subcollection)) "" else paste0(" / ", subcollection)
    subcollection_label <- if (is.null(subcollection)) "" else subcollection

    message(sprintf(
      "Loading MSigDB collection %s%s",
      spec$collection,
      subcollection_suffix
    ))

    msigdb_df <- tryCatch(
      load_msigdb_table(spec$collection, subcollection),
      error = function(error) {
        message(sprintf("Skipping unavailable collection candidate: %s", conditionMessage(error)))
        NULL
      }
    )

    if (!is.null(msigdb_df) && nrow(msigdb_df) > 0) {
      msigdb_df <- msigdb_df[!is.na(msigdb_df$gene_symbol) & nzchar(msigdb_df$gene_symbol), , drop = FALSE]
      if (nrow(msigdb_df) > 0) {
        return(list(
          name = collection_name,
          data = msigdb_df,
          subcollection = subcollection_label
        ))
      }
    }
  }

  warning(sprintf("No gene sets loaded for collection %s", collection_name))
  NULL
}

build_pathway_bundle <- function(collection_name, spec) {
  collection <- load_collection(collection_name, spec)
  if (is.null(collection)) {
    return(NULL)
  }

  msigdb_df <- collection$data
  pathways <- split(msigdb_df$gene_symbol, msigdb_df$gs_name)
  pathways <- lapply(pathways, unique)
  pathway_sizes <- lengths(pathways)
  pathways <- pathways[pathway_sizes >= MIN_PATHWAY_SIZE & pathway_sizes <= MAX_PATHWAY_SIZE]

  metadata_columns <- c("gs_name", "gs_description", "gs_collection", "gs_subcollection", "gs_exact_source", "gs_url")
  metadata_columns <- intersect(metadata_columns, names(msigdb_df))
  pathway_metadata <- unique(msigdb_df[, metadata_columns, drop = FALSE])

  db_version <- ""
  if ("db_version" %in% names(msigdb_df)) {
    db_version <- unique(msigdb_df$db_version)[[1]]
  }

  list(
    name = collection_name,
    pathways = pathways,
    metadata = pathway_metadata,
    db_version = db_version
  )
}

annotate_pathways <- function(result, pathway_bundle) {
  if (!nrow(result)) {
    return(result)
  }

  metadata <- pathway_bundle$metadata
  metadata_index <- match(result$pathway, metadata$gs_name)
  result$pathway_description <- safe_column(metadata, "gs_description")[metadata_index]
  result$pathway_description[is.na(result$pathway_description)] <- ""
  result$db_version <- pathway_bundle$db_version
  result
}

read_ranked_statistics <- function(rnk_path) {
  ranked <- read.table(
    rnk_path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(ranked) < 2) {
    stop(sprintf("Ranked file has fewer than two columns: %s", rnk_path))
  }

  ranked <- ranked[, 1:2, drop = FALSE]
  names(ranked) <- c("gene_symbol", "rank_score")
  ranked$gene_symbol <- trimws(ranked$gene_symbol)
  ranked$rank_score <- as.numeric(ranked$rank_score)
  ranked <- ranked[!is.na(ranked$rank_score) & !is.na(ranked$gene_symbol) & nzchar(ranked$gene_symbol), , drop = FALSE]

  ranked$abs_rank_score <- abs(ranked$rank_score)
  ranked <- ranked[order(-ranked$abs_rank_score, ranked$gene_symbol), , drop = FALSE]
  ranked <- ranked[!duplicated(ranked$gene_symbol), , drop = FALSE]

  stats <- ranked$rank_score
  names(stats) <- ranked$gene_symbol
  sort(stats, decreasing = TRUE)
}

run_gsea <- function(contrast, stats, pathway_bundle) {
  if (!length(stats) || !length(pathway_bundle$pathways)) {
    return(empty_gsea_table())
  }

  set.seed(1)
  result <- as.data.frame(fgsea::fgsea(
    pathways = pathway_bundle$pathways,
    stats = stats,
    minSize = MIN_PATHWAY_SIZE,
    maxSize = MAX_PATHWAY_SIZE,
    eps = 0
  ), stringsAsFactors = FALSE)

  if (!nrow(result)) {
    return(empty_gsea_table())
  }

  result$leading_edge <- vapply(result$leadingEdge, paste, character(1), collapse = ";")
  result$leadingEdge <- NULL
  result$contrast <- contrast
  result$collection <- pathway_bundle$name
  result <- annotate_pathways(result, pathway_bundle)
  result <- result[order(result$padj, result$pval), , drop = FALSE]

  result[, c(
    "contrast",
    "collection",
    "pathway",
    "pathway_description",
    "pval",
    "padj",
    "log2err",
    "ES",
    "NES",
    "size",
    "leading_edge",
    "db_version"
  )]
}

read_significant_genes <- function(pathway_input_path) {
  if (!file.exists(pathway_input_path)) {
    return(character())
  }

  table <- read.csv(pathway_input_path, stringsAsFactors = FALSE)
  if (!("gene_symbol" %in% names(table))) {
    return(character())
  }

  genes <- unique(trimws(table$gene_symbol))
  genes[!is.na(genes) & nzchar(genes)]
}

run_ora <- function(contrast, stats, significant_genes, threshold, pathway_bundle) {
  universe <- unique(names(stats))
  significant_genes <- intersect(unique(significant_genes), universe)

  if (!length(universe) || !length(significant_genes) || !length(pathway_bundle$pathways)) {
    return(empty_ora_table())
  }

  rows <- lapply(names(pathway_bundle$pathways), function(pathway_name) {
    pathway_genes <- intersect(unique(pathway_bundle$pathways[[pathway_name]]), universe)
    pathway_size <- length(pathway_genes)

    if (pathway_size < MIN_PATHWAY_SIZE || pathway_size > MAX_PATHWAY_SIZE) {
      return(NULL)
    }

    overlap_genes <- intersect(significant_genes, pathway_genes)
    overlap_size <- length(overlap_genes)

    if (overlap_size == 0) {
      return(NULL)
    }

    contingency <- matrix(
      c(
        overlap_size,
        length(significant_genes) - overlap_size,
        pathway_size - overlap_size,
        length(universe) - length(significant_genes) - pathway_size + overlap_size
      ),
      nrow = 2
    )

    data.frame(
      contrast = contrast,
      collection = pathway_bundle$name,
      threshold = threshold,
      pathway = pathway_name,
      overlap_size = overlap_size,
      significant_size = length(significant_genes),
      pathway_size = pathway_size,
      universe_size = length(universe),
      gene_ratio = overlap_size / length(significant_genes),
      background_ratio = pathway_size / length(universe),
      pval = fisher.test(contingency, alternative = "greater")$p.value,
      overlap_genes = paste(sort(overlap_genes), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(result) || !nrow(result)) {
    return(empty_ora_table())
  }

  result$padj <- p.adjust(result$pval, method = "BH")
  result <- annotate_pathways(result, pathway_bundle)
  result <- result[order(result$padj, result$pval), , drop = FALSE]

  result[, c(
    "contrast",
    "collection",
    "threshold",
    "pathway",
    "pathway_description",
    "overlap_size",
    "significant_size",
    "pathway_size",
    "universe_size",
    "gene_ratio",
    "background_ratio",
    "pval",
    "padj",
    "overlap_genes",
    "db_version"
  )]
}

write_result <- function(result, output_path) {
  write.csv(result, output_path, row.names = FALSE)
}

summarize_result <- function(result, analysis, contrast, collection, threshold, n_input_genes, n_gene_sets) {
  data.frame(
    analysis = analysis,
    contrast = contrast,
    collection = collection,
    threshold = threshold,
    n_input_genes = n_input_genes,
    n_gene_sets = n_gene_sets,
    n_results = nrow(result),
    n_fdr_lt_0_05 = sum(result$padj < 0.05, na.rm = TRUE),
    n_fdr_lt_0_25 = sum(result$padj < FDR_CUTOFF, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

message("Loading MSigDB pathway collections")
pathway_bundles <- lapply(names(pathway_collection_specs), function(collection_name) {
  build_pathway_bundle(collection_name, pathway_collection_specs[[collection_name]])
})
names(pathway_bundles) <- names(pathway_collection_specs)
pathway_bundles <- pathway_bundles[!vapply(pathway_bundles, is.null, logical(1))]

if (!length(pathway_bundles)) {
  stop("No MSigDB pathway collections could be loaded.")
}

ranked_files <- list.files(limma_results_dir, pattern = "^de_.*_ranked\\.rnk$", full.names = TRUE)
if (!length(ranked_files)) {
  stop(sprintf("No ranked limma files found in %s", limma_results_dir))
}

gsea_results <- list()
ora_results <- list()
summary_rows <- list()

for (ranked_file in ranked_files) {
  contrast <- sub("_ranked\\.rnk$", "", basename(ranked_file))
  message(sprintf("Running pathway analysis for %s", contrast))
  stats <- read_ranked_statistics(ranked_file)

  p05_genes <- read_significant_genes(file.path(limma_results_dir, paste0(contrast, "_pathway_input_p05.csv")))
  fdr015_genes <- read_significant_genes(file.path(limma_results_dir, paste0(contrast, "_pathway_input_fdr015.csv")))

  for (pathway_bundle in pathway_bundles) {
    gsea <- run_gsea(contrast, stats, pathway_bundle)
    gsea_output_path <- file.path(output_dir, paste0("gsea_", contrast, "_", pathway_bundle$name, ".csv"))
    write_result(gsea, gsea_output_path)
    gsea_results[[paste(contrast, pathway_bundle$name, sep = ":")]] <- gsea
    summary_rows[[paste("gsea", contrast, pathway_bundle$name, sep = ":")]] <- summarize_result(
      gsea,
      "gsea",
      contrast,
      pathway_bundle$name,
      "ranked_all_genes",
      length(stats),
      length(pathway_bundle$pathways)
    )

    ora_inputs <- list(p05 = p05_genes, fdr015 = fdr015_genes)
    for (threshold in names(ora_inputs)) {
      ora <- run_ora(contrast, stats, ora_inputs[[threshold]], threshold, pathway_bundle)
      ora_output_path <- file.path(output_dir, paste0("ora_", contrast, "_", pathway_bundle$name, "_", threshold, ".csv"))
      write_result(ora, ora_output_path)
      ora_results[[paste(contrast, pathway_bundle$name, threshold, sep = ":")]] <- ora
      summary_rows[[paste("ora", contrast, pathway_bundle$name, threshold, sep = ":")]] <- summarize_result(
        ora,
        "ora",
        contrast,
        pathway_bundle$name,
        threshold,
        length(ora_inputs[[threshold]]),
        length(pathway_bundle$pathways)
      )
    }
  }
}

combined_gsea <- if (length(gsea_results)) do.call(rbind, gsea_results) else empty_gsea_table()
combined_ora <- if (length(ora_results)) do.call(rbind, ora_results) else empty_ora_table()
summary_table <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()

write_result(combined_gsea, file.path(output_dir, "gsea_all_results.csv"))
write_result(combined_ora, file.path(output_dir, "ora_all_results.csv"))
write_result(summary_table, summary_csv_path)

message("Pathway-analysis workflow complete")



message("Generating pathway-analysis interpretation table")

pathway_dir <- file.path("artifacts", "pathway")

gsea <- read.csv(
  file.path(pathway_dir, "gsea_all_results.csv"),
  stringsAsFactors = FALSE
)

ora <- read.csv(
  file.path(pathway_dir, "ora_all_results.csv"),
  stringsAsFactors = FALSE
)

gsea$direction <- ifelse(
  gsea$NES > 0,
  "higher in Topiramate-vs-Placebo delta",
  "lower in Topiramate-vs-Placebo delta"
)

gsea$evidence_level <- ifelse(
  gsea$padj < 0.05,
  "strong",
  ifelse(gsea$padj < 0.25, "exploratory", "not_significant")
)

gsea_hits <- gsea[gsea$padj < 0.25, ]

gsea_hits <- gsea_hits[order(
  gsea_hits$contrast,
  gsea_hits$padj,
  -abs(gsea_hits$NES)
), ]

gsea_hits <- gsea_hits[, c(
  "contrast",
  "collection",
  "pathway",
  "pathway_description",
  "NES",
  "pval",
  "padj",
  "direction",
  "evidence_level",
  "size",
  "leading_edge"
)]

write.csv(
  gsea_hits,
  file.path(pathway_dir, "top_gsea_interpretation_candidates.csv"),
  row.names = FALSE
)

ora_hits <- ora[ora$padj < 0.05, ]

ora_hits <- ora_hits[order(
  ora_hits$contrast,
  ora_hits$threshold,
  ora_hits$padj
), ]

write.csv(
  ora_hits,
  file.path(pathway_dir, "top_ora_supporting_candidates.csv"),
  row.names = FALSE
)


message("Pathway-analysis interpretation table complete")