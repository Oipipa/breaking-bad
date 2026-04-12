suppressPackageStartupMessages({
  library(affy)
  library(AnnotationDbi)
  library(Biobase)
  library(hgu133plus2.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    paste(
      "Usage: quality_control.R",
      "<master_sheet.csv>",
      "<raw_dir>",
      "<expression_matrix.csv>",
      "<preprocessing_qc.csv>",
      "<expression_bundle.rds>"
    )
  )
}

master_sheet_path <- args[[1]]
raw_dir <- args[[2]]
expression_matrix_path <- args[[3]]
preprocessing_qc_path <- args[[4]]
expression_bundle_path <- args[[5]]

collapse_annotation <- function(values) {
  values <- unique(trimws(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) {
    return("")
  }
  paste(values, collapse = ";")
}

message("Reading master sample sheet")
metadata <- read.csv(master_sheet_path, stringsAsFactors = FALSE)
raw_files <- file.path(raw_dir, metadata$raw_filename)

message(sprintf("Reading %d CEL files", length(raw_files)))
raw_data <- ReadAffy(filenames = raw_files)
sampleNames(raw_data) <- metadata$geo_sample_id

message("Computing MAS5 present calls")
present_calls <- mas5calls(raw_data)
present_pct <- colMeans(exprs(present_calls) == "P") * 100

message("Running RMA normalization")
eset <- rma(raw_data)
sampleNames(eset) <- metadata$geo_sample_id
expression_matrix <- exprs(eset)

message("Looking up probe annotations")
probe_ids <- rownames(expression_matrix)
annotation_raw <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "PROBEID"
)
annotation_summary <- aggregate(
  annotation_raw[c("SYMBOL", "ENTREZID", "GENENAME")],
  list(probe_set_id = annotation_raw$PROBEID),
  collapse_annotation
)

annotation_df <- data.frame(
  probe_set_id = probe_ids,
  gene_symbol = annotation_summary$SYMBOL[match(probe_ids, annotation_summary$probe_set_id)],
  entrez_id = annotation_summary$ENTREZID[match(probe_ids, annotation_summary$probe_set_id)],
  gene_name = annotation_summary$GENENAME[match(probe_ids, annotation_summary$probe_set_id)],
  stringsAsFactors = FALSE,
  row.names = probe_ids
)
annotation_df[is.na(annotation_df)] <- ""

preprocessing_qc <- data.frame(
  geo_sample_id = metadata$geo_sample_id,
  present_pct = present_pct,
  stringsAsFactors = FALSE
)

message("Writing preprocessing artifacts")
expression_df <- cbind(
  annotation_df,
  as.data.frame(expression_matrix, check.names = FALSE)
)
write.csv(expression_df, expression_matrix_path, row.names = FALSE)
write.csv(preprocessing_qc, preprocessing_qc_path, row.names = FALSE)

pheno_data <- data.frame(
  metadata,
  present_pct = present_pct,
  stringsAsFactors = FALSE
)
row.names(pheno_data) <- pheno_data$geo_sample_id

feature_metadata <- annotation_df[, c("gene_symbol", "entrez_id", "gene_name"), drop = FALSE]
row.names(feature_metadata) <- annotation_df$probe_set_id
phenoData(eset) <- new("AnnotatedDataFrame", data = pheno_data)
featureData(eset) <- new("AnnotatedDataFrame", data = feature_metadata)
saveRDS(eset, expression_bundle_path)

message("Quality-control workflow complete")
