args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: _deps.R <install_missing:true|false> <package> [<package> ...]")
}

install_missing <- tolower(args[[1]]) == "true"
packages <- args[-1]

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  if (!install_missing) {
    stop("BiocManager is required but not installed.")
  }
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  if (!install_missing) {
    stop(sprintf("Missing required R/Bioconductor packages: %s", paste(missing, collapse = ", ")))
  }
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
