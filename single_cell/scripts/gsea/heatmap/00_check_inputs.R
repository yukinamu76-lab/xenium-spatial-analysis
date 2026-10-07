suppressPackageStartupMessages({library(readr); library(dplyr); library(purrr); library(tibble)})
source("config.R")
required_cols <- c("pathway", "padj", "NES")
checks <- purrr::imap_dfr(input_files, function(path, cluster) {
  if (!file.exists(path)) {
    return(tibble(cluster = cluster, file = path, exists = FALSE, n_rows = NA_integer_, missing_columns = "file missing"))
  }
  x <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  miss <- setdiff(required_cols, names(x))
  tibble(cluster = cluster, file = path, exists = TRUE, n_rows = nrow(x), missing_columns = paste(miss, collapse = "; "))
})
print(checks, n = Inf)
if (any(!checks$exists) || any(nzchar(checks$missing_columns))) stop("Input check failed. See table above.")
if (!file.exists(classification_file)) stop("Classification file not found: ", classification_file)
message("All input files and the classification file passed basic checks.")
