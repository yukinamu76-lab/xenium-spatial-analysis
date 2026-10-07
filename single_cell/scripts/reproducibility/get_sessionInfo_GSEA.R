# Session information for the GSEA analysis environment
# Run this script in the same R environment that was used for scripts/gsea/.
# It does not load any input tables or run GSEA.

packages <- c(
  "dplyr",
  "fgsea",
  "ggplot2",
  "igraph",
  "msigdbr",
  "pheatmap",
  "purrr",
  "readr",
  "scales",
  "stringr",
  "tibble",
  "tidyr"
)

missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "The following packages are not installed in this environment: ",
    paste(missing_packages, collapse = ", ")
  )
}

# Attach the packages used directly by the public GSEA scripts.
for (pkg in packages) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

out_file <- "sessionInfo_GSEA.txt"

sink(out_file)
cat("GSEA analysis environment\n")
cat("Generated:", format(Sys.time(), tz = "Asia/Tokyo"), "\n\n")

cat("Package versions explicitly used by scripts/gsea/\n")
for (pkg in packages) {
  cat(sprintf("%-20s %s\n", pkg, as.character(packageVersion(pkg))))
}
cat("\n")

sessionInfo()
sink()

message("Wrote: ", normalizePath(out_file, mustWork = FALSE))
