# Session information for the scRNA-seq analysis environment
# Run this script in the same R environment that was used for scripts/scRNA/.
# It does not load any data files or run the analysis.

packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "harmony",
  "patchwork",
  "readr",
  "RColorBrewer",
  "slingshot",
  "speckle",
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

# Attach the packages used directly by the public scRNA scripts.
for (pkg in packages) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

out_file <- "sessionInfo_scRNA.txt"

sink(out_file)
cat("scRNA-seq analysis environment\n")
cat("Generated:", format(Sys.time(), tz = "Asia/Tokyo"), "\n\n")

cat("Package versions explicitly used by scripts/scRNA/\n")
for (pkg in packages) {
  cat(sprintf("%-20s %s\n", pkg, as.character(packageVersion(pkg))))
}
cat("\n")

sessionInfo()
sink()

message("Wrote: ", normalizePath(out_file, mustWork = FALSE))
