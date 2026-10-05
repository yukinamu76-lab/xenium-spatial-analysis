# ============================================================
# Q-value calculation
#
# This script calculates q-values from nominal P values using
# the qvalue R package.
#
# Input:
#   A CSV file containing a column named "p".
#   Each row corresponds to one tested feature.
#
# Output:
#   The original input table with appended q-value and
#   local false discovery rate (lfdr) columns.
# ============================================================

library(qvalue)

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

input_file <- "input_pvalues.csv"
output_file <- "qvalue_results.csv"

# ------------------------------------------------------------
# Read input
# ------------------------------------------------------------

data <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE
)

if (!"p" %in% colnames(data)) {
  stop("Input CSV must contain a column named 'p'.")
}

pvalues <- data$p

# ------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------

if (any(is.na(pvalues))) {
  stop("The p-value column contains NA values.")
}

if (any(pvalues < 0 | pvalues > 1)) {
  stop("All p-values must be between 0 and 1.")
}

# ------------------------------------------------------------
# Q-value calculation
# ------------------------------------------------------------

qobj <- qvalue(
  p = pvalues
)

data$qvalue <- qobj$qvalues
data$lfdr <- qobj$lfdr

pi0 <- qobj$pi0

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("Number of tested features:", length(pvalues), "\n")
cat("Estimated pi0:", pi0, "\n")

print(summary(qobj))

# ------------------------------------------------------------
# Optional diagnostic plots
# ------------------------------------------------------------

hist(qobj)
plot(qobj)

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------

write.csv(
  data,
  file = output_file,
  row.names = FALSE,
  quote = TRUE
)

# ------------------------------------------------------------
# Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo_qvalue.txt"
)
