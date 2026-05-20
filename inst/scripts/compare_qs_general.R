args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript compare_qs_general.R before.rds current.rds")
}

before <- readRDS(args[[1]])
current <- readRDS(args[[2]])

qs_before <- before$quality_score
qs_current <- current$quality_score

cat("BEFORE\n")
cat("Branch:", before$git_branch, "\n")
cat("Commit:", before$git_hash, "\n")
cat("QScore column:", before$qscore_column, "\n")
cat("Detection:", before$qscore_detection_method, "\n\n")

cat("CURRENT\n")
cat("Branch:", current$git_branch, "\n")
cat("Commit:", current$git_hash, "\n")
cat("QScore column:", current$qscore_column, "\n")
cat("Detection:", current$qscore_detection_method, "\n\n")

cat("Lengths:\n")
cat("before:", length(qs_before), "\n")
cat("current:", length(qs_current), "\n\n")

if (length(qs_before) != length(qs_current)) {
    stop("Different number of Quality Score values. Cannot compare directly.")
}

comparison <- data.frame(
    index=seq_along(qs_before),
    quality_score_before=qs_before,
    quality_score_current=qs_current,
    difference=qs_current - qs_before
)

cat("Summary before:\n")
print(summary(qs_before))

cat("\nSummary current:\n")
print(summary(qs_current))

cat("\nSummary difference current - before:\n")
print(summary(comparison$difference))

cat("\nall.equal:\n")
print(all.equal(qs_before, qs_current))

cat("\nExactly different values:\n")
print(sum(qs_before != qs_current, na.rm=TRUE))

cat("\nDifferent above tolerance 1e-8:\n")
print(sum(abs(qs_current - qs_before) > 1e-8, na.rm=TRUE))

cat("\nDifferent above tolerance 1e-6:\n")
print(sum(abs(qs_current - qs_before) > 1e-6, na.rm=TRUE))

cat("\nDifferent above tolerance 1e-4:\n")
print(sum(abs(qs_current - qs_before) > 1e-4, na.rm=TRUE))

out_csv <- file.path(
    dirname(args[[2]]),
    "quality_score_comparison.csv"
)

out_rds <- file.path(
    dirname(args[[2]]),
    "quality_score_comparison.rds"
)

write.csv(comparison, out_csv, row.names=FALSE)
saveRDS(comparison, out_rds)

cat("\nSaved:\n")
cat(out_csv, "\n")
cat(out_rds, "\n")
