# Track R second code path (base R): MPI, coverage, tier, grade from a rows table's pillar points; compares with the Python values.
# usage: Rscript second_path_mpi.R ROWS_CSV OUT_JSON
a <- commandArgs(TRUE); x <- read.csv(a[1], stringsAsFactors = FALSE)
w <- c(P1 = 30, P2 = 20, P3 = 20, P4 = 15, P5 = 15)
mp <- cv <- rep(NA_real_, nrow(x))
for (i in seq_len(nrow(x))) {
  p <- unlist(x[i, names(w)]); ok <- is.finite(p)
  if (any(ok)) { mp[i] <- sum(w[ok] * p[ok]) / sum(w[ok]); cv[i] <- sum(w[ok]) / sum(w) } else { cv[i] <- 0 }
}
gr <- ifelse(!is.finite(mp), "", ifelse(mp >= 85, "A", ifelse(mp >= 70, "B", ifelse(mp >= 55, "C", ifelse(mp >= 40, "D", "F")))))
ti <- ifelse(cv >= 0.65 - 1e-12, 1L, 2L)
d <- list(rows = nrow(x), max_abs_MPI = max(abs(mp - x$MPI), na.rm = TRUE), max_abs_cov = max(abs(cv - x$coverage), na.rm = TRUE),
          tier_equal = all(ti == x$tier), grade_equal = all(gr == ifelse(is.na(x$grade), "", x$grade)))
writeLines(sprintf('{"rows": %d, "max_abs_MPI": %.3g, "max_abs_cov": %.3g, "tier_equal": %s, "grade_equal": %s}', d$rows, d$max_abs_MPI, d$max_abs_cov,
                   tolower(d$tier_equal), tolower(d$grade_equal)), a[2])
cat(readLines(a[2]), "\n")
