# F15 second path, base R + jsonlite: point S of the engine-refined rows (C2 by domain on P*, C5 net statewide) from the
# exported pair values, written independently of e_score.py / imf2_lib.comp_S. Rscript second.R RESTRICTED_TAG_DIR SCORES_CSV OUT_CSV
a <- commandArgs(TRUE); X <- read.csv(file.path(a[1], "second_in.csv"), colClasses = c(pair_id = "character", dom = "character"))
D <- jsonlite::fromJSON(file.path(a[1], "second_doms.json")); S <- read.csv(a[2]); MINU <- 30; FL <- 0.01
X$cov <- as.logical(X$cov)
lnr <- function(m, o) ifelse(o > 0, log(pmax(m, FL * o) / o), NA)
mz <- function(v, c) { v[is.na(v)] <- 0; sum(v[c]) / sum(c) }
out <- NULL
for (md in unique(X$version)) {
  c2 <- X[X$version == md & X$set == "c2", ]; c5 <- X[X$version == md & X$set == "c5", ]
  e <- sapply(c("g", "m", "n"), function(q) {
    v <- sapply(seq_along(D$doms), function(i) { cc <- c2$cov & c2$dom == D$doms[i]
      if (sum(cc) < MINU) return(NA); abs(lnr(mz(c2[[paste0("pred_", q)]], cc), mz(c2[[paste0("obs_", q)]], cc))) })
    ok <- is.finite(v); sum(v[ok] * D$aw[ok]) / sum(D$aw[ok]) })
  s2 <- mean(e); s5 <- abs(lnr(mz(c5$pred_n, c5$cov), mz(c5$obs_n, c5$cov)))
  r <- S[S$version == md, ]
  out <- rbind(out, data.frame(version = md, S_C2_R = s2, S_C2_py = r$S_C2, d_C2 = s2 - r$S_C2, S_C5_R = s5, S_C5_py = r$S_C5net, d_C5 = s5 - r$S_C5net))
}
write.csv(out, a[3], row.names = FALSE); print(out, digits = 15)
