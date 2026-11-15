## trkC_gain PREREG 2.2: run the registered NCASI builder (md5 af95202619a1353c5d9b31bd0a9fbc61, file unchanged) on the GA/IN FVS
## default growth points. Substitutions applied in memory with asserted counts: input file, output folder, state list of the figure,
## arm label calibrated -> default (figure, summary, sensitivity). usage: Rscript c_fit_sit.R BUILDER POINTS_CSV OUT_DIR
a <- commandArgs(TRUE); src <- readLines(a[1])
stopifnot(unname(tools::md5sum(a[1])) == "af95202619a1353c5d9b31bd0a9fbc61")
sub1 <- function(s, from, to, n) { k <- sum(grepl(from, s, fixed = TRUE)); if (k != n) stop(sprintf("expected %d lines with '%s', found %d", n, from, k)); gsub(from, to, s, fixed = TRUE) }
src <- sub1(src, 'IN_CSV    <- "fvs_growth_points_DATA_2026-09-17.csv"', sprintf('IN_CSV    <- "%s"', a[2]), 1)
src <- sub1(src, 'OUT_DIR   <- "output"', sprintf('OUT_DIR   <- "%s"', a[3]), 1)
src <- sub1(src, 'sts <- c("ME", "NH", "VT", "MI", "WI", "MN")', 'sts <- c("GA", "IN")', 1)
src <- sub1(src, 'arm == "calibrated"', 'arm == "default"', 4)
src <- sub1(src, 'mfrow = c(2, 3)', 'mfrow = c(1, 2)', 1)
if (length(a) >= 4 && a[4] == "weighted") {   # trkC_gain stress variant cpw: nlsLM weights = column w (CONDPROP_UNADJ)
  src <- sub1(src, 'fit_cr <- function(age, vol, b0_mult = B0_UPPER_MULT) {', 'fit_cr <- function(age, vol, b0_mult = B0_UPPER_MULT, w = NULL) {', 1)
  src <- sub1(src, 'nlsLM(vol ~ b0 * (1 - exp(-b1 * age))^b2, start = st,', 'nlsLM(vol ~ b0 * (1 - exp(-b1 * age))^b2, start = st, weights = if (is.null(w)) rep(1, length(vol)) else w,', 1)
  src <- sub1(src, 'm <- fit_cr(d$stand_age, d$vol)', 'm <- fit_cr(d$stand_age, d$vol, w = if ("w" %in% names(d)) d$w else NULL)', 1)
}
eval(parse(text = src))
