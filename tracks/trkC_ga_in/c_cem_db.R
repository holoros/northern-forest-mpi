#!/usr/bin/env Rscript
# Multistate track (Amendment 1, A1.4): CEM FIA database fia_db_<ST>.rds for GA and IN. Port of tools/ncasi_rc/eng/cem_build_db.R
# (same six tables, same class and column handling) with the donor list taken from the engine rule (engines/cem/engine_seeded/
# run_projection.R get_donor_states): GA = GA, FL, SC, NC, TN, AL; IN has no entry, engine default = IN only.
# Writes multistate/restricted/cem_rds/fia_db_<ST>.rds (600). Usage: ms_cem_db.R ST
suppressPackageStartupMessages({library(data.table); library(tibble)})
ST <- commandArgs(TRUE)[1]; Sys.umask("077")
don <- list(GA = c("GA", "FL", "SC", "NC", "TN", "AL"), IN = c("IN"))[[ST]]
R <- path.expand("~/restricted/ncasi-modeleval/fia_data_rc"); OUT <- path.expand("~/perseus_me/tracks/trkC_gain/restricted/cem/rds")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE, mode = "0700")
dir_of <- function(s) if (s %in% c("GA", "IN")) file.path(R, s) else file.path(R, "donors", s)
tabs <- c("PLOT","COND","TREE","TREE_GRM_COMPONENT","TREE_GRM_MIDPT","SUBP_COND_CHNG_MTRX")
for (s in don) for (t in tabs) if (!file.exists(file.path(dir_of(s), sprintf("%s_%s.csv", s, t)))) stop("missing donor table ", s, "_", t)
db <- lapply(tabs, function(t) as.data.frame(rbindlist(lapply(don, function(s) fread(file.path(dir_of(s), sprintf("%s_%s.csv", s, t)), showProgress = FALSE, integer64 = "double")), fill = TRUE)))
names(db) <- tabs; class(db) <- "FIA.Database"
f <- file.path(OUT, sprintf("fia_db_%s.rds", ST)); saveRDS(db, f); Sys.chmod(f, "0600")
cat(ST, "donors", paste(don, collapse = ","), "\n"); for (t in tabs) cat(t, nrow(db[[t]]), "\n")
